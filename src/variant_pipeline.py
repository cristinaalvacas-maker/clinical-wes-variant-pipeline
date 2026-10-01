#!/usr/bin/env python3
"""
variant_pipeline.py — Clinical WES variant filtering & prioritization pipeline.

Implements the filtering workflow described in a clinical WES protocol
(master's thesis, Precision Medicine & Clinical Genetics):
  1. Quality filtering (QUAL / depth / genotype quality)
  2. Population-frequency filtering (gnomAD MAF thresholds by inheritance model)
  3. Functional-consequence filtering (VEP IMPACT)
  4. Inheritance-model filtering (dominant / recessive)
  5. Phenotype-driven prioritization (HPO candidate-gene panel)
  6. Ranked candidate report with supporting evidence

Input:  VEP-annotated VCF (CSQ field in INFO), single sample.
Output: ranked candidate table (TSV) + human-readable summary report.

Pure Python 3, standard library only.
"""

import argparse
import csv
import sys
from dataclasses import dataclass, field

# VEP CSQ sub-fields (must match the ##INFO=<ID=CSQ,...> header of the input VCF)
CSQ_FIELDS = [
    "Allele", "Consequence", "IMPACT", "SYMBOL", "HGVSc", "HGVSp",
    "gnomAD_AF", "CADD_PHRED", "SIFT", "PolyPhen", "REVEL", "ClinVar_CLNSIG",
]

# Consequences we keep at the consequence-filtering step.
# (Synonymous variants are noted separately: they can affect splicing and
# should not be discarded blindly in a real clinical review.)
KEEP_IMPACTS = {"HIGH", "MODERATE"}


@dataclass
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str
    qual: float
    filt: str
    dp: int
    gt: str
    gq: int
    csq: dict
    dropped_by: str = ""          # filter name that removed it, "" if kept
    drop_reason: str = ""
    score: int = 0
    evidence: list = field(default_factory=list)


def parse_info(info_str):
    info = {}
    for item in info_str.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            info[k] = v
        else:
            info[item] = True
    return info


def parse_csq(csq_raw):
    """Parse the first CSQ annotation (most severe consequence, VEP default order)."""
    first = csq_raw.split(",")[0]
    parts = first.split("|")
    # pad in case the annotation has fewer fields
    parts += [""] * (len(CSQ_FIELDS) - len(parts))
    return dict(zip(CSQ_FIELDS, parts))


def load_vcf(path):
    variants = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 10:
                continue
            chrom, pos, _id, ref, alt, qual, filt, info_s, fmt, sample = cols[:10]
            info = parse_info(info_s)
            if "CSQ" not in info:
                print(f"warning: variant {chrom}:{pos} has no CSQ annotation, skipped",
                      file=sys.stderr)
                continue
            fmt_keys = fmt.split(":")
            fmt_vals = sample.split(":")
            s = dict(zip(fmt_keys, fmt_vals))
            try:
                dp = int(s.get("DP", info.get("DP", 0)) or 0)
            except ValueError:
                dp = 0
            try:
                gq = int(s.get("GQ", 0) or 0)
            except ValueError:
                gq = 0
            try:
                q = float(qual)
            except ValueError:
                q = 0.0
            variants.append(Variant(
                chrom=chrom, pos=int(pos), ref=ref, alt=alt,
                qual=q, filt=filt, dp=dp, gt=s.get("GT", "./."), gq=gq,
                csq=parse_csq(info["CSQ"]),
            ))
    return variants


def load_panel(path):
    with open(path) as fh:
        return {ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")}


def fnum(s, default=0.0):
    try:
        return float(s)
    except (ValueError, TypeError):
        return default


# ---------------------------------------------------------------- filters

def qc_filter(v, min_qual, min_dp, min_gq):
    if v.filt not in ("PASS", "."):
        return False, f"FILTER={v.filt}"
    if v.qual < min_qual:
        return False, f"QUAL {v.qual} < {min_qual}"
    if v.dp < min_dp:
        return False, f"DP {v.dp} < {min_dp}"
    if v.gq < min_gq:
        return False, f"GQ {v.gq} < {min_gq}"
    return True, ""


def frequency_filter(v, maf_cutoff):
    af = fnum(v.csq.get("gnomAD_AF"), None)
    if af is None:
        return True, "gnomAD_AF=UNKNOWN"
    if af >= maf_cutoff:
        return False, f"gnomAD_AF {af} >= {maf_cutoff}"
    return True, ""

def consequence_filter(v):
    impact = v.csq.get("IMPACT", "")
    if impact not in KEEP_IMPACTS:
        return False, f"IMPACT={impact or 'missing'} ({v.csq.get('Consequence','')})"
    return True, ""


def inheritance_filter(v, model):
    gt = v.gt.replace("|", "/")
    if model == "dominant":
        if gt in ("0/1", "1/0"):
            return True, ""
        return False, f"genotype {v.gt} not heterozygous (dominant model)"
    if model == "recessive":
        if gt in ("1/1",):
            return True, ""
        return False, f"genotype {v.gt} not homozygous (recessive model)"
    return True, ""


FILTERS = [
    ("QC", qc_filter),
    ("frequency", frequency_filter),
    ("consequence", consequence_filter),
    ("inheritance", inheritance_filter),
]


def apply_filters(variants, args):
    kept = []
    for v in variants:
        ok, reason = qc_filter(v, args.min_qual, args.min_dp, args.min_gq)
        if not ok:
            v.dropped_by, v.drop_reason = "QC", reason
            continue
        maf = args.maf_recessive if args.model == "recessive" else args.maf_dominant
        ok, reason = frequency_filter(v, maf)
        if not ok:
            v.dropped_by, v.drop_reason = "frequency", reason
            continue
        ok, reason = consequence_filter(v)
        if not ok:
            v.dropped_by, v.drop_reason = "consequence", reason
            continue
        ok, reason = inheritance_filter(v, args.model)
        if not ok:
            v.dropped_by, v.drop_reason = "inheritance", reason
            continue
        kept.append(v)
    return kept


# ------------------------------------------------------- prioritization

def prioritize(variants, panel):
    """Score candidates; higher = stronger candidate. Evidence is recorded
    per variant so the ranking is fully explainable (as in a clinical report)."""
    for v in variants:
        c = v.csq
        score, ev = 0, []

        gene = c.get("SYMBOL", "")
        if gene in panel:
            score += 2
            ev.append(f"in HPO candidate-gene panel ({gene})")

        af = fnum(c.get("gnomAD_AF"), None)
        if af is None:
            ev.append("gnomAD: frequency unavailable (UNKNOWN)")
        elif af == 0:
            score += 2
            ev.append("gnomAD AF = 0 (reported value)")
        elif af < 1e-4:
            score += 1
            ev.append(f"ultra-rare (gnomAD AF {af})")

        impact = c.get("IMPACT", "")
        if impact == "HIGH":
            score += 2
            ev.append(f"high-impact consequence ({c.get('Consequence')})")
        elif impact == "MODERATE":
            score += 1
            ev.append(f"moderate-impact consequence ({c.get('Consequence')})")

        revel = fnum(c.get("REVEL"), 0.0)
        cadd = fnum(c.get("CADD_PHRED"), 0.0)
        if revel >= 0.7:
            score += 1
            ev.append(f"REVEL {revel} (damaging)")
        if cadd >= 20:
            score += 1
            ev.append(f"CADD {cadd} (top 1% deleterious)")
        if c.get("SIFT") == "deleterious":
            score += 1
            ev.append("SIFT: deleterious")
        if c.get("PolyPhen") in ("probably_damaging", "possibly_damaging"):
            score += 1
            ev.append(f"PolyPhen: {c.get('PolyPhen')}")

        clinvar = (c.get("ClinVar_CLNSIG") or "").lower()
        if "pathogenic" in clinvar and "likely" not in clinvar and "benign" not in clinvar:
            score += 2
            ev.append("ClinVar: pathogenic")
        elif "likely_pathogenic" in clinvar:
            score += 1
            ev.append("ClinVar: likely pathogenic")
        elif "uncertain" in clinvar:
            ev.append("ClinVar: uncertain significance (VUS)")

        v.score, v.evidence = score, ev

    return sorted(variants, key=lambda v: (-v.score, v.chrom, v.pos))


# ------------------------------------------------------------- reporting

def write_tsv(variants, path):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["rank", "chrom", "pos", "ref", "alt", "gene", "HGVSc", "HGVSp",
                    "consequence", "IMPACT", "genotype", "gnomAD_AF",
                    "CADD", "REVEL", "SIFT", "PolyPhen", "ClinVar", "score"])
        for i, v in enumerate(variants, 1):
            c = v.csq
            w.writerow([i, v.chrom, v.pos, v.ref, v.alt, c.get("SYMBOL"),
                        c.get("HGVSc"), c.get("HGVSp"), c.get("Consequence"),
                        c.get("IMPACT"), v.gt, c.get("gnomAD_AF"),
                        c.get("CADD_PHRED"), c.get("REVEL"), c.get("SIFT"),
                        c.get("PolyPhen"), c.get("ClinVar_CLNSIG"), v.score])


def write_report(variants, dropped, args, path):
    total = len(variants) + len(dropped)
    with open(path, "w") as fh:
        fh.write("VARIANT FILTERING & PRIORITIZATION REPORT\n")
        fh.write("=" * 50 + "\n")
        fh.write(f"Input VCF:            {args.vcf}\n")
        fh.write(f"Candidate gene panel: {args.panel}\n")
        fh.write(f"Inheritance model:    {args.model}\n")
        fh.write(f"MAF cutoff:           {args.maf_dominant} (dominant) / "
                 f"{args.maf_recessive} (recessive)\n")
        fh.write(f"QC cutoffs:           QUAL>={args.min_qual} DP>={args.min_dp} "
                 f"GQ>={args.min_gq}\n\n")
        fh.write(f"Variants processed:   {total}\n")
        fh.write(f"Variants retained:    {len(variants)}\n")
        by_filter = {}
        for v in dropped:
            by_filter[v.dropped_by] = by_filter.get(v.dropped_by, 0) + 1
        fh.write("Removed by filter:    " +
                 ", ".join(f"{k} ({n})" for k, n in sorted(by_filter.items())) + "\n")
        fh.write("\nCANDIDATE VARIANTS (ranked)\n")
        fh.write("-" * 50 + "\n")
        for i, v in enumerate(variants, 1):
            c = v.csq
            fh.write(f"\n[{i}] {v.chrom}:{v.pos} {v.ref}>{v.alt} "
                     f"({c.get('SYMBOL')}, {c.get('HGVSc')}, {c.get('HGVSp')}) "
                     f"score={v.score}\n")
            fh.write(f"    Consequence: {c.get('Consequence')} [{c.get('IMPACT')}] | "
                     f"genotype {v.gt} | gnomAD AF {c.get('gnomAD_AF')}\n")
            fh.write("    Evidence:\n")
            for e in v.evidence:
                fh.write(f"      - {e}\n")
        fh.write("\nNOTE: computational prioritization supports, but does not replace, "
                 "expert clinical interpretation (ACMG/AMP).\n")


# ------------------------------------------------------------------ main

def main():
    ap = argparse.ArgumentParser(
        description="Clinical WES variant filtering & prioritization pipeline")
    ap.add_argument("--vcf", required=True, help="VEP-annotated VCF (CSQ in INFO)")
    ap.add_argument("--panel", required=True, help="candidate gene list (one per line)")
    ap.add_argument("--out-prefix", default="results/candidates",
                    help="output prefix for .tsv and .report.txt")
    ap.add_argument("--model", choices=["dominant", "recessive"], default="dominant")
    ap.add_argument("--maf-dominant", type=float, default=1e-5)
    ap.add_argument("--maf-recessive", type=float, default=0.05)
    ap.add_argument("--min-qual", type=float, default=30.0)
    ap.add_argument("--min-dp", type=int, default=20)
    ap.add_argument("--min-gq", type=int, default=20)
    args = ap.parse_args()

    variants = load_vcf(args.vcf)
    panel = load_panel(args.panel)
    kept = apply_filters(variants, args)
    dropped = [v for v in variants if v.dropped_by]
    ranked = prioritize(kept, panel)

    write_tsv(ranked, args.out_prefix + ".tsv")
    write_report(ranked, dropped, args, args.out_prefix + ".report.txt")

    print(f"processed {len(variants)} variants → "
          f"{len(ranked)} candidates "
          f"({len(dropped)} filtered out)")
    print(f"wrote {args.out_prefix}.tsv and {args.out_prefix}.report.txt")
    if ranked:
        top = ranked[0].csq
        print(f"top candidate: {ranked[0].chrom}:{ranked[0].pos} "
              f"{top.get('SYMBOL')} {top.get('HGVSp')} (score {ranked[0].score})")


if __name__ == "__main__":
    main()
