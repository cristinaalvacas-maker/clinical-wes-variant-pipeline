#!/usr/bin/env python3
"""Exploratory, single-sample VEP-annotated WES variant filtering pipeline.

Not a clinical classifier. Keeps an auditable record of every parsed ALT allele.
"""
import argparse
import csv
import math
import sys
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

KEEP_IMPACTS = {"HIGH", "MODERATE"}
OUTPUT_FIELDS = ["chrom", "pos", "ref", "alt", "genotype", "gene", "transcript",
                 "selection", "consequence", "IMPACT", "gnomAD_AF", "frequency_status",
                 "CADD", "REVEL", "SIFT", "PolyPhen", "ClinVar", "score", "evidence", "status", "reason"]
REVIEW_FIELDS = ["chrom", "pos", "ref", "alt", "genotype", "variant_status", "reason",
                 "transcript", "gene", "CANONICAL", "MANE_SELECT", "Consequence", "IMPACT",
                 "gnomAD_AF", "selected", "selection_status"]

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
    csq: dict = field(default_factory=dict)
    all_csq: list = field(default_factory=list)
    transcript_status: str = ""
    status: str = "unprocessed"
    reason: str = ""
    af: float | None = None
    frequency_status: str = ""
    score: int | None = None
    evidence: list = field(default_factory=list)
    alt_index: int = 1

    @property
    def key(self):
        return (self.chrom, str(self.pos), self.ref, self.alt)


def parse_info(raw):
    result = {}
    for item in raw.split(';'):
        if '=' in item:
            k, val = item.split('=', 1)
            result[k] = val
        elif item:
            result[item] = True
    return result


def parse_csq(raw, fields):
    annotations = []
    for annotation in raw.split(','):
        values = annotation.split('|')
        values += [''] * max(0, len(fields) - len(values))
        annotations.append(dict(zip(fields, values)))
    return annotations


def vep_allele(ref, alt):
    """VEP's commonly used minimal allele representation for anchored indels."""
    r, a = ref, alt
    while r and a and r[0] == a[0]:
        r, a = r[1:], a[1:]
    return a or '-'


def match_annotations(annotations, ref, alts):
    """Match CSQ Allele to ALT. Ambiguous matches are flagged, never guessed."""
    matches = [[] for _ in alts]
    uncertain = [False] * len(alts)
    representations = [{alt, vep_allele(ref, alt)} for alt in alts]
    for c in annotations:
        allele = c.get('Allele', '')
        indices = [i for i, reps in enumerate(representations) if allele and allele in reps]
        if len(indices) == 1:
            matches[indices[0]].append(c)
        else:
            # An ambiguous/unmatched annotation cannot safely be attributed to an ALT.
            for i in (indices if indices else range(len(alts))):
                uncertain[i] = True
    return matches, uncertain


def number(value, default=None):
    try:
        result = float(value)
        return result if math.isfinite(result) else default
    except (ValueError, TypeError):
        return default


def load_overrides(path):
    if not path:
        return {}
    overrides = {}
    with open(path, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        needed = {'chrom', 'pos', 'ref', 'alt', 'transcript'}
        if not needed.issubset(reader.fieldnames or []):
            raise ValueError('Override TSV needs columns: chrom, pos, ref, alt, transcript')
        for row in reader:
            key = tuple(row[k].strip() for k in ('chrom', 'pos', 'ref', 'alt'))
            value = row['transcript'].strip()
            if not all(key) or not value or key in overrides:
                raise ValueError(f'Invalid or duplicate transcript override: {key}')
            overrides[key] = value
    return overrides


def load_vcf(path, overrides=None):
    overrides = overrides or {}
    variants = []
    csq_fields = None
    seen_overrides = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith('##INFO=<ID=CSQ,'):
                if 'Format: ' not in line:
                    raise ValueError('CSQ header found without Format definition')
                csq_fields = line.split('Format: ', 1)[1].split('"', 1)[0].rstrip('>\n').split('|')
                required = {'Allele', 'Feature', 'CANONICAL', 'IMPACT', 'Consequence', 'gnomAD_AF'}
                missing = required - set(csq_fields)
                if missing:
                    raise ValueError(f'CSQ Format missing required fields: {sorted(missing)}')
                if len(set(csq_fields)) != len(csq_fields):
                    raise ValueError('Duplicate CSQ field names in VCF header')
                continue
            if line.startswith('#'):
                continue
            if not line.strip():
                continue
            if csq_fields is None:
                raise ValueError('CSQ Format definition not found in VCF header')
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 10:
                raise ValueError('Expected single-sample VCF with at least 10 columns')
            chrom, pos, _, ref, alt_raw, qual, filt, info_raw, fmt, sample = cols[:10]
            alts = alt_raw.split(',')
            info = parse_info(info_raw)
            sample_data = dict(zip(fmt.split(':'), sample.split(':')))
            dp = number(sample_data.get('DP', info.get('DP')), 0)
            gq = number(sample_data.get('GQ'), 0)
            q = number(qual, 0)
            all_annotations = parse_csq(info['CSQ'], csq_fields) if isinstance(info.get('CSQ'), str) else []
            grouped, uncertain = match_annotations(all_annotations, ref, alts)
            for i, alt in enumerate(alts):
                key = (chrom, str(pos), ref, alt)
                own = grouped[i]
                requested = overrides.get(key)
                selected = None
                if requested:
                    seen_overrides.add(key)
                    selected = next((c for c in own if c.get('Feature') == requested), None)
                    selection = 'manual' if selected else 'manual_not_found'
                else:
                    selected = next((c for c in own if c.get('CANONICAL') == 'YES'), None)
                    selection = 'canonical' if selected else 'review'
                v = Variant(chrom, int(pos), ref, alt, q, filt, int(dp),
                            sample_data.get('GT', './.'), int(gq), selected or {}, own,
                            selection, alt_index=i+1)
                if uncertain[i]:
                    v.status, v.reason = 'review', 'Ambiguous/unmatched CSQ allele mapping; inspect original VCF'
                elif not own:
                    v.status, v.reason = 'review', 'No CSQ annotation for this ALT'
                elif selection == 'manual_not_found':
                    v.status, v.reason = 'review', f'Requested transcript {requested} not found for this ALT'
                elif selection == 'review':
                    v.status, v.reason = 'review', 'No canonical transcript; select an ENST in overrides TSV'
                variants.append(v)
    if csq_fields is None:
        raise ValueError('CSQ Format definition not found in VCF header')
    unused = set(overrides) - seen_overrides
    if unused:
        raise ValueError(f'Overrides not matched to VCF variants: {sorted(unused)}')
    return variants


def load_panel(path):
    with open(path) as fh:
        return {line.strip() for line in fh if line.strip() and not line.startswith('#')}


def qc_filter(v, args):
    if v.filt != 'PASS':
        return f'FILTER={v.filt} (only PASS accepted)'
    if v.qual < args.min_qual:
        return f'QUAL {v.qual} < {args.min_qual}'
    if v.dp < args.min_dp:
        return f'DP {v.dp} < {args.min_dp}'
    if v.gq < args.min_gq:
        return f'GQ {v.gq} < {args.min_gq}'
    return ''


def allele_genotype(v):
    """Normalize the selected ALT allele as 0/0, 0/1 or 1/1."""
    raw = v.gt.replace('|', '/').split('/')
    if len(raw) != 2 or '.' in raw:
        return None
    try:
        alleles = [int(x) for x in raw]
    except ValueError:
        return None
    copies = alleles.count(v.alt_index)
    return copies


def allele_frequency(v):
    """Require concordant known AF across ALT-matched annotations; never use arbitrary transcript."""
    raw = [c.get('gnomAD_AF', '') for c in v.all_csq]
    known = [number(x) for x in raw if x not in ('', '.', None)]
    if any(x is None or x < 0 or x > 1 for x in known):
        return None, 'INVALID_AF'
    if not known:
        return None, 'UNKNOWN'
    if len(known) != len(raw):
        return None, 'INCOMPLETE_AF'
    if max(known) - min(known) > 1e-12:
        return None, 'CONFLICTING_AF'
    return known[0], 'KNOWN'


def apply_filters(variants, args):
    for v in variants:
        # QC and allele-specific inheritance apply even to annotations requiring review.
        reason = qc_filter(v, args)
        if reason:
            v.status, v.reason = 'excluded_QC', reason
            continue
        copies = allele_genotype(v)
        if copies is None:
            v.status, v.reason = 'review', 'Missing/unsupported genotype; inspect ALT dosage'
            continue
        if args.model == 'dominant' and copies != 1:
            v.status, v.reason = 'excluded_inheritance', f'ALT dosage {copies}, expected 1 (simple dominant model)'
            continue
        if args.model == 'recessive' and copies != 2:
            v.status, v.reason = 'excluded_inheritance', f'ALT dosage {copies}, expected 2 (simple recessive model)'
            continue
        if v.status == 'review':
            continue
        v.af, v.frequency_status = allele_frequency(v)
        if v.frequency_status in ('INVALID_AF', 'INCOMPLETE_AF', 'CONFLICTING_AF'):
            v.status, v.reason = 'review', f'Frequency needs review: {v.frequency_status}'
            continue
        maf = args.maf_recessive if args.model == 'recessive' else args.maf_dominant
        if v.af is not None and v.af >= maf:
            v.status, v.reason = 'excluded_frequency', f'gnomAD_AF {v.af} >= {maf}'
            continue
        impacts = {c.get('IMPACT', '') for c in v.all_csq}
        if not (impacts & KEEP_IMPACTS):
            v.status, v.reason = 'excluded_consequence', f'No HIGH/MODERATE annotation: {sorted(impacts)}'
            continue
        selected_impact = v.csq.get('IMPACT', '')
        if selected_impact not in KEEP_IMPACTS:
            v.status, v.reason = 'review', 'Selected transcript is not HIGH/MODERATE; alternative transcript is'
            continue
        v.status = 'candidate'
        v.reason = ('Alternative transcript has HIGH impact; inspect alongside selected annotation'
                    if selected_impact != 'HIGH' and 'HIGH' in impacts else '')
    return [v for v in variants if v.status == 'candidate']


def prioritize(variants, panel):
    for v in variants:
        c = v.csq
        score, evidence = 0, []
        gene = c.get('SYMBOL', '')
        if gene in panel:
            score += 2
            evidence.append(f'in HPO candidate-gene panel ({gene})')
        if v.af is None:
            evidence.append('gnomAD AF unavailable (UNKNOWN)')
        elif v.af == 0:
            score += 2
            evidence.append('gnomAD AF = 0 (reported value)')
        elif v.af < 1e-4:
            score += 1
            evidence.append(f'ultra-rare (gnomAD AF {v.af})')
        impact = c.get('IMPACT', '')
        if impact == 'HIGH':
            score += 2
            evidence.append(f'high-impact consequence ({c.get("Consequence")})')
        elif impact == 'MODERATE':
            score += 1
            evidence.append(f'moderate-impact consequence ({c.get("Consequence")})')
        revel, cadd = number(c.get('REVEL'), 0), number(c.get('CADD_PHRED'), 0)
        if revel >= .7:
            score += 1
            evidence.append(f'REVEL {revel} (supporting computational evidence)')
        if cadd >= 20:
            score += 1
            evidence.append(f'CADD {cadd} (supporting computational evidence)')
        if c.get('SIFT', '').split('(')[0] == 'deleterious':
            score += 1
            evidence.append('SIFT: deleterious')
        if c.get('PolyPhen', '').split('(')[0] in ('probably_damaging', 'possibly_damaging'):
            score += 1
            evidence.append(f'PolyPhen: {c.get("PolyPhen")}')
        clinvar = (c.get('ClinVar_CLNSIG') or '').lower()
        if clinvar in ('pathogenic', 'pathogenic/likely_pathogenic'):
            score += 2
            evidence.append(f'ClinVar: {clinvar}')
        elif clinvar == 'likely_pathogenic':
            score += 1
            evidence.append('ClinVar: likely pathogenic')
        elif clinvar:
            evidence.append(f'ClinVar: {clinvar} (not automatically scored)')
        v.score, v.evidence = score, evidence
    return sorted(variants, key=lambda v: (-v.score, v.chrom, v.pos, v.alt))


def candidate_row(v, rank=''):
    c = v.csq
    return dict(rank=rank, chrom=v.chrom, pos=v.pos, ref=v.ref, alt=v.alt,
                genotype=v.gt, gene=c.get('SYMBOL', ''), transcript=c.get('Feature', ''),
                selection=v.transcript_status, consequence=c.get('Consequence', ''),
                IMPACT=c.get('IMPACT', ''), gnomAD_AF='' if v.af is None else v.af,
                frequency_status=v.frequency_status, CADD=c.get('CADD_PHRED', ''),
                REVEL=c.get('REVEL', ''), SIFT=c.get('SIFT', ''),
                PolyPhen=c.get('PolyPhen', ''), ClinVar=c.get('ClinVar_CLNSIG', ''),
                score='' if v.score is None else v.score, evidence='; '.join(v.evidence),
                status=v.status, reason=v.reason)


def write_table(path, columns, rows):
    with open(path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def write_outputs(variants, ranked, args):
    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    write_table(str(prefix) + '.tsv', ['rank'] + OUTPUT_FIELDS,
                [candidate_row(v, i) for i, v in enumerate(ranked, 1)])
    # Every parsed ALT is included here, including excluded variants and all alternative transcripts.
    review_rows = []
    for v in variants:
        annotations = v.all_csq or [{}]
        for c in annotations:
            review_rows.append(dict(chrom=v.chrom, pos=v.pos, ref=v.ref, alt=v.alt,
                                    genotype=v.gt, variant_status=v.status, reason=v.reason,
                                    transcript=c.get('Feature', ''), gene=c.get('SYMBOL', ''),
                                    CANONICAL=c.get('CANONICAL', ''), MANE_SELECT=c.get('MANE_SELECT', ''),
                                    Consequence=c.get('Consequence', ''), IMPACT=c.get('IMPACT', ''),
                                    gnomAD_AF=c.get('gnomAD_AF', ''),
                                    selected='YES' if c is v.csq else 'NO',
                                    selection_status=v.transcript_status))
    write_table(str(prefix) + '.transcript_review.tsv', REVIEW_FIELDS, review_rows)
    counts = Counter(v.status for v in variants)
    with open(str(prefix) + '.report.txt', 'w') as fh:
        fh.write('EXPLORATORY WES VARIANT FILTERING REPORT\n' + '=' * 48 + '\n')
        fh.write(f'Input: {args.vcf}\nPanel: {args.panel}\nModel: {args.model}\n')
        fh.write(f'Cutoffs: dominant AF < {args.maf_dominant}, recessive AF < {args.maf_recessive}; '
                 f'QUAL >= {args.min_qual}, DP >= {args.min_dp}, GQ >= {args.min_gq}; FILTER=PASS\n')
        fh.write(f'ALT alleles processed: {len(variants)}\nCandidates: {len(ranked)}\n')
        for key, value in sorted(counts.items()):
            fh.write(f'  {key}: {value}\n')
        fh.write('\nCANDIDATES (selected-transcript score only)\n')
        for i, v in enumerate(ranked, 1):
            fh.write(f'[{i}] {v.chrom}:{v.pos} {v.ref}>{v.alt} '
                     f'{v.csq.get("SYMBOL", "")} {v.csq.get("Feature", "")} '
                     f'[{v.csq.get("IMPACT", "")}] score={v.score}\n')
            for item in v.evidence:
                fh.write(f'  - {item}\n')
            if v.reason:
                fh.write(f'  REVIEW FLAG: {v.reason}\n')
        fh.write('\nREVIEW QUEUE\n')
        for v in variants:
            if v.status == 'review':
                fh.write(f'{v.chrom}:{v.pos} {v.ref}>{v.alt}: {v.reason}\n')
                fh.write('  Available ENSTs: ' + ', '.join(c.get('Feature', '') for c in v.all_csq) + '\n')
        fh.write('\nEXCLUSIONS (auditable in transcript_review.tsv)\n')
        for v in variants:
            if v.status.startswith('excluded_'):
                fh.write(f'{v.chrom}:{v.pos} {v.ref}>{v.alt}: {v.status}: {v.reason}\n')
        fh.write('\nLIMITATIONS: Research/exploratory tool, not a clinical classification. '
                 'No tissue-expression or gene-disease mechanism inference. '
                 'Only diploid single-sample SNV/small-indel genotypes and simple heterozygous dominant / '
                 'homozygous recessive models; no compound heterozygosity, X-linked handling, segregation or CNV analysis. '
                 'Unmatched/ambiguous CSQ allele mapping requires manual inspection. '
                 'Canonical is a default, not a guarantee of clinical relevance. '
                 'All parsed ALT/transcript annotations, including excluded variants, are retained in the review TSV.\n')
    return counts


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--vcf', required=True)
    parser.add_argument('--panel', required=True)
    parser.add_argument('--out-prefix', default='results/candidates')
    parser.add_argument('--model', choices=['dominant', 'recessive'], default='dominant')
    parser.add_argument('--transcript-override', help='TSV: chrom pos ref alt transcript (tab-separated)')
    parser.add_argument('--maf-dominant', type=float, default=1e-5)
    parser.add_argument('--maf-recessive', type=float, default=.05)
    parser.add_argument('--min-qual', type=float, default=30)
    parser.add_argument('--min-dp', type=int, default=20)
    parser.add_argument('--min-gq', type=int, default=20)
    args = parser.parse_args(argv)
    overrides = load_overrides(args.transcript_override)
    variants = load_vcf(args.vcf, overrides)
    panel = load_panel(args.panel)
    ranked = prioritize(apply_filters(variants, args), panel)
    counts = write_outputs(variants, ranked, args)
    print(f'Processed {len(variants)} ALT alleles: {len(ranked)} candidates, '
          f'{counts.get("review", 0)} review; wrote {args.out_prefix}.*')
    return variants


if __name__ == '__main__':
    main()
