#!/usr/bin/env python3
"""Exploratory, single-sample VEP-annotated WES variant filtering pipeline.

Not a clinical classifier. Keeps an auditable record of every parsed ALT allele.
"""
import argparse
import csv
import math
import json
import sys
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

KEEP_IMPACTS = {"HIGH", "MODERATE"}
OUTPUT_FIELDS = ["chrom", "pos", "ref", "alt", "genotype", "gene", "transcript",
                 "selection", "consequence", "IMPACT", "gnomAD_AF", "frequency_status",
                 "CADD", "REVEL", "SIFT", "PolyPhen", "ClinVar", "score", "evidence", "status", "reason",
                 "alt_dosage", "inheritance_compatibility", "all_csq_json"]
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
    alt_dosage: int | None = None
    inheritance_compatibility: str = ""

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
                required = {'Allele', 'IMPACT', 'Consequence'}
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
                    selection = 'canonical' if selected else 'all_transcripts'
                v = Variant(chrom, int(pos), ref, alt, q, filt, int(dp),
                            sample_data.get('GT', './.'), int(gq), selected or {}, own,
                            selection, alt_index=i+1)
                if uncertain[i]:
                    v.status, v.reason = 'review', 'Ambiguous/unmatched CSQ allele mapping; inspect original VCF'
                elif not own:
                    v.status, v.reason = 'review', 'No CSQ annotation for this ALT'
                elif selection == 'manual_not_found':
                    v.reason = f'Requested transcript {requested} not found; all transcripts remain in analysis'
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
    """Filter ALT alleles, not transcripts. No canonical annotation is required."""
    for v in variants:
        reason = qc_filter(v, args)
        if reason:
            v.status, v.reason = 'excluded_QC', reason
            continue
        # Inheritance is metadata, NEVER a reason to discard an allele.
        copies = allele_genotype(v)
        v.alt_dosage = copies
        if copies is None:
            v.inheritance_compatibility = 'unknown_genotype'
        elif args.model == 'dominant':
            v.inheritance_compatibility = 'compatible_simple_dominant' if copies == 1 else 'other_dosage_review'
        else:
            v.inheritance_compatibility = 'compatible_simple_homozygous_recessive' if copies == 2 else 'other_dosage_review_compound_het_not_assessed'
        if v.status == 'review':  # Only uncertain CSQ-to-ALT mapping prevents safe filtering.
            continue
        v.af, v.frequency_status = allele_frequency(v)
        if v.frequency_status in ('INVALID_AF', 'INCOMPLETE_AF', 'CONFLICTING_AF'):
            v.status, v.reason = 'review', f'Frequency needs review: {v.frequency_status}'
            continue
        # Independent, mutually exclusive frequency lanes; genotype never excludes.
        # Unknown AF remains visible in the main lane, explicitly flagged as unknown.
        if v.af is not None and v.af >= args.maf_recessive:
            v.status, v.reason = 'excluded_frequency', f'gnomAD_AF {v.af} >= {args.maf_recessive}'
            continue
        impacts = {c.get('IMPACT', '') for c in v.all_csq}
        if not (impacts & KEEP_IMPACTS):
            v.status, v.reason = 'excluded_consequence', f'No HIGH/MODERATE annotation: {sorted(impacts)}'
            continue
        v.status = ('additional_candidate' if v.af is not None and v.af >= args.maf_dominant
                    else 'candidate')
        notes = []
        if v.transcript_status == 'all_transcripts':
            notes.append('No canonical annotation: assessed all ALT-matched transcripts')
        if v.transcript_status == 'manual_not_found':
            notes.append('Requested transcript not found; assessed all ALT-matched transcripts')
        if v.csq and v.csq.get('IMPACT') not in KEEP_IMPACTS:
            notes.append('Selected/canonical transcript LOW or MODIFIER; another transcript is HIGH/MODERATE')
        if len(impacts & KEEP_IMPACTS) > 1 or ('HIGH' in impacts and any(x not in KEEP_IMPACTS for x in impacts)):
            notes.append('Transcript impact discordance: inspect all annotations')
        v.reason = '; '.join(notes)
    return [v for v in variants if v.status in ('candidate', 'additional_candidate')]


def prioritize(variants, panel):
    """Score once per ALT allele. Each evidence category contributes at most once."""
    for v in variants:
        annotations = v.all_csq
        score, evidence = 0, []
        genes = sorted({c.get('SYMBOL', '') for c in annotations if c.get('SYMBOL') in panel})
        if genes:
            score += 2
            evidence.append(f'in HPO candidate-gene panel ({", ".join(genes)})')
        if v.af is None:
            evidence.append(f'gnomAD AF unavailable ({v.frequency_status})')
        elif v.af == 0:
            score += 2
            evidence.append('gnomAD AF = 0 (reported value)')
        elif v.af < 1e-4:
            score += 1
            evidence.append(f'ultra-rare (gnomAD AF {v.af})')
        impacts = {c.get('IMPACT', '') for c in annotations}
        if 'HIGH' in impacts:
            score += 2
            evidence.append('HIGH consequence in at least one transcript')
        elif 'MODERATE' in impacts:
            score += 1
            evidence.append('MODERATE consequence in at least one transcript')
        revel = [number(c.get('REVEL')) for c in annotations]
        if any(x is not None and x >= .7 for x in revel):
            score += 1
            evidence.append('REVEL >= 0.7 in at least one annotation')
        cadd = [number(c.get('CADD_PHRED')) for c in annotations]
        if any(x is not None and x >= 20 for x in cadd):
            score += 1
            evidence.append('CADD >= 20 in at least one annotation')
        if any(c.get('SIFT', '').split('(')[0] == 'deleterious' for c in annotations):
            score += 1
            evidence.append('SIFT: deleterious in at least one annotation')
        if any(c.get('PolyPhen', '').split('(')[0] in ('probably_damaging', 'possibly_damaging') for c in annotations):
            score += 1
            evidence.append('PolyPhen damaging in at least one annotation')
        clinvars = {(c.get('ClinVar_CLNSIG') or '').lower() for c in annotations}
        if 'pathogenic' in clinvars or 'pathogenic/likely_pathogenic' in clinvars:
            score += 2
            evidence.append('ClinVar: pathogenic annotation (verify clinical context)')
        elif 'likely_pathogenic' in clinvars:
            score += 1
            evidence.append('ClinVar: likely pathogenic annotation (verify clinical context)')
        elif any(clinvars):
            evidence.append('ClinVar annotations require review; not automatically scored')
        v.score, v.evidence = score, evidence
    return sorted(variants, key=lambda v: (-v.score, v.chrom, v.pos, v.alt))


def candidate_row(v, rank=''):
    annotations = v.all_csq
    def unique(field):
        return '; '.join(sorted({c.get(field, '') for c in annotations if c.get(field, '')}))
    return dict(rank=rank, chrom=v.chrom, pos=v.pos, ref=v.ref, alt=v.alt,
                genotype=v.gt, gene=unique('SYMBOL'), transcript=unique('Feature'),
                selection=v.transcript_status, consequence=unique('Consequence'),
                IMPACT=unique('IMPACT'), gnomAD_AF='' if v.af is None else v.af,
                frequency_status=v.frequency_status, CADD=unique('CADD_PHRED'),
                REVEL=unique('REVEL'), SIFT=unique('SIFT'), PolyPhen=unique('PolyPhen'),
                ClinVar=unique('ClinVar_CLNSIG'), score='' if v.score is None else v.score,
                evidence='; '.join(v.evidence), status=v.status, reason=v.reason,
                alt_dosage='' if v.alt_dosage is None else v.alt_dosage,
                inheritance_compatibility=v.inheritance_compatibility,
                all_csq_json=json.dumps(annotations, ensure_ascii=False, separators=(',', ':')))


def write_table(path, columns, rows):
    with open(path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def write_outputs(variants, ranked, args):
    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    primary = [v for v in ranked if v.status == 'candidate']
    additional = [v for v in ranked if v.status == 'additional_candidate']
    write_table(str(prefix) + '.tsv', ['rank'] + OUTPUT_FIELDS,
                [candidate_row(v, i) for i, v in enumerate(primary, 1)])
    write_table(str(prefix) + '.additional_candidates.tsv', ['rank'] + OUTPUT_FIELDS,
                [candidate_row(v, i) for i, v in enumerate(additional, 1)])
    # Separate audit output: excluded/review alleles do not clutter candidate TSV.
    excluded = []
    for v in variants:
        if v.status not in ('candidate', 'additional_candidate'):
            row = candidate_row(v)
            row['reason'] = v.reason
            excluded.append(row)
    write_table(str(prefix) + '.excluded.tsv', ['rank'] + OUTPUT_FIELDS, excluded)
    counts = Counter(v.status for v in variants)
    with open(str(prefix) + '.report.txt', 'w') as fh:
        fh.write('EXPLORATORY WES VARIANT FILTERING REPORT\n' + '=' * 48 + '\n')
        fh.write(f'Input: {args.vcf}\nPanel: {args.panel}\nModel: {args.model}\n')
        fh.write(f'Parallel AF lanes: primary AF < {args.maf_dominant}, additional {args.maf_dominant} <= AF < {args.maf_recessive}; '
                 f'QUAL >= {args.min_qual}, DP >= {args.min_dp}, GQ >= {args.min_gq}; FILTER=PASS\n')
        fh.write(f'ALT alleles processed: {len(variants)}\n'
                 f'Primary candidates: {len(primary)}\nAdditional candidates: {len(additional)}\n'
                 f'Excluded or review: {len(variants)-len(primary)-len(additional)}\n')
        for key, value in sorted(counts.items()):
            fh.write(f'  {key}: {value}\n')
        fh.write('\nCANDIDATES (one score per ALT, using all matched transcripts)\n')
        for i, v in enumerate(primary, 1):
            fh.write(f'[{i}] {v.chrom}:{v.pos} {v.ref}>{v.alt} '
                     f'{v.csq.get("SYMBOL", "")} {v.csq.get("Feature", "")} '
                     f'[{v.csq.get("IMPACT", "")}] score={v.score}\n')
            for item in v.evidence:
                fh.write(f'  - {item}\n')
            if v.reason:
                fh.write(f'  REVIEW FLAG: {v.reason}\n')
        fh.write('\nADDITIONAL CANDIDATES (strict AF failed; extended AF passed)\n')
        for i, v in enumerate(additional, 1):
            fh.write(f'[{i}] {v.chrom}:{v.pos} {v.ref}>{v.alt} AF={v.af} score={v.score}\n')
            for item in v.evidence:
                fh.write(f'  - {item}\n')
            if v.reason:
                fh.write(f'  REVIEW FLAG: {v.reason}\n')
        fh.write('\nREVIEW QUEUE (full details in excluded.tsv)\n')
        for v in variants:
            if v.status == 'review':
                fh.write(f'{v.chrom}:{v.pos} {v.ref}>{v.alt}: {v.reason}\n')
                fh.write('  Available ENSTs: ' + ', '.join(c.get('Feature', '') for c in v.all_csq) + '\n')
        fh.write('\nEXCLUSIONS (auditable in excluded.tsv)\n')
        for v in variants:
            if v.status.startswith('excluded_'):
                fh.write(f'{v.chrom}:{v.pos} {v.ref}>{v.alt}: {v.status}: {v.reason}\n')
        fh.write('\nLIMITATIONS: Research/exploratory tool, not a clinical classification. '
                 'No tissue-expression or gene-disease mechanism inference. '
                 'Inheritance compatibility is descriptive, never an exclusion; no compound heterozygosity, '
                 'X-linked handling, segregation or CNV analysis. '
                 'Unmatched/ambiguous CSQ allele mapping requires manual inspection. '
                 'Canonical is informational and never determines filtering or scoring. '
                 'Unknown AF stays in the primary file flagged UNKNOWN (not assumed rare). '
                 'All ALT-matched CSQ annotations are retained per allele as JSON in the corresponding TSV.\n')
    return counts


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--vcf', required=True)
    parser.add_argument('--panel', required=True)
    parser.add_argument('--out-prefix', default='results/candidates')
    parser.add_argument('--model', choices=['dominant', 'recessive'], default='dominant')
    parser.add_argument('--transcript-override', help='TSV: chrom pos ref alt transcript (tab-separated)')
    parser.add_argument('--maf-dominant', type=float, default=1e-5)
    parser.add_argument('--maf-recessive', type=float, default=.05, help='Upper AF boundary of additional lane (exclusive)')
    parser.add_argument('--min-qual', type=float, default=30)
    parser.add_argument('--min-dp', type=int, default=20)
    parser.add_argument('--min-gq', type=int, default=20)
    args = parser.parse_args(argv)
    if not (0 <= args.maf_dominant < args.maf_recessive <= 1):
        parser.error('Require 0 <= --maf-dominant < --maf-recessive <= 1')
    overrides = load_overrides(args.transcript_override)
    variants = load_vcf(args.vcf, overrides)
    panel = load_panel(args.panel)
    ranked = prioritize(apply_filters(variants, args), panel)
    counts = write_outputs(variants, ranked, args)
    print(f'Processed {len(variants)} ALT alleles: {counts.get("candidate", 0)} primary, '
          f'{counts.get("additional_candidate", 0)} additional, '
          f'{counts.get("review", 0)} review; wrote {args.out_prefix}.*')
    return variants


if __name__ == '__main__':
    main()
