# Clinical WES Variant Filtering & Prioritization Pipeline

A lightweight, dependency-free Python pipeline that implements a clinical
whole-exome sequencing (WES) variant filtering workflow: quality filtering,
population-frequency filtering, functional-consequence filtering,
inheritance-model filtering, and phenotype-driven prioritization with an
explainable evidence trail per candidate variant.

Built as a portfolio project applying the end-to-end WES workflow designed in
my master's thesis (*Precision Medicine & Clinical Genetics*): sample →
variant calling → annotation → filtering → prioritization (HPO) → ACMG/AMP
interpretation → clinical reporting.

## What it does

Given a VEP-annotated VCF (single sample) and a candidate-gene panel, it:

1. **QC filter** — drops variants failing `FILTER`, or below QUAL / depth / GQ cutoffs
2. **Frequency filter** — drops variants above the gnomAD MAF cutoff for the
   chosen inheritance model (default 1e-4 dominant, 0.05 recessive)
3. **Consequence filter** — keeps HIGH / MODERATE VEP IMPACT; reports what was removed and why
4. **Inheritance filter** — dominant (heterozygous) or recessive (homozygous) genotype logic
5. **Prioritization** — scores each surviving variant on: gene-panel membership,
   rarity/novelty, consequence severity, in-silico predictors (REVEL, CADD, SIFT,
   PolyPhen), and ClinVar significance — with every point of score backed by a
   human-readable evidence line
6. **Reporting** — ranked candidate table (TSV) + summary report

## Usage

```bash
python3 src/variant_pipeline.py \
  --vcf data/demo_annotated.vcf \
  --panel data/hpo_gene_panel.txt \
  --out-prefix results/candidates \
  --model dominant
```

Options: `--maf-dominant`, `--maf-recessive`, `--min-qual`, `--min-dp`,
`--min-gq`, `--model {dominant,recessive}`.

Try it on the included synthetic demo (23 variants, neurodevelopmental
gene panel) — no installation needed, standard library only:

```bash
python3 src/variant_pipeline.py --vcf data/demo_annotated.vcf \
  --panel data/hpo_gene_panel.txt --out-prefix results/candidates
```

Expected: 9 candidates ranked, top hit a ClinVar-pathogenic stop-gain in
*SMARCB1* with full evidence trail in `results/candidates.report.txt`.

## Project structure

```
src/variant_pipeline.py   # the pipeline (pure Python 3, no dependencies)
data/demo_annotated.vcf   # synthetic VEP-annotated demo VCF (23 variants)
data/hpo_gene_panel.txt   # candidate genes for the demo phenotype
results/                  # example output (TSV + report)
```

## Notes

- The demo VCF is fully synthetic, created to exercise every filter branch
  (low quality, common benign, synonymous, intronic, VUS, pathogenic).
- In a real setting, step 0 is annotation with Ensembl VEP; the pipeline
  consumes its standard `CSQ` output.
- Computational prioritization supports — but never replaces — expert
  ACMG/AMP interpretation.
