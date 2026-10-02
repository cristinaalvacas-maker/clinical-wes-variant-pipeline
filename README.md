# Clinical WES Variant Filtering & Prioritization Pipeline

Dependency-free exploratory Python pipeline for single-sample, VEP-annotated WES VCFs. It filters **ALT alleles, not transcripts**, preserves all matched CSQ annotations and scores each candidate only once. It is a research/portfolio tool, not a validated clinical classifier or ACMG/AMP interpretation engine.

## Workflow

1. Parse CSQ fields from the VCF header and separate multi-ALT records.
2. Apply QC (FILTER=PASS, QUAL >= 30, DP >= 20, GQ >= 20).
3. Apply gnomAD AF threshold (default dominant < 1e-5, recessive < 0.05). Missing AF is UNKNOWN, not zero; inconsistent or invalid annotations go to review.
4. Retain a variant if **any** matched transcript has HIGH or MODERATE impact; keep *all* matched transcript annotations, even LOW/MODIFIER. Neither Feature nor CANONICAL is required.
5. Annotate ALT dosage and simple inheritance compatibility **without excluding by genotype or inheritance**. This is descriptive only: no compound-heterozygosity, segregation, X-linked or gene-specific mechanism inference.
6. Prioritize each candidate ALT once using panel membership, rarity, impact and available prediction/ClinVar annotations. Scores are exploratory, not pathogenicity probabilities.

## Run

From the repository root (Python 3.10+; no external packages):

```bash
python3 src/variant_pipeline.py \
  --vcf data/demo_annotated.vcf \
  --panel data/hpo_gene_panel.txt \
  --out-prefix results/candidates \
  --model dominant
```

`--model recessive` changes the AF cutoff and the **informational** inheritance compatibility label; it does not discard alleles on genotype grounds. AF thresholds may be changed using `--maf-dominant` and `--maf-recessive`. Other parameters: `--min-qual`, `--min-dp`, `--min-gq`, `--transcript-override` (TSV columns `chrom pos ref alt transcript`). Transcript overrides affect reference metadata, not filtering or scoring.

## Three separate outputs

| Output | Contents |
| --- | --- |
| `results/candidates.tsv` | **Only** ranked candidates, one row per ALT; all matched CSQ transcript annotations in `all_csq_json`; genotype, ALT dosage and descriptive inheritance compatibility. |
| `results/candidates.excluded.tsv` | QC, frequency and consequence exclusions **plus unresolved review cases**, each with explicit `status` and `reason`; full annotations retained in `all_csq_json`. Review is not synonymous with exclusion. |
| `results/candidates.report.txt` | Run settings, counts, candidate evidence, audit summary and limitations. |

No `excluded_inheritance` status exists. An allele with ALT dosage different from the simple selected model remains a candidate if it passes QC, AF and consequence. The report does not automatically designate incidental/secondary findings; those require expert interpretation and appropriate clinical context.

**Important:** `--model dominant` still uses the *dominant AF threshold*. A rare homozygous allele below that threshold is retained, but an allele above it can still be excluded by frequency. Choose frequency settings consciously when looking for findings outside the primary hypothesis.

## Tests and GitHub Actions

```bash
python3 -m unittest discover -s src -p 'test_*.py' -v
```

The suite tests multi-ALT parsing, all-transcript handling, absent CANONICAL/Feature, frequency edge cases, genotype retention in both modes, output separation and score-once behavior. GitHub Actions runs these tests and the repository's original 23-ALT demo VCF. A successful workflow is a reproducibility check, **not** clinical validation.

## Limitations

Single-sample diploid exploratory analysis; simple SNV/indel allele matching requires further validation on representative VEP output. Ambiguous CSQ-to-ALT mapping and invalid/conflicting frequency are marked for review instead of guessed. This pipeline does not establish transcript tissue expression, disease mechanism, phase, compound heterozygosity, segregation, CNVs or pathogenicity classification.
