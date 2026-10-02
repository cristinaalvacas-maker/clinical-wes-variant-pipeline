# Clinical WES Variant Filtering & Prioritization Pipeline

A lightweight, dependency-free Python pipeline for **exploratory, single-sample whole-exome sequencing (WES) variant filtering and prioritization**. It processes VEP-annotated VCF data, preserves transcript-level annotations, applies quality, allele-frequency, consequence and simple inheritance filters, and ranks candidates with a human-readable evidence trail.

This portfolio project implements part of the end-to-end WES workflow designed in my master's thesis (*Precision Medicine & Clinical Genetics*): sample → variant calling → annotation → filtering → phenotype-informed prioritization (HPO) → expert ACMG/AMP interpretation → clinical reporting.

**Scope:** This program supports exploratory research and review. It is not a validated clinical diagnostic tool, does not classify variants under ACMG/AMP, and does not replace expert interpretation.

## What the pipeline does

1. **Reads VEP annotations** from the VCF `CSQ` header instead of assuming a fixed field order. Each ALT allele is processed separately; unmatched or ambiguous allele-to-CSQ mappings are flagged for review rather than guessed.
2. **Preserves every matched transcript annotation**, including those for excluded variants. The canonical transcript is the default annotation for scoring, not necessarily the clinically relevant transcript.
3. **Applies QC**: only `FILTER=PASS`, with configurable minimum QUAL, depth (DP) and genotype quality (GQ).
4. **Applies allele-frequency filtering** using gnomAD AF and model-specific thresholds. Unknown frequency is distinguished from AF = 0; invalid, incomplete or conflicting transcript AF values are flagged for review.
5. **Checks consequences across all ALT-matched transcripts**. A LOW/MODIFIER canonical transcript cannot silently remove an allele with a HIGH/MODERATE alternative. Such a case enters the review queue without an invented score. A selected MODERATE transcript with a HIGH alternative is retained with a review flag.
6. **Applies a simple inheritance model**: one ALT copy for the dominant/heterozygous mode, two copies for the recessive/homozygous mode.
7. **Prioritizes candidates** using candidate-gene panel membership, rarity, selected-transcript impact, REVEL, CADD, SIFT, PolyPhen and ClinVar annotations. Each awarded point has a readable evidence line. The panel may be informed by HPO-based phenotype assessment; the pipeline itself does not compute HPO similarity.
8. **Reports** ranked candidates, a full transcript review table and a summary of review cases and exclusions.

## Requirements

- Python **3.10 or newer**; Python standard library only (no `pip install` required).
- A **single-sample, VEP-annotated VCF** with a `CSQ` header defining at least: `Allele`, `Feature`, `CANONICAL`, `IMPACT`, `Consequence` and `gnomAD_AF`.
- A plain-text candidate-gene panel containing one gene symbol per line.

Other CSQ fields, including `SYMBOL`, `MANE_SELECT`, `REVEL`, `CADD_PHRED`, `SIFT`, `PolyPhen` and `ClinVar_CLNSIG`, are used or displayed when supplied. Missing optional annotations cannot contribute evidence.

## Quick start

Run commands from the **repository root**:

```bash
python3 src/variant_pipeline.py \
  --vcf data/demo_annotated.vcf \
  --panel data/hpo_gene_panel.txt \
  --out-prefix results/candidates \
  --model dominant
```

For a recessive/homozygous exploratory pass:

```bash
python3 src/variant_pipeline.py \
  --vcf data/demo_annotated.vcf \
  --panel data/hpo_gene_panel.txt \
  --out-prefix results/recessive \
  --model recessive
```

**Important:** the included 23-variant synthetic demonstration VCF was originally created for the earlier version of the pipeline. Before relying on its output with this revision, check that its VEP `CSQ` header contains the newly required transcript and allele fields. The former README's result of nine candidates and a top SMARCB1 hit has **not been revalidated** against this revised implementation.

### Defaults and options

| Option | Default | Meaning |
| --- | --- | --- |
| `--model` | `dominant` | `dominant` (one ALT copy) or `recessive` (two ALT copies) |
| `--maf-dominant` | `1e-5` | Dominant AF cutoff; values at or above it are excluded |
| `--maf-recessive` | `0.05` | Recessive AF cutoff; values at or above it are excluded |
| `--min-qual` | `30` | Minimum variant QUAL |
| `--min-dp` | `20` | Minimum read depth |
| `--min-gq` | `20` | Minimum genotype quality |
| `--out-prefix` | `results/candidates` | Prefix for output files |
| `--transcript-override` | none | Optional TSV with manual ENST selections per variant |

Only `FILTER=PASS` is accepted. The output directory is created automatically.

## Output files

With `--out-prefix results/candidates`, the pipeline produces:

| File | Contents |
| --- | --- |
| `results/candidates.tsv` | Ranked candidates, selected transcript, heuristic score and supporting evidence |
| `results/candidates.transcript_review.tsv` | Every parsed ALT and every matched transcript annotation, **including excluded alleles**, with selection status, consequence and review/exclusion reason |
| `results/candidates.report.txt` | Counts, settings, candidate evidence, review queue, exclusions and limitations |

The transcript review TSV can be filtered by `variant_status`, `IMPACT`, `CANONICAL`, `selected` and `selection_status`. An allele requiring transcript review is retained there and is **not assigned a candidate score** until a suitable annotation is selected and the pipeline is rerun.

## Reviewing alternative transcripts and rerunning

The pipeline initially selects the VEP canonical transcript when available, but **retains all annotations** so that a reviewer can examine alternative ENSTs, consequences and optional MANE Select information. It does **not** infer tissue expression, gene–disease mechanism or the most clinically relevant transcript automatically.

If the canonical transcript is missing, or a relevant alternative warrants investigation, review `results/candidates.transcript_review.tsv`. After evaluating the evidence externally, create a **tab-separated** file, for example `overrides.tsv`:

```text
chrom	pos	ref	alt	transcript
1	101	A	G	ENST00000000002
```

The row above is illustrative, **not a real recommended transcript**. Use the exact chromosome, position, REF, ALT and ENST values from your VCF/review table. Each row specifies one allele-level variant; different variants may have different selected transcripts.

Rerun automatically using the chosen annotation:

```bash
python3 src/variant_pipeline.py \
  --vcf data/demo_annotated.vcf \
  --panel data/hpo_gene_panel.txt \
  --model dominant \
  --transcript-override overrides.tsv \
  --out-prefix results/reviewed
```

A requested ENST absent from the annotations matched to that ALT remains in the review queue, with an explanation. An override key that matches no VCF allele causes an explicit input error to help catch mistakes.

## Project structure

```text
src/
  variant_pipeline.py       # Transcript-aware exploratory WES pipeline
  test_variant_pipeline.py  # Automated synthetic tests
data/
  demo_annotated.vcf        # Original synthetic demo; verify revised CSQ compatibility
  hpo_gene_panel.txt        # Demo candidate-gene panel
results/                   # Generated candidate, transcript-review and report files
README.md                  # Project documentation
```

## Automated tests

From the repository root:

```bash
cd src
python3 -m unittest -v test_variant_pipeline.py
```

The synthetic test suite covers CSQ field-order validation, canonical/alternative transcript discrepancies, missing canonical annotations, manual overrides, multi-ALT mapping and genotype dosage, dominant/recessive thresholds, AF edge cases, indels, ambiguous annotations and output accounting. Passing synthetic tests does **not** constitute clinical validation; representative real VEP outputs and complex indels need further assessment.

## Scientific and implementation limitations

- **Transcript choice:** VEP `CANONICAL` is only a default. MANE, tissue-specific expression, phenotype, gene–disease associations and disease mechanism require expert review. Alternative HIGH impact does not automatically imply clinical relevance or pathogenicity.
- **Inheritance:** Only simple diploid, single-sample heterozygous dominant and homozygous recessive ALT dosage are modeled. No compound heterozygosity, specific X-linked handling, family segregation, CNVs or automatic gene-specific inheritance assignment. Running both modes is not comprehensive inheritance analysis.
- **Allele mapping:** Multi-ALT records are separated and VEP allele representations are matched conservatively for simple SNVs/anchored indels. Ambiguous or unmatched annotations are flagged; complex normalization should be checked independently.
- **Frequency:** AF must be concordant across ALT-matched annotations. Completely missing AF is reported as `UNKNOWN`, not zero; incomplete, invalid or conflicting values trigger review.
- **Filtering and scoring:** These are exploratory heuristics, not pathogenicity probabilities or ACMG/AMP classifications. A variant excluded under a selected model should not be interpreted as clinically irrelevant. All parsed ALT/transcript annotations remain auditable in the transcript review output.
- **Validation:** Synthetic tests exercise selected scenarios but do not establish diagnostic sensitivity, specificity or fitness for clinical deployment.

The intended workflow is **automated filtering → transparent transcript review → optional documented ENST override → automatic rerun → expert variant interpretation**.
