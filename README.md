# Clinical WES Variant Filtering & Prioritization Pipeline

Dependency-free exploratory, single-sample VEP-annotated WES variant filtering. **Research/portfolio demonstration only; not clinically validated and does not make ACMG/AMP classifications.**

## Processing

- Each ALT allele is handled independently; all matched CSQ/transcript annotations are retained in `all_csq_json` in its output row. CANONICAL and Feature are optional metadata and never filter an allele.
- QC requires FILTER=PASS, QUAL >= 30, DP >= 20, GQ >= 20 by default. At least one transcript must have HIGH/MODERATE impact.
- **Both AF lanes run in one pass, regardless of `--model` or genotype:**
  - Primary: known gnomAD AF < `1e-5` (strict).
  - Additional: `1e-5` <= known AF < `0.05` (extended). No duplicate between lanes.
  - AF >= `0.05`: excluded. Unknown AF remains in primary with `frequency_status=UNKNOWN` and is **not assumed rare**; invalid, incomplete or discordant annotation AF goes to review in the excluded/audit TSV.
- Inheritance compatibility and ALT dosage are descriptive metadata only. Neither heterozygous nor homozygous genotypes are excluded based on `--model`. Compound heterozygosity, segregation and X-linked rules are not assessed.
- Prioritization awards evidence once per ALT allele, never once per transcript. A HIGH/MODERATE annotation on any matched transcript may retain the allele; all LOW/MODIFIER annotations remain visible too.

## Run from repository root

```bash
python3 src/variant_pipeline.py \
  --vcf data/demo_annotated.vcf \
  --panel data/hpo_gene_panel.txt \
  --out-prefix results/candidates \
  --model dominant
```

`--model recessive` changes the descriptive `inheritance_compatibility` label, **not** which AF lane receives a variant. Optional `--maf-dominant` (default `1e-5`) sets the primary boundary and `--maf-recessive` (default `0.05`) sets the additional lane's upper boundary; both must be ordered correctly. These option names are retained for CLI compatibility.

## Independent outputs

| File | Content |
|---|---|
| `results/candidates.tsv` | Primary strict-AF candidates only (and explicitly flagged UNKNOWN AF candidates) |
| `results/candidates.additional_candidates.tsv` | Additional AF-lane candidates only; no duplicates with primary |
| `results/candidates.excluded.tsv` | QC, consequence, AF exclusions and annotation review cases, with reasons |
| `results/candidates.report.txt` | Counts for both lanes, exclusions and evidence for ranked candidates |

All three TSVs preserve every ALT-matched transcript annotation in `all_csq_json`. Output rows are one per ALT allele, not one per transcript. Accounting should satisfy: primary + additional + excluded/review = total ALT alleles processed.

## Tests

```bash
python3 -m unittest discover -s src -p "test_*.py" -v
```

The existing GitHub Actions workflow executes the same tests and runs the original synthetic 23-ALT VCF. It uploads the generated `results/candidates.*` artifact. Passing synthetic tests is not clinical validation. See the report and excluded TSV for any uncertainties; this program does not infer tissue relevance, gene-disease mechanism, clinical actionability or ACMG secondary findings.
