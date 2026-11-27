# hrdness

**English** | [日本語](README.ja.md)

Estimate **HRDness** (homologous-recombination-deficiency-like expression) in high-grade serous ovarian carcinoma (HGSOC)
from gene expression, using the 249-gene TCGA-OV HRDness signature.

> Signature: Takamatsu *et al.*, *Br J Cancer* 2023;128:1095–1104, doi:[10.1038/s41416-022-02122-9](https://doi.org/10.1038/s41416-022-02122-9)
> (UP 173 genes + DN 75 genes). `hrdness` is a re-implementation for applying the signature to new cohorts.

**For academic research use only (non-commercial). Not a clinical diagnostic.**

## What it does

| Output | Method | Use |
|---|---|---|
| **1. HRDness score** (continuous) | ssGSEA (rank-normalised, weight 0.75) of the UP and DN gene sets; score = UP − DN. Python port of the original R library (matches the R output). | Comparing samples **within a cohort**. Reference only when n < 50. |
| **2. HRDness probability** (continuous) + **annotation** (`HRD-high` / `HRD-low` / `undetermined`) | "Double z" (2Z): log2 expression of the 249 genes → z-score **within each sample** → z-score **per gene across samples** → logistic regression + random forest trained on TCGA-OV (HRD-high vs HRD-low extreme groups). | Per-sample call, relative to the cohort. Works for a **single sample** if saved cohort statistics are supplied. |

Annotation thresholds: probability ≥ 0.65 → `HRD-high`, ≤ 0.35 → `HRD-low`, otherwise `undetermined`.

Supported platforms: **RNA-seq** (default) and **Affymetrix** (`--platform affy`). Agilent is not supported yet.

## Read this before using the calls

The 2Z method makes the annotation **relative to the input cohort**.

- `HRD-high` means "more HRD-like than the rest of *this* cohort", not "HRD-positive". The high fraction is pulled toward
  the TCGA training mix (~50%) and is **not** a prevalence estimate. It does not replace an HRD score cut-off (e.g. 42) or a BRCA/HRD test.
- **Input must be HGSOC.** Other histologies, normal tissue, or very low-purity / stroma-rich samples distort the gene-wise
  statistics and therefore the calls.
- **Use ≥ 50 samples** (a practical guideline, not an empirically derived threshold). Below that, `hrdness` falls back to a
  sample-z-only mode (`1Z-single`) whose absolute probabilities shift between cohorts and platforms; treat it as reference only.
- **The cohort must be representative.** BRCA-mutant-only, HRD-high-only, responders-only, etc. will be forced toward a ~50/50 split.
  Analyse batches / sites separately, or correct them, when they are mixed.
- A single sample can be called **only** against statistics saved from ≥ 50 samples processed with the **same pipeline**
  (see [Single-sample calls](#single-sample-calls)). Statistics from a different platform or pipeline must not be reused.

## Installation

Requires Python 3 with `numpy`, `pandas`, `scipy`, `scikit-learn` (no R, no ComBat).
Developed and tested with Python 3.11, numpy 2.4, pandas 2.3, scipy 1.17, scikit-learn 1.8.

```bash
git clone <this repository>
cd hrd_gexp_signatures
pip install -r requirements.txt
python3 -m hrdness --help
```

The package must be run from (or installed alongside) the repository, because the TCGA reference tables are read from `data/`.

## Input

A tab-separated matrix, **genes × samples**:

- first column = HGSOC **gene symbols** (hg38 / HGNC), first row = sample IDs;
- **all genes**, not only the 249 signature genes (needed for ssGSEA; with < 5000 genes output 1 is skipped);
- duplicated symbols are collapsed to the row with the highest mean expression; missing values are filled with the row median.

| Platform | What to provide |
|---|---|
| RNA-seq (`--platform rnaseq`, default) | **TPM or FPKM** (linear scale). The tool applies log2(x + 1). If the matrix is already log2(x + 1), pass `--log`. Raw counts (integers) are rejected: normalise first. |
| Affymetrix (`--platform affy`) | **RMA-like log2** expression collapsed to gene symbols. Do **not** apply log2(x + 1). MAS5 values trigger a warning (negative values) and values > 50 are rejected. |

TPM vs FPKM, constant scaling and pseudocounts (0.01–1) have little effect because of the within-sample z-score.

## Usage

```bash
python3 -m hrdness <expr.tsv> <out_prefix> [--platform rnaseq|affy] [--log] [--save-stats FILE] [--stats FILE]
```

Example (RNA-seq TPM, 100 samples):

```bash
python3 -m hrdness my_cohort_tpm.tsv results/my_cohort --save-stats results/my_cohort_stats.tsv
# -> results/my_cohort_HRDness.tsv, results/my_cohort_stats.tsv
```

### Single-sample calls

1. Run a cohort of ≥ 50 samples with `--save-stats` (same platform, same quantification pipeline).
2. Call a new sample against those statistics:

```bash
python3 -m hrdness new_sample_tpm.tsv results/new_sample --stats results/my_cohort_stats.tsv
```

The result is identical to running the sample inside the reference cohort (`HRDness_mode = 2Z-stats`).
The call is relative to the reference cohort.

## Output: `<out_prefix>_HRDness.tsv`

| Column | Meaning |
|---|---|
| `ssGSEA_UP`, `ssGSEA_DN` | ssGSEA enrichment scores of the UP / DN gene sets |
| `HRDness_score` | UP − DN (output 1) |
| `HRDness_score_z` | z-score of `HRDness_score` within the cohort |
| `HRDness_score_reliable` | `True` if n ≥ 50 |
| `HRDness_score_group` | `high` / `low` split at the cohort median (n ≥ 50 only) |
| `LR_prob`, `RF_prob` | logistic regression / random forest (mean of 5 seeds) probabilities |
| `HRDness_probability` | mean of `LR_prob` and `RF_prob` (output 2) |
| `HRDness_annotation` | `HRD-high` (≥ 0.65) / `HRD-low` (≤ 0.35) / `undetermined` |
| `HRDness_mode` | `2Z-cohort` (n ≥ 50), `2Z-stats` (`--stats`), or `1Z-single` (n < 50 without stats; reference only) |
| `HRDness_probability_1Z` | probability from the within-sample z-score only; reference, shifts between cohorts |
| `n_signature_genes` | signature genes found in the input (of 249) |

Warnings (small n, odd scale, low gene overlap, extreme 1Z mean probability) are printed to stderr.

## How it works

1. Take the 249 signature genes present in both the input and the TCGA reference (≥ 200 required).
2. z-score across those genes within each sample.
3. z-score each gene across samples (cohort statistics, or saved statistics with `--stats`).
4. Score with LR + RF trained on TCGA-OV primary tumours (`-01` samples, one per case; HRD-high vs HRD-low extreme groups) that went
   through the *same* transformations. RNA-seq reference: GDC log2(FPKM + 1). Affy reference: Xena `HT_HG-U133A`.

Reference HRD labels: HRD-high = both HRD score and SigMA (signature 3) high; HRD-low = both low; discordant cases excluded
(thresholds from ROC against BRCA1/2 LOH). See the paper for details.

## Performance and caveats

- **TCGA internal CV (5×5)**: accuracy ≈ 0.88, AUC ≈ 0.94. This is **optimistic**: the 249 genes were selected as DEGs from the same
  TCGA HRD-high/low samples.
- **Independent cohorts (RNA-seq, 4 cohorts, trained on TCGA only)**: within-cohort Spearman with the cohort's HRD score ≈ 0.53, AUC for
  above-median HRD score ≈ 0.80. Discrimination is moderate. A low call does not guarantee the tumour is HR-proficient (in these cohorts
  25–45% of `HRD-low` calls still had HRD score < 42).
- Compared with the original ComBat-to-TCGA + RF/SVM/KNN/LR scheme (re-implemented), 2Z was equivalent or better on these cohorts and
  yields fewer `undetermined` calls. This is an exploratory comparison (4 cohorts, same TCGA label definition).
- Adding more classifiers (SVM, KNN, MLP, gradient boosting) did not improve external performance.
- **Affymetrix** support is verified for behaviour only (RMA vs MAS5 inputs, single-sample reproduction); accuracy against an
  independent HRD measure has **not** been validated. The Xena normalisation of the TCGA Affy reference is inferred from its value range (≈ 2.6–13).
- Platform-specific statistics: `--stats` files are specific to the platform and pipeline that produced them.

## Repository layout

```
hrdness/              package (core.py: 2Z + LR/RF, ssgsea.py: ssGSEA, __main__.py: CLI)
data/signature/       249-gene signature (UP / DN lists, gmt)
data/reference_TCGA/  TCGA-OV reference expression (249 genes) and HRD labels
requirements.txt      Python dependencies
NOTICE.md             third-party notices (ssGSEA license, TCGA data)
```

## Citation

If you use this tool, please cite the signature paper:

> Takamatsu S, Yoshihara K, Baba T, Shimada M, Yoshida H, Kajiyama H, Oda K, Mandai M, Okamoto A, Enomoto T, Matsumura N.
> Prognostic relevance of HRDness gene expression signature in ovarian high-grade serous carcinoma; JGOG3025-TR2 study.
> *Br J Cancer* 2023;128(6):1095–1104. doi:[10.1038/s41416-022-02122-9](https://doi.org/10.1038/s41416-022-02122-9)

## AI assistance

This software was developed with the assistance of Claude (Anthropic) through Claude Code, including code writing, evaluation scripts and
documentation drafting. The author (Shiro Takamatsu) directed the work, made the design decisions, and is responsible for the released content.
Commits made with Claude's assistance carry a `Co-Authored-By` trailer.

## License

**Academic / non-commercial use only.** The code is released under the
[PolyForm Noncommercial License 1.0.0](LICENSE). Research, education, and use by non-profit, public-research, educational and government
organisations are permitted; **commercial use is not permitted** without separate written permission from the copyright holder.
Note that this is a source-available licence, not an OSI-approved open-source licence.

- `hrdness/ssgsea.py` is a Python port of the Broad Institute ssGSEA code (BSD-style license). Its license text is reproduced in
  [NOTICE.md](NOTICE.md). Method: Barbie *et al.*, *Nature* 2009;462:108–112.
- The reference data in `data/reference_TCGA/` are derived from open-access TCGA data (GDC and UCSC Xena). Do not attempt to identify
  participants, and acknowledge TCGA / GDC / UCSC Xena as data sources in publications (see [NOTICE.md](NOTICE.md)).
