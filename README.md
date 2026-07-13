<div align="center">

# 🧬 PLit & ReThiN

### Stability-Based and Information-Theoretic Unsupervised Feature Selection Methods for Single-Cell RNA Sequencing

**Maitreya Sameer Ganu**
*Indian Institute of Science Education and Research (IISER), Thiruvananthapuram*

*Co-author & Advisor: Dr. Clint P. George — Indian Institute of Technology (IIT) Goa*

<br>

![Language](https://img.shields.io/badge/Language-R-276DC3?style=for-the-badge&logo=r&logoColor=white)
![Field](https://img.shields.io/badge/Field-Bioinformatics-00758F?style=for-the-badge)
![Topic](https://img.shields.io/badge/Topic-Feature%20Selection-orange?style=for-the-badge)
![Models](https://img.shields.io/badge/Models-Poisson%20%7C%20Negative%20Binomial-2E8B57?style=for-the-badge)
![Status](https://img.shields.io/badge/Status-Manuscript%20in%20Preparation-yellow?style=for-the-badge)

<br>

[Abstract](#project-abstract) • [Core Hypothesis](#core-hypothesis) • [Methodology](#methodology) • [Results](#results) • [Repository](#repository) • [Usage](#usage) • [Citation](#citation)

</div>

---

## Project Abstract

Feature selection is a critical preprocessing step in single-cell RNA sequencing (scRNA-seq), shaping downstream clustering, cell-type annotation, and every analysis built on top of it. This project introduces two unsupervised feature selection methods for count data:

- **PLit** (*Parametric Length Information Test*) — ranks genes by comparing the description length of their empirical count distribution against a fitted parametric null, using the Minimum Description Length (MDL) principle.
- **ReThiN** (*Reproducibility via Thinning*) — ranks genes by how reproducible their expression profile is under data thinning, measured via a split-half correlation.

Both methods are implemented for the **Poisson** and **Negative Binomial (NB)** count models — PLit extends to any parametric count family, and ReThiN to any convolution-closed count family — and both are deployed inside a single subsampling-based stability wrapper. We benchmark all four instances (PLit-Poisson, PLit-NB, ReThiN-Poisson, ReThiN-NB) against five established feature selectors (scran HVG, Seurat VST, Pearson Residuals, M3Drop, scry Deviance) and a random baseline, across **seven public scRNA-seq datasets** with ground-truth cell labels, scoring downstream k-means clustering with Adjusted Rand Index (ARI) and Normalized Mutual Information (NMI) across four feature budgets (K = 100, 200, 500, 1000). Both proposed methods are highly competitive with state-of-the-art variance-based approaches, achieving **top-3 or first-place** performance across the majority of dataset–metric combinations tested, despite using simpler models and fewer hyperparameters than any competing method.

## Core Hypothesis

Existing scRNA-seq feature selectors typically bury their assumptions inside preprocessing choices — a normalization strategy, a trend-fitting bandwidth, a clipping threshold, a fixed overdispersion constant — so that two "standard" selectors can disagree substantially even on identical data, and it's often unclear *which* assumption is driving *which* selected gene.

The hypothesis behind this project is that a feature selector should instead make exactly **one** explicit, statable, and swappable assumption about the count-generating process, and should return a score that estimates a **fixed population quantity** — not a score defined only relative to a trend that shifts with the dataset's gene composition.

- **PLit** operationalizes this via a *parametric null model*: it scores each gene by how much better its empirical distribution is explained by itself than by that null — a two-part MDL codelength gap that reduces to a penalized empirical KL divergence.
- **ReThiN** operationalizes this via a *count family closed under thinning*: it scores each gene by how reproducible its expression pattern is across two independently thinned halves of the data — a statistic that estimates a variance-components ratio.

Both are deployed inside the identical subsampling stability wrapper, so any performance difference observed in benchmarking reflects the discriminative power of the two proposed estimands, not an artifact of the evaluation protocol.

## Methodology

### Proposed Methods

**PLit (Parametric Length Information Test).** For each gene, PLit fits a parametric null model by maximum likelihood — Poisson (1 parameter) or Negative Binomial (2 parameters: mean + dispersion) — and compares its fitted log-likelihood against a fully empirical (saturated) log-likelihood built from the gene's observed count frequencies. Following Rissanen's (1978) two-part MDL construction, the empirical model is charged a complexity penalty proportional to its number of distinct observed count values. The unpenalized score is exactly `n × KL(empirical distribution ‖ fitted null)` — a penalized estimate of how far a gene's expression departs from its best-fitting null.

**ReThiN (Reproducibility via Thinning).** For each gene, ReThiN randomly splits every observed count into two independent halves using a fair thinning operator — Binomial(·, 0.5) for Poisson, Beta-Binomial for Negative Binomial — normalizes each half within-cell, and computes the correlation between the two halves across cells. This split-half correlation is a finite-sample estimate of a population variance-components ratio, `σ² / (σ² + 4w̄)`, that is exactly zero for genes whose variability is pure sampling noise and increases monotonically with genuine biological signal.

Both statistics generalize beyond Poisson and NB — PLit to any parametric count model, ReThiN to any convolution-closed count family.

### Benchmarking Protocol

Every method — proposed and baseline — is passed through an identical subsampling-and-aggregation wrapper, so any performance difference reflects the underlying score rather than the evaluation procedure:

1. Filter genes expressed in < 10 cells; compute library-size normalization factors once on the filtered matrix.
2. For 5 independent seeds, run 20 rounds of 80%-cell subsampling; apply each method's core score per round and aggregate the 20 rankings by average rank.
3. For each feature budget K ∈ {100, 200, 500, 1000}, select the top-K genes, log-normalize, run PCA (15 PCs), and cluster with k-means (30 seeds × 25 restarts) using the true number of populations.
4. Report mean ± SD of ARI and NMI over the 5 subsampling seeds.

### Methods Compared

| Method | Statistical Model | Input | Hyperparameters | Estimand |
|---|---|---|---|---|
| **PLit** † | Parametric null vs. empirical (MDL) | Raw counts | None | Penalised KL divergence from best-fit null |
| **ReThiN** † | Thinning + split-half correlation | Raw counts | n_thin | Variance ratio σ² / (σ² + 4w̄) |
| scran HVG | Mean–variance trend (LOESS) | Log-normalised | Trend span | None (residual above fitted trend) |
| Seurat VST | Mean–variance trend (LOESS) | Raw counts | Span, clip threshold | None (clipped standardised variance) |
| scry Deviance | Constant-proportion multinomial null | Raw counts | None | Multinomial deviance (unpenalised LR) |
| Pearson Residuals | Poisson/NB null, fixed overdispersion | Raw counts | Overdispersion const. | Residual variance (depth-corrected) |
| M3Drop | Michaelis–Menten dropout | Raw counts | None | None (dropout-vs-mean deviation) |

† Proposed method.

### Benchmark Datasets

| # | Dataset | System | Cells | Cell Types | Platform |
|---|---|---|---|---|---|
| 1 | Baron | Human pancreas | 8,569 | 13 | inDrop |
| 2 | Tian CellBench | Human cell lines | 895 | 3 | 10x Chromium |
| 3 | Zhengmix4eq | PBMC | 3,994 | 4 | 10x Chromium |
| 4 | Zhengmix8eq | PBMC | 3,994 | 8 | 10x Chromium |
| 5 | Segerstolpe | Human pancreas | 2,209 | 14 | Smart-seq2 |
| 6 | Darmanis | Human brain | 420 | 8 | Fluidigm C1 |
| 7 | Zeisel | Mouse cortex & hippocampus | 3,005 | 9 | STRT-Seq |

## Results

For each dataset, the table below shows which methods achieve the best **ARI** (top 3) and best **NMI** (top 1), using each method's **best score across all four feature budgets** (K = 100, 200, 500, 1000) tested, excluding the random baseline. The four rightmost columns show where PLit and ReThiN — under their Poisson and NB instances — land in that ranking.

**Legend:** ✅✅✅ = ranked **#1** on ARI or NMI for that dataset · ✅ = ranked in the **top 3** on ARI (or tied for the dataset's best) · — = outside the top 3 on both metrics

| Dataset | 🥇🥈🥉 Top 3 ARI | 🏆 Top NMI | PLit (Poisson) | PLit (NB) | ReThiN (Poisson) | ReThiN (NB) |
|---|---|---|---|---|---|---|
| **Segerstolpe**<br>(Human Pancreas) | 1. PLit (NB) — 0.625<br>2. ReThiN (Poisson) — 0.564<br>3. Seurat VST — 0.534 | ReThiN (Poisson) — 0.707 | — | ✅✅✅<br>0.625 ARI (K=1000) | ✅✅✅<br>0.707 NMI (K=1000) | — |
| **Darmanis**<br>(Human Brain) | 1. ReThiN (Poisson) — 0.921<br>2. PLit (NB) — 0.902<br>3. M3Drop — 0.881 | ReThiN (Poisson) — 0.882 | — | ✅<br>0.902 ARI (K=1000) | ✅✅✅<br>0.921 ARI / 0.882 NMI (K=1000) | — |
| **Tian CellBench**† | 8-way tie at 0.997 ARI (near-saturated dataset) | 8-way tie at 0.993 NMI | ✅<br>0.997 ARI, tied (K=500) | — | ✅<br>0.997 ARI, tied (K=500) | ✅<br>0.997 ARI, tied (K=500) |
| **Zhengmix4eq** | 1. scry Deviance — 0.930<br>2. PLit (Poisson) — 0.926<br>3. ReThiN (Poisson) — 0.907 | scry Deviance — 0.921 | ✅<br>0.926 ARI (K=500/1000) | — | ✅<br>0.907 ARI (K=1000) | — |
| **Zhengmix8eq** | 1. scry Deviance — 0.666<br>2. PLit (Poisson) — 0.664<br>3. ReThiN (Poisson) — 0.651 | scry Deviance — 0.752 | ✅<br>0.664 ARI (K=500/1000) | — | ✅<br>0.651 ARI (K=1000) | — |
| **Baron**<br>(Human Pancreas) | 1. Seurat VST — 0.619<br>2. ReThiN (Poisson) — 0.583<br>3. PLit (Poisson) — 0.579 | Seurat VST — 0.703 | ✅<br>0.579 ARI (K=500) | — | ✅<br>0.583 ARI (K=500) | — |
| **Zeisel**<br>(Mouse Brain) | 1. Pearson Residuals — 0.885<br>2. ReThiN (Poisson) — 0.882<br>3. M3Drop — 0.871 | ReThiN (Poisson) — 0.839 | — | — | ✅✅✅<br>0.839 NMI (K=1000) | — |

† *Tian CellBench saturates quickly: by K ≥ 200, nearly every method (proposed and baseline) reaches ARI ≈ 0.997 / NMI ≈ 0.993, so this dataset does not meaningfully discriminate between feature selectors.*

### Aggregate Ranking Across All 7 Datasets

Averaging each method's ARI/NMI across feature budgets, then ranking the 7 methods per dataset (1 = best), and averaging ranks across all 7 datasets:

| Method | Mean Rank (ARI) | Mean Rank (NMI) |
|---|---|---|
| **ReThiN** † | **2.57** 🥇 | **2.71** 🥇 |
| M3Drop | 3.14 | 3.00 |
| **PLit** † | 3.57 | 3.79 |
| scry Deviance | 3.71 | 3.43 |
| scran HVG | 4.64 | 4.21 |
| Pearson Residuals | 4.79 | 4.57 |
| Seurat VST | 5.57 | 6.29 |

† Proposed method.

> **A note on significance:** with only 7 benchmark datasets, pairwise significance tests are underpowered. Dataset-level bootstrap confidence intervals (10,000 resamples) show only two comparisons excluding zero on both metrics: ReThiN over Seurat VST, and ReThiN over Pearson Residuals. All other differences reported above — including PLit's and ReThiN's other wins — are directionally consistent point estimates, not certified effects at this sample size. See the manuscript for full details.

## Repository

All code, benchmarking scripts, and processed datasets used in this project are available at:

🔗 **https://github.com/MaitreyaGanu/PLit-ReThiN**

## Usage

```r
# Clone the repository
git clone https://github.com/MaitreyaGanu/PLit-ReThiN.git
cd PLit-ReThiN

# Install R dependencies
install.packages(c("Matrix", "MASS"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("scran", "scry", "M3Drop", "Seurat"))

# Run the benchmarking pipeline
# (see repository for dataset-specific entry points and the shared subsampling wrapper)
```

## Citation

```bibtex
@unpublished{ganu2026plitrethin,
  title        = {Stability-Based and Information-Theoretic Unsupervised Feature Selection Methods for Single-Cell RNA Sequencing},
  author       = {Ganu, Maitreya Sameer and George, Clint P.},
  year         = {2026},
  note         = {Manuscript in preparation},
  howpublished = {\url{https://github.com/MaitreyaGanu/PLit-ReThiN}}
}
```
