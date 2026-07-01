# PLit-ReThiN

Two stability-based unsupervised feature selectors for single-cell RNA-seq (scRNA-seq) count data, benchmarked against six established methods across seven annotated datasets.

Status: **draft / pre-publication**. Authors and code availability statement to be added.

---

## Table of Contents

- [Overview](#overview)
- [Methods](#methods)
  - [PLit — Poisson Length Information Test](#plit--poisson-length-information-test)
  - [ReThiN — Reproducibility via Thinning](#rethin--reproducibility-via-thinning)
- [Repository Structure](#repository-structure)
- [Requirements](#requirements)
- [Usage](#usage)
- [Benchmarking Framework](#benchmarking-framework)
- [Datasets](#datasets)
- [Results Summary](#results-summary)
- [Limitations](#limitations)
- [Citation](#citation)
- [License](#license)
- [Contact](#contact)

---

## Overview

Established scRNA-seq feature selectors (scran HVG, Seurat VST, Pearson residuals) rank genes via a mean–variance trend or a fixed dispersion constant, decisions that are somewhat arbitrary and can silently shift the selected gene set. PLit and ReThiN instead operate directly on raw counts with at most one non-defining hyperparameter.

| Property | PLit | ReThiN |
|---|---|---|
| Core idea | MDL code-length gap: single-rate Poisson vs. empirical count distribution | Split-half correlation under Poisson (binomial) thinning |
| Input | Raw counts | Raw counts |
| Hyperparameters | None | `n_thin` (repetitions; does not change what is measured) |
| Closest baseline | scry Deviance | Molecular cross-validation / count splitting (Batson et al. 2019) |
| Stability mechanism | Bootstrap subsampling (80% cells, no replacement) + rank averaging | Same, plus averaging over `n_thin` thinning draws |

---

## Methods

### PLit — Poisson Length Information Test

For each gene, compares two-part MDL description lengths:

- **Null model** `M0`: single-rate Poisson, MLE `λ_j = (1/n) Σ_i x_ji`
- **Alternative** `M1`: fully empirical (saturated) categorical distribution of observed counts, `V_j` distinct values, `V_j − 1` free parameters

```
S_j = L̃_1,j − L_0,j − ((V_j − 2) / 2) · ln(n)
```

Genes well-fit by a single Poisson rate score near zero/negative; genes with structured (multimodal, bursty) expression score positively and rank higher. Ranks are averaged over `B` bootstrap rounds (80% cell subsampling, no replacement).

### ReThiN — Reproducibility via Thinning

Each count `x_ji ~ Poisson(λ_ji)` is split via `A_ji ~ Binomial(x_ji, 0.5)`, `B_ji = x_ji − A_ji`. Under Poisson thinning, `A_ji ⊥ B_ji`, each marginally `Poisson(λ_ji/2)`.

- **Lemma:** thinning splits are independent Poisson(λ/2) variables.
- **Proposition:** `Corr(A_ji, B_ji) = σ_j² / (σ_j² + 2μ_j)`, where `μ_j, σ_j²` are the mean/variance of `λ_ji` across cells.

Correlation is 0 when a gene's rate is constant across cells (pure Poisson noise) and increases with the biological signal-to-noise ratio `σ_j²/μ_j`. The algorithm correlates within-cell-normalized splits, averaged over `n_thin` thinning repeats and `B` bootstrap rounds.

> **Known open gap:** the identity above is exact for raw counts; a finite-sample bound on the approximation error introduced by within-cell normalization is not yet derived (see [Limitations](#limitations)).

---

## Repository Structure

```
PLit-ReThiN/
├── BaronPancreas/
│   ├── BaronPancreas.R
│   ├── benchmark_summary_BaronPancreas.csv
│   ├── fig_ARI_BaronPancreas.pdf
│   ├── fig_NMI_BaronPancreas.pdf
│   ├── fig_Runtime_BaronPancreas.pdf
│   ├── raw_results_BaronPancreas.csv
│   ├── runtime_BaronPancreas.csv
│   └── session_info_BaronPancreas.txt
│
├── DarmanisHumanBrain/
│   ├── DarmanisBrain.R
│   ├── benchmark_summary_DarmanisBrain.csv
│   ├── fig_ARI_DarmanisBrain.pdf
│   ├── fig_NMI_DarmanisBrain.pdf
│   ├── fig_Runtime_DarmanisBrain.pdf
│   ├── raw_results_DarmanisBrain.csv
│   ├── runtime_DarmanisBrain.csv
│   └── session_info_DarmanisBrain.txt
│
├── SegerstolpePancreas/
│   ├── SegerstolpePancreas.R
│   ├── benchmark_summary_SegerstolpePancreas.csv
│   ├── fig_ARI_SegerstolpePancreas.pdf
│   ├── fig_NMI_SegerstolpePancreas.pdf
│   ├── fig_Runtime_SegerstolpePancreas.pdf
│   ├── raw_results_SegerstolpePancreas.csv
│   ├── runtime_SegerstolpePancreas.csv
│   └── session_info_SegerstolpePancreas.txt
│
├── TianCellBench/
│   ├── TianCellBench.R
│   ├── benchmark_summary_TianCellBench.csv
│   ├── fig_ARI_TianCellBench.pdf
│   ├── fig_NMI_TianCellBench.pdf
│   ├── fig_Runtime_TianCellBench.pdf
│   ├── raw_results_TianCellBench.csv
│   ├── runtime_TianCellBench.csv
│   └── session_info_TianCellBench.txt
│
├── ZeiselBrain/
│   ├── ZeiselBrain.R
│   ├── benchmark_summary_ZeiselBrain.csv
│   ├── fig_ARI_ZeiselBrain.pdf
│   ├── fig_NMI_ZeiselBrain.pdf
│   ├── fig_Runtime_ZeiselBrain.pdf
│   ├── raw_results_ZeiselBrain.csv
│   ├── runtime_ZeiselBrain.csv
│   └── session_info_ZeiselBrain.txt
│
├── ZhengMix4eq/
│   ├── ZhengMix4eq.R
│   ├── benchmark_summary_ZhengMix4eq.csv
│   ├── fig_ARI_ZhengMix4eq.pdf
│   ├── fig_NMI_ZhengMix4eq.pdf
│   ├── fig_Runtime_ZhengMix4eq.pdf
│   ├── raw_results_ZhengMix4eq.csv
│   ├── runtime_ZhengMix4eq.csv
│   └── session_info_ZhengMix4eq.txt
│
├── ZhengMix8eq/
│   ├── ZhengMix8eq.R
│   ├── benchmark_summary_ZhengPBMC.csv
│   ├── fig_ARI_ZhengPBMC.pdf
│   ├── fig_NMI_ZhengPBMC.pdf
│   ├── fig_Runtime_ZhengPBMC.pdf
│   ├── raw_results_ZhengPBMC.csv
│   ├── runtime_ZhengPBMC.csv
│   └── session_info_ZhengPBMC.txt
│
├── LICENSE
└── README.md
```



---

## Requirements

- R ≥ 4.x
- Packages: `scran`, `Seurat`, `scry`, `M3Drop`, standard stats/plotting stack
- Note: BLAS threading is uncontrolled across methods — a known limitation of the current runtime comparisons (not correctness-affecting).

---

## Benchmarking Framework

Identical protocol applied to every dataset/method (no label leakage — ground truth used only for evaluation):

| Setting | Value |
|---|---|
| Gene filtering | Expressed in ≥ 10 cells |
| Normalization | Library-size factors, computed once on full filtered matrix |
| Bootstrap seeds | 5 |
| Bootstrap rounds/seed | 20 |
| Cell subsampling | 80%, no replacement |
| Rank aggregation | Average rank over 20 rounds/seed |
| Feature budgets (K) | 100, 200, 500, 1000 |
| Dimension reduction | PCA, first 15 PCs |
| Clustering | k-means, K_true centers, 25 restarts, 30 eval seeds |
| Reported stats | Mean ± SD over 5 bootstrap seeds |
| Metrics | ARI, NMI |
| n_thin (ReThiN) | 5 |

---

## Datasets

| Dataset | System | Cells | Types | Platform | Ground Truth |
|---|---|---|---|---|---|
| Baron | Human pancreas | 8,569 | 13 | inDrop | Expert-curated (clustering + marker annotation) |
| Tian CellBench | Human cell lines | 895 | 3 | 10X | Experimentally fixed (cultured separately, mixed) |
| Zhengmix4eq | PBMC | 3,994 | 4 | 10X | Kit-purified, computationally mixed |
| Zhengmix8eq | PBMC | 3,994 | 8 | 10X | Kit-purified, computationally mixed |
| Segerstolpe | Human pancreas | 2,209 | 14 | Smart-seq2 | Expert-curated |
| Darmanis | Human brain | 420 | 8 | Fluidigm C1 | Expert-curated |
| Zeisel | Mouse cortex/hippocampus | 3,005 | 9 | STRT-Seq | Expert-curated |

---

## Results Summary

Mean rank across 7 datasets (1 = best, 7 = worst; Random baseline excluded):

| Method | Mean rank (ARI) | Mean rank (NMI) |
|---|---|---|
| **ReThiN** † | **2.57** | **2.71** |
| M3Drop | 3.14 | 3.00 |
| **PLit** † | 3.57 | 3.79 |
| scry Deviance | 3.71 | 3.43 |
| scran HVG | 4.64 | 4.21 |
| Pearson Residuals | 4.79 | 4.57 |
| Seurat VST | 5.57 | 6.29 |

† proposed method.

At `n = 7` datasets, Wilcoxon/Friedman/sign tests are underpowered; only two proposed-vs-baseline comparisons have bootstrap 95% CIs excluding zero (ReThiN vs. Seurat VST, ReThiN vs. Pearson Residuals). All other comparisons are directionally consistent, not statistically confirmed because of low power of statistical tests at `n = 7` datasets.

---

## Limitations

- No pairwise significance testing at `n = 7` datasets (see Cross-Dataset Statistics above) — mean rank / median effect size reported instead, with bootstrap CIs where available.
- ReThiN's within-cell-normalized split-half correlation is an asymptotic approximation to the exact unnormalized identity; a finite-sample error bound is not yet derived.
- PLit's Poisson null has no explicit per-cell depth offset.
- BLAS threading uncontrolled across methods (runtime comparisons only, not accuracy).
- Benchmarked on 7 datasets; broader coverage would strengthen generalizability claims.

---

## Citation

```
TBD — to be added upon publication.
```

---

## License

MIT License

Copyright (c) 2026 MaitreyaGanu

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

---

## Contact

TBD — authors and correspondence to be added.
