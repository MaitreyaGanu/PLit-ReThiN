<div align="center">

# 🧬 PLit & ReThiN 🔬

### Unsupervised Feature Selection for scRNA-seq Count Data via Description Length and Data Thinning

**<u>Maitreya Sameer Ganu</u>**<br>
*Indian Institute of Science Education and Research (IISER) Thiruvananthapuram*<br>

Advisor: **Dr. Clint P. George**<br>
*Indian Institute of Technology (IIT) Goa*

<br>

<a href="https://github.com/MaitreyaGanu/PLit-ReThiN/stargazers"><img src="https://img.shields.io/github/stars/MaitreyaGanu/PLit-ReThiN?style=for-the-badge&logo=github&color=gold" alt="GitHub stars"/></a>
<img src="https://img.shields.io/badge/Language-R-276DC3?style=for-the-badge&logo=r&logoColor=white" alt="R"/>
<img src="https://img.shields.io/badge/Field-Bioinformatics-00758F?style=for-the-badge" alt="Bioinformatics"/>
<img src="https://img.shields.io/badge/Topic-Unsupervised%20Feature%20Selection-orange?style=for-the-badge" alt="Unsupervised feature selection"/>
<img src="https://img.shields.io/badge/Status-Manuscript%20in%20Preparation-yellow?style=for-the-badge" alt="Status"/>

<br><br>

[Overview](#overview) • [Key findings](#key-findings) • [Methods](#methods) • [Theory](#theory) • [Benchmark](#benchmark-setup) • [Results](#results) • [Limitations](#limitations) • [Installation](#installation-and-usage)

</div>

## Overview

Feature selection is a standard step before dimensionality reduction and clustering of single-cell RNA sequencing (scRNA-seq) data. Most unsupervised selectors score a gene against a trend fitted across all genes, or reward any departure from a null model without accounting for the complexity needed to describe it. This project proposes two feature-selection methods for count data, each of which targets an explicitly defined per-gene quantity:

- **PLit** (*Parametric Length Information Test*) uses the minimum description length (MDL) principle to compare a Poisson or negative-binomial (NB) null with the empirical distribution of a gene's counts. Its score is a complexity-penalized empirical Kullback–Leibler divergence.
- **ReThiN** (*Reproducibility via Thinning*) splits each count into two independent halves by Poisson data thinning and scores a gene by how well its variation between cells reproduces across the two halves.

We compare three proposed instances, **PLit (Poisson)**, **PLit (NB)** and **ReThiN**, with five established selectors (scran HVG, Seurat VST, scry Deviance, analytic Pearson residuals, M3Drop) and a random baseline. The benchmark covers **seven annotated scRNA-seq datasets** and four feature budgets (K = 100, 200, 500, 1000), and scores *k*-means clustering by ARI and NMI. Every method is run through the same subsampling-and-rank-aggregation wrapper, on the same cell subsamples.

## Key findings

> These are stated plainly, including where the proposed methods do **not** win. With seven datasets, all differences are descriptive; no significance test has been applied.

**ReThiN is consistently near the top.**
- It has the **best mean rank of all eight methods** in both ARI (3.20) and NMI (3.12), ahead of scry Deviance (3.91 / 3.73).
- It is in the **top three on all seven datasets** (mean ARI rank over the four budgets), but **first on none**.
- Across the 28 dataset–budget combinations, it is in the top three in 19 for ARI and never ranks below fifth.
- Its lead grows with K. Its average ARI gain over random selection is the largest of all methods at K = 500 and K = 1000 (0.219 and 0.188), and within 0.001 of the largest at K = 100 and 200.

**PLit (Poisson) behaves like scry Deviance.**
- Its ARI is within 0.037 of scry's in 27 of 28 dataset–budget combinations.
- It is first or second at every budget on Zhengmix4eq, and at K ≥ 200 on Zhengmix8eq.

**PLit's null must match the technical noise.**
- On the **four informative UMI datasets**, PLit (Poisson) has a mean ARI rank of 3.38. PLit (NB) has 7.56 and falls below random on Zeisel at every budget.
- On the **two read-count datasets** (Segerstolpe, Darmanis) the order reverses: PLit (NB) ranks **1.25**, the best of all methods, while PLit (Poisson) ranks 5.62.
- UMI counts are close to Poisson, so a negative binomial with its own dispersion absorbs real variation between cells. Read counts carry extra overdispersion from amplification, so there the negative binomial is the better reference. This rests on only two read-count datasets, which also differ in platform, depth and cell number.

**Informed selection matters most at small budgets.** The ARI gain over random selection falls from 0.25–0.37 at K = 100 to 0.04–0.19 at K = 1000.

## Methods

**PLit.** For each gene, PLit fits a parametric null by maximum likelihood: Poisson (rate) or negative binomial (mean and size, with the size found numerically on [10⁻³, 10⁶]). It compares the null's log-likelihood with the log-likelihood of the gene's empirical count distribution. Following Rissanen's two-part MDL code, the empirical model is charged (½ ln n) per extra parameter, where it has V − 1 parameters for V distinct count values. Before the penalty, the score is exactly n × KL(empirical ‖ fitted null).

**ReThiN.** Each count is split into two halves with Binomial(x, ½), so that under Poisson noise the halves are independent given the cell's expected expression. Each half is divided by its cell total, and the score is the Pearson correlation of the two halves across cells, averaged over n_thin = 5 thinning repeats. Genes whose variation is only sampling noise score about zero; genes whose expected expression differs between cells score higher.

> An earlier version also had a negative-binomial instance of ReThiN. It has been removed: when the dispersion is estimated from the same data, it absorbs the biological variation the score is meant to detect.

### Methods compared

| Method | Statistical model | Input | Tuning parameters | What the score targets |
|---|---|---|---|---|
| **PLit (Poisson)** † | Poisson null vs. empirical distribution (MDL) | Raw counts | None | Penalized KL divergence from the fitted null |
| **PLit (NB)** † | NB null vs. empirical distribution (MDL) | Raw counts | None | Penalized KL divergence from the fitted null |
| **ReThiN** † | Poisson thinning + split-half correlation | Raw counts | n_thin | σ² / (σ² + 2μ) |
| scran HVG | Mean–variance trend | Log-normalized | Trend settings | Variance above the fitted trend |
| Seurat VST | LOESS of log-variance on log-mean | Raw counts | Span, clip value | Clipped standardized variance |
| scry Deviance | Constant-proportion (binomial) null | Raw counts | None | Unpenalized deviance |
| Pearson residuals | NB null with depth offset, fixed θ = 100 | Raw counts | θ | Residual variance |
| M3Drop | Michaelis–Menten dropout model | Raw counts | None | Excess zeros for a gene's mean |

† Proposed method.

## Theory

**PLit.** With $V_j$ distinct observed counts and a null with $d_0$ parameters,

$$S_j = n\,\widehat{\mathrm{KL}}\left(\hat p_j \,\middle\|\, f(\cdot;\hat\theta_j)\right) - \frac{(V_j-1)-d_0}{2}\ln n .$$

| Instance | Penalty | Null parameters |
|---|---|---|
| Poisson ($d_0 = 1$) | $\tfrac{V_j-2}{2}\ln n$ | rate $\lambda_j$ |
| Negative binomial ($d_0 = 2$) | $\tfrac{V_j-3}{2}\ln n$ | mean $\mu_j$, size $r_j$ |

**ReThiN.** If $X \mid Z \sim \mathrm{Poisson}(\mu_Z)$ and the two halves are obtained by Binomial(·, ½) thinning, then, with $\mu_j$ and $\sigma_j^2$ the mean and between-cell variance of a gene's expected expression,

$$\mathrm{Corr}(A_{ji}, B_{ji}) = \frac{\sigma_j^2}{\sigma_j^2 + 2\mu_j}.$$

This is zero exactly when the gene's expected expression does not vary between cells ($\sigma_j^2 = 0$), and it increases with $\sigma_j^2/\mu_j$. The identity holds for the halves before normalization under constant expected depth. The implemented score normalizes within each cell, so it approximates this quantity; no finite-depth error bound is given.

Full derivations are in the appendix of the manuscript.

## Benchmark setup

1. Remove genes expressed in fewer than 10 cells, then compute library-size factors once on the filtered matrix.
2. For each of 5 subsampling seeds, draw 20 subsamples of 80% of the cells **without replacement**. Score genes on each subsample and aggregate the 20 rankings by average rank. Every method receives the same subsamples.
3. For each budget K ∈ {100, 200, 500, 1000}, take the top-K genes, log-normalize, keep 15 principal components, and run *k*-means with the true number of populations (30 seeds × 25 restarts).
4. Report the mean ± SD of ARI and NMI over the 5 subsampling seeds.

> The wrapper is inspired by stability selection, but it aggregates ranks over 80% subsamples. The false-selection guarantees of stability selection therefore do **not** apply; the wrapper is a variance-reduction and fair-comparison device.

### Datasets

| # | Dataset | System | Genes* | Cells | Types | Platform | Counts | Ground truth |
|---|---|---|---|---|---|---|---|---|
| 1 | Baron | Human pancreas | 15,117 | 8,569 | 14 | inDrop | UMI | Marker-validated clustering |
| 2 | Tian (CellBench) | Human cell lines | 16,208 | 902 | 3 | 10x Chromium | UMI | Pure cell-line identity |
| 3 | Zhengmix4eq | PBMC | 10,434 | 3,994 | 4 | 10x Chromium | UMI | Kit-purified, computationally mixed |
| 4 | Zhengmix8eq | PBMC | 10,600 | 3,994 | 8 | 10x Chromium | UMI | Kit-purified, computationally mixed |
| 5 | Zeisel | Mouse cortex & hippocampus | 16,484 | 3,005 | 7 | STRT-Seq | UMI | Marker-based expert annotation |
| 6 | Segerstolpe | Human pancreas | 18,992 | 2,209 | 14 | Smart-seq2 | Reads | Marker-based expert annotation |
| 7 | Darmanis | Human brain | 15,102 | 285 | 6 | Fluidigm C1 | Reads | Marker-based expert annotation |

\* After removing genes expressed in fewer than 10 cells.

## Results

<div align="center">
<img src="https://raw.githubusercontent.com/MaitreyaGanu/PLit-ReThiN/main/combined%20summary/fig_rank_summary.png"
     alt="Rank of each method within each dataset and budget"
     width="100%"></div>

*Rank of each method (1 = best; random baseline excluded) within each dataset and budget, for ARI (top) and NMI (bottom). Datasets left of the black line contain UMI counts, those to the right read counts. The last column is the mean rank over all 28 dataset–budget combinations.*

### Mean rank across datasets and budgets

| Method | ARI (all 7) | NMI (all 7) | ARI, UMI (4 informative) | ARI, read counts (2) |
|---|---|---|---|---|
| **ReThiN** † | **3.20** 🥇 | **3.12** 🥇 | 3.31 | 2.69 |
| scry Deviance | 3.91 | 3.73 | 3.31 | 4.38 |
| **PLit (Poisson)** † | 4.09 | 4.21 | 3.38 | 5.62 |
| M3Drop | 4.21 | 3.84 | 4.81 | 3.81 |
| Pearson residuals | 4.61 | 4.82 | 4.06 | 6.06 |
| scran HVG | 5.00 | 4.84 | 4.19 | 7.19 |
| Seurat VST | 5.27 | 5.52 | 5.38 | 5.00 |
| **PLit (NB)** † | 5.71 | 5.91 | 7.56 | **1.25** 🥇 |

† Proposed method. Ranks are among the eight methods (1 = best), with ties given the average rank. Leaving out Tian, where all methods are at ceiling, does not change the order.

### Per dataset (mean ARI rank over the four budgets)

| Dataset | Counts | Top three | PLit (Poisson) | PLit (NB) | ReThiN |
|---|---|---|---|---|---|
| Zhengmix4eq | UMI | scry Deviance, PLit (Poisson), ReThiN | 1.75 | 6.50 | 3.25 |
| Zhengmix8eq | UMI | scry Deviance, PLit (Poisson), ReThiN | 2.00 | 8.00 | 3.50 |
| Baron | UMI | Pearson residuals, Seurat VST, ReThiN | 4.75 | 7.75 | 3.25 |
| Zeisel | UMI | M3Drop, scran HVG, ReThiN | 5.00 | 8.00 | 3.25 |
| Tian‡ | UMI | M3Drop, ReThiN, then a three-way tie (PLit (Poisson), Pearson residuals, scran HVG) | 3.88 | 7.25 | 3.75 |
| Segerstolpe | Reads | PLit (NB), ReThiN, M3Drop | 5.50 | **1.00** | 2.62 |
| Darmanis | Reads | PLit (NB), ReThiN, Seurat VST | 5.75 | **1.50** | 2.75 |

‡ Tian is at ceiling: every informed method reaches ARI ≥ 0.987 at every budget, and even random selection reaches ARI ≥ 0.969, so its ranks reflect differences of at most 0.02.

Per-dataset tables (mean ± SD at every budget) and ARI/NMI curves are in the manuscript and in the repository's output folders.

## Limitations

- **Sequencing depth.** Both methods assume the same expected total count in every cell. Neither models depth explicitly, and ReThiN's within-cell normalization is an approximation without a finite-depth error bound.
- **Choice of null.** The Poisson null worked best on UMI counts and the NB null on read counts, so the null must be chosen from the sequencing protocol. The read-count evidence rests on two datasets.
- **ReThiN on read counts.** ReThiN's zero-correlation guarantee requires Poisson noise, which read counts violate. It still ranked among the best methods on both read-count datasets, but this is not guaranteed in general.
- **Batch and donor effects.** Neither method separates biological variation from differences between batches or donors.
- **Evaluation scope.** Methods are assessed only through *k*-means clustering with the true number of populations (ARI/NMI), on seven datasets. Differences are descriptive and not statistically tested.

To reproduce the benchmark:

```bash
git clone https://github.com/MaitreyaGanu/PLit-ReThiN.git
cd PLit-ReThiN
```

```r
# Core dependencies (see the repository for the complete list)
install.packages(c("Matrix", "MASS", "dplyr", "ggplot2"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("SingleCellExperiment", "scuttle", "scran", "scry", "M3Drop"))
install.packages("Seurat")
```

Each dataset has a Poisson script (PLit (Poisson), ReThiN and the baselines) and an NB script (PLit (NB) and the baselines). Each script writes `benchmark_summary_<dataset>.csv`, `raw_results_<dataset>.csv` and `runtime_<dataset>.csv`.

## Metrics

- **ARI (Adjusted Rand Index):** agreement between predicted clusters and ground-truth labels, corrected for chance. 1 = perfect, about 0 = random.
- **NMI (Normalized Mutual Information):** information shared between clusters and labels, normalized to [0, 1]. Higher is better.

## ⭐ Star history

<a href="https://star-history.com/#MaitreyaGanu/PLit-ReThiN&Date">
  <img src="https://api.star-history.com/svg?repos=MaitreyaGanu/PLit-ReThiN&type=Date" alt="Star history chart" width="600"/>
</a>

---

<div align="center">
<sub><b>Status:</b> manuscript in preparation • results reproduce from the scripts in this repository.</sub>
</div>
