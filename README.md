<div align="center">

# 🧬 PLit & ReThiN

### Stability-Based and Information-Theoretic Unsupervised Feature Selection Methods for Single-Cell RNA Sequencing

**Maitreya Sameer Ganu**
*Indian Institute of Science Education and Research (IISER), Thiruvananthapuram*
*Advisor: Dr. Clint P. George — Indian Institute of Technology (IIT) Goa*

<br>

![Language](https://img.shields.io/badge/Language-R-276DC3?style=for-the-badge&logo=r&logoColor=white)
![Field](https://img.shields.io/badge/Field-Bioinformatics-00758F?style=for-the-badge)
![Topic](https://img.shields.io/badge/Topic-Feature%20Selection-orange?style=for-the-badge)
![Models](https://img.shields.io/badge/Models-Poisson%20%7C%20Negative%20Binomial-2E8B57?style=for-the-badge)
![Status](https://img.shields.io/badge/Status-Manuscript%20in%20Preparation-yellow?style=for-the-badge)

<br>

[Abstract](#project-abstract) • [Key Findings](#key-findings-honest-summary) • [Core Hypothesis](#core-hypothesis) • [Methodology](#methodology) • [Theory](#theoretical-results) • [Results](#results) • [Limitations](#limitations) • [Usage](#usage) • [Citation](#citation)

</div>

---

## Project Abstract

Feature selection is a critical preprocessing step in single-cell RNA sequencing (scRNA-seq), shaping downstream clustering, cell-type annotation, and every analysis built on top of it. This project introduces two unsupervised feature selection methods for count data:

- **PLit** (*Parametric Length Information Test*) — ranks genes by comparing the description length of their empirical count distribution against a fitted parametric null, using the Minimum Description Length (MDL) principle.
- **ReThiN** (*Reproducibility via Thinning*) — ranks genes by how reproducible their expression profile is under data thinning, measured via a split-half correlation.

Both methods are implemented for the **Poisson** and **Negative Binomial (NB)** count models — PLit extends to any parametric count family, and ReThiN to any convolution-closed count family — and both are deployed inside a single subsampling-based stability wrapper. We benchmark all four instances (PLit-Poisson, PLit-NB, ReThiN-Poisson, ReThiN-NB) against five established feature selectors (scran HVG, Seurat VST, Pearson Residuals, M3Drop, scry Deviance) and a random baseline, across **seven public scRNA-seq datasets** with ground-truth cell labels, scoring downstream *k*-means clustering with Adjusted Rand Index (ARI) and Normalized Mutual Information (NMI) across four feature budgets (K = 100, 200, 500, 1000).

In their **Poisson** instance, both methods are highly competitive with state-of-the-art variance-based approaches: ReThiN attains the **best average rank of all seven methods**, PLit ranks third, and a proposed method lands in the **top three on ARI for every dataset** (best score across budgets). The **Negative Binomial** instance is a clean theoretical extension but an empirically *mixed* one — it helps PLit on some datasets (notably Segerstolpe and Darmanis) yet does not uniformly improve over the Poisson version, and ReThiN-NB in particular underperforms (see [Key Findings](#key-findings-honest-summary)). Both methods use **at most one hyperparameter** (none for PLit, one for ReThiN) — fewer than the trend- and residual-based baselines (scran, Seurat, Pearson) and on par with the parameter-free ones (scry, M3Drop).

## Key Findings (Honest Summary)

> This section states the results as plainly as possible, including where the methods **do not** win. The benchmark uses only 7 datasets, so most differences are directional, not statistically certified.

**What works (Poisson instances):**
- **ReThiN (Poisson) has the best mean rank of all 7 methods** on both ARI (2.57) and NMI (2.71), ahead of every established baseline.
- **PLit (Poisson) ranks 3rd on ARI** (3.57), ahead of four of the five baselines.
- A proposed method appears in the **top-3 ARI on all 7 datasets**, and ReThiN (Poisson) has the **single best NMI** on 3 of 7 datasets (Segerstolpe, Darmanis, Zeisel).

**What is only weakly supported:**
- With 7 datasets, dataset-level bootstrap CIs exclude zero on **both** metrics for exactly **two** comparisons: **ReThiN > Seurat VST** and **ReThiN > Pearson Residuals**. Every other win (including PLit's) is a directional point estimate, not a certified effect.
- PLit shows a **small but statistically supported deficit vs. scry Deviance on NMI** — an honest negative for PLit.

**What does *not* work (Negative Binomial instances):**
- The NB extension is theoretically clean but **empirically a net negative**. Under the NB instance, **ReThiN falls to the worst mean rank of the seven methods** (6.57 ARI / 6.57 NMI) and **PLit falls to second-to-last** (4.71 / 4.86); the parameter-free M3Drop becomes the top-ranked method.
- NB helps only in specific places — chiefly **PLit-NB on Segerstolpe and Darmanis**. It is presented as an initial extension that needs refinement, **not** a replacement for the Poisson formulation.

**Bottom line:** the contribution is the pair of *Poisson* estimators, which are competitive-and-simple; the NB generalization is a correct derivation whose empirical payoff is dataset-dependent and, on average, unfavorable.

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
2. For 5 independent seeds, run 20 rounds of 80%-cell subsampling (**without replacement**); apply each method's core score per round and aggregate the 20 rankings by average rank.
3. For each feature budget K ∈ {100, 200, 500, 1000}, select the top-K genes, log-normalize, run PCA (15 PCs), and cluster with *k*-means (30 seeds × 25 restarts) using the true number of populations.
4. Report mean ± SD of ARI and NMI over the 5 subsampling seeds.

> **Note on terminology:** the wrapper uses **80% subsampling without replacement**, which is *subsampling* (not bootstrap resampling). It is inspired by stability selection (Meinshausen & Bühlmann, 2010; Shah & Samworth, 2013) but does **not** reproduce their complementary-pairs procedure, so their formal false-selection guarantees do **not** apply to these benchmark numbers. The stability layer here is a variance-reduction/fair-comparison device, not an error-control guarantee.

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

*Ground-truth quality varies: Tian CellBench (cell lines mixed by design) is experimentally fixed; Zhengmix4eq/8eq are kit-purified PBMCs computationally mixed in known proportions (not FACS-sorted); Baron, Segerstolpe, Darmanis, and Zeisel provide expert-curated marker-based annotations, so agreement on these should be read against a curated computational reference rather than an independently measured identity.*

## Theoretical Results

Both scores reduce to closed-form population quantities. (Renders as math on GitHub.)

**PLit** — penalized empirical KL divergence from the fitted null, where $V_j$ = number of distinct observed counts and $d_0$ = number of null parameters:

$$S_j = n\,\widehat{\mathrm{KL}}\!\left(\hat p \;\middle\|\; f(\cdot;\hat\theta_j)\right) - \frac{(V_j-1)-d_0}{2}\ln n$$

| Instance | Penalty term | Null parameters |
|---|---|---|
| Poisson ($d_0=1$) | $\tfrac{V_j-2}{2}\ln n$ | rate $\lambda_j$ |
| Negative Binomial ($d_0=2$) | $\tfrac{V_j-3}{2}\ln n$ | mean $\mu_j$, dispersion $r_j$ |

**ReThiN** — split-half correlation of two thinned halves, an estimate of a variance-components ratio ($\tilde\sigma_j^2$ = biological signal, $\bar w_j$ = mean thinning-noise floor):

$$\mathrm{Corr}(A_{ji}, B_{ji}) = \frac{\tilde\sigma_j^2}{\tilde\sigma_j^2 + 4\bar w_j}$$

| Instance | Closed form |
|---|---|
| Poisson | $\dfrac{\sigma_j^2}{\sigma_j^2 + 2\mu_j}$ |
| Negative Binomial | $\dfrac{\tilde\sigma_j^2}{\tilde\sigma_j^2 + 2\mu_j + \frac{2}{r_j}\,\mathbb{E}[\mu_{ji}^2]}$ |

The correlation is exactly **0** when a gene's variability is pure sampling noise ($\sigma_j^2 = 0$) and rises monotonically with the biological signal-to-noise ratio. The ReThiN identities are derived for unnormalized counts under an idealized uniform-depth assumption; in practice ReThiN normalizes within-cell, so the implemented statistic is an **approximation** to these population identities (see [Limitations](#limitations)).

## Results

For each dataset, the table below shows which methods achieve the best **ARI** (top 3) and best **NMI** (top 1), using each method's **best score across all four feature budgets** (K = 100, 200, 500, 1000) tested, excluding the random baseline. The four rightmost columns show where PLit and ReThiN — under their Poisson and NB instances — land in that ranking.

> ⚠️ **Read this as a "best-case per method" view.** Taking each method's best budget flatters *every* method equally (it is applied to baselines too), so it is a fair *relative* comparison but an optimistic *absolute* one. The [Aggregate Ranking](#aggregate-ranking-across-all-7-datasets) below, which averages *across* budgets, is the more conservative summary.

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

† *Tian CellBench saturates quickly: by K ≥ 200, nearly every method (proposed and baseline) reaches ARI ≈ 0.997 / NMI ≈ 0.993, so this dataset does not meaningfully discriminate between feature selectors. PLit-NB plateaus at 0.990 and is the one method outside the 8-way tie.*

### Aggregate Ranking Across All 7 Datasets

Averaging each method's ARI/NMI across feature budgets, then ranking the 7 methods per dataset (1 = best), and averaging ranks across all 7 datasets. **The Poisson and NB columns each rank the five baselines together with the corresponding instance of the proposed methods** — so this shows, honestly, how the two proposed instances fare against the field:

| Method | Poisson — ARI | Poisson — NMI | NB — ARI | NB — NMI |
|---|---|---|---|---|
| **ReThiN** † | **2.57** 🥇 | **2.71** 🥇 | 6.57 | 6.57 |
| M3Drop | 3.14 | 3.00 | **2.71** 🥇 | **2.57** 🥇 |
| **PLit** † | 3.57 | 3.79 | 4.71 | 4.86 |
| scry Deviance | 3.71 | 3.43 | 3.00 | 2.86 |
| scran HVG | 4.64 | 4.21 | 3.21 | 3.07 |
| Pearson Residuals | 4.79 | 4.57 | 3.36 | 3.36 |
| Seurat VST | 5.57 | 6.29 | 4.43 | 4.71 |

† Proposed method. *(Rank 1 = best; the seven ranks sum to 28 in each column.)*

**Reading it straight:** under the **Poisson** instance ReThiN is 1st and PLit is 3rd; under the **NB** instance the same two methods drop to **last (6.57)** and **second-to-last (4.71)**, and M3Drop takes the top spot. The Poisson methods are the result to take away; the NB instances are an honest negative.

> **A note on significance:** with only 7 benchmark datasets, pairwise significance tests are underpowered. Dataset-level bootstrap confidence intervals (10,000 resamples) show only **two** comparisons excluding zero on both metrics: **ReThiN over Seurat VST** and **ReThiN over Pearson Residuals**. PLit additionally shows a small but statistically supported **deficit vs. scry Deviance on NMI**. All other differences reported above — including PLit's and ReThiN's other wins — are directionally consistent point estimates, not certified effects at this sample size. See the manuscript for full CI tables.

## Limitations

- **Approximate normalization in ReThiN.** The split-half correlation identities assume unnormalized counts under uniform sequencing depth; the implementation normalizes within-cell, so the theory is an approximation and no finite-sample error bound is established yet.
- **The NB extension underperforms.** It is a clean derivation but, empirically, usually worse than the Poisson version (ReThiN-NB especially). It should be read as a proof-of-concept for the generalization, not a recommended default.
- **Benchmarking ≠ error-controlled deployment.** The 80% subsampling wrapper enables fair comparison but does not carry the formal false-selection guarantees of complementary-pairs stability selection.
- **Evaluation scope.** Methods are assessed only through downstream *k*-means clustering (ARI/NMI). Trajectory inference, differential expression, cell-type annotation, and deep-learning pipelines are untested.
- **Small benchmark.** Seven datasets limit statistical power; most cross-method differences are directional rather than certified.

## Metrics

- **ARI (Adjusted Rand Index):** agreement between predicted clusters and ground-truth labels, corrected for chance. 1 = perfect, ~0 = random.
- **NMI (Normalized Mutual Information):** shared information between clusters and labels, normalized to [0, 1]. Higher = better.

## Repository

All code, benchmarking scripts, and processed datasets used in this project are available at:

🔗 **https://github.com/MaitreyaGanu/PLit-ReThiN**

See the repository for the per-dataset entry points, the shared subsampling wrapper, and the exact package/version dependencies.

## Usage

```r
# Clone the repository
git clone https://github.com/MaitreyaGanu/PLit-ReThiN.git
cd PLit-ReThiN

# Core R dependencies (see the repo for the authoritative, complete list)
install.packages(c("Matrix", "MASS"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("scran", "scry", "M3Drop", "Seurat"))

# Run the benchmarking pipeline
# (see repository for dataset-specific entry points and the shared subsampling wrapper)
```
---

<div align="center">
<sub><b>Status:</b> manuscript in preparation • author affiliations to be finalized • results reproduce from the scripts in this repository.</sub>
</div>
