# StaBITUFS

**Sta**bility-**B**ased **I**nformation-**T**heoretic **U**nsupervised **F**eature **S**election for single-cell RNA sequencing.

Two assumption-light feature selectors that operate directly on raw count matrices:

- **`PLit()`** — *Parametric Length Information Test.* Ranks genes by a two-part minimum-description-length code-length gap against a fitted parametric null. Before the penalty, the score equals `n × KL(empirical ‖ fitted null)`.
- **`ReThiN()`** — *Reproducibility via Thinning.* Ranks genes by the split-half correlation of two independently thinned halves of the counts.

Both support a **Poisson** and a **Negative Binomial** model.

## Installation

```r
install.packages("remotes")
remotes::install_github("MaitreyaGanu/PLit-ReThiN", subdir = "StaBITUFS")
library(StaBITUFS)
```

## Usage

`X` is a raw count matrix with **genes in rows and cells in columns** (pass `genes_are_rows = FALSE` if yours is cells × genes).

```r
# PLit  (d0 = 1 for Poisson, d0 = 2 for Negative Binomial)
plit_pois <- PLit(X, model = "poisson")
plit_nb   <- PLit(X, model = "nb")
head(plit_pois)          # data.frame: feature, score, rank  (rank 1 = best)

# ReThiN
set.seed(1)              # ReThiN is stochastic; seed for reproducibility
rethin_pois <- ReThiN(X, model = "poisson", n_thin = 5)
rethin_nb   <- ReThiN(X, model = "nb", n_thin = 5)

# select the top-K genes
top_genes <- head(rethin_pois$feature, 500)
```

Each function returns a `data.frame` ordered best-first with columns `feature`, `score`, `rank`.

## Notes

- Input should be **raw integer counts**, not log-normalised values; non-integers are rounded with a warning.
- Filter very low-expression genes first (e.g. expressed in fewer than 10 cells), as in the paper.
- `ReThiN()` thins at random, so set a seed (`set.seed(...)`) for reproducible rankings.
- Both are `O(mn)` per pass; `ReThiN()` additionally scales with `n_thin`.

## Method & paper

These methods are introduced in *"Stability-Based and Information-Theoretic Unsupervised Feature Selection Methods for Single-Cell RNA Sequencing"* (Ganu & George). Full derivations and benchmarks: https://github.com/MaitreyaGanu/PLit-ReThiN

## License

MIT © 2026 Maitreya Sameer Ganu
