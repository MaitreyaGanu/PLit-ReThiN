#' ReThiN: Reproducibility via Thinning
#'
#' Ranks genes by the split-half correlation of two independently thinned halves
#' of the count matrix. Each count is split with a fair thinning operator
#' (Binomial(x, 0.5) for the Poisson model; a Beta-Binomial split of the size
#' parameter for the Negative Binomial model), the two halves are normalised
#' within each cell, and their correlation across cells is computed per gene and
#' averaged over \code{n_thin} repetitions. High correlation indicates
#' reproducible biological structure; correlation near zero indicates pure
#' sampling noise.
#'
#' @param X Raw count matrix; genes in rows and cells in columns by default.
#' @param model Thinning model: \code{"poisson"} (fair binomial split) or
#'   \code{"nb"} (Beta-Binomial split with per-gene dispersion).
#' @param n_thin Number of thinning repetitions to average over (default 5).
#' @param genes_are_rows Logical; set \code{FALSE} if your matrix is cells x genes.
#' @return A \code{data.frame} with columns \code{feature}, \code{score} and
#'   \code{rank}, ordered best-first (highest mean split-half correlation = rank 1).
#' @examples
#' \dontrun{
#'   set.seed(1)
#'   res <- ReThiN(X, model = "poisson", n_thin = 5)
#'   top <- head(res$feature, 500)
#' }
#' @export
ReThiN <- function(X, model = c("poisson", "nb"), n_thin = 5, genes_are_rows = TRUE) {
  model <- match.arg(model)
  X <- .prepare_counts(X, genes_are_rows)
  m <- nrow(X); n <- ncol(X)

  r_gene <- NULL
  if (model == "nb") {                      # per-gene dispersion (size)
    r_gene <- vapply(seq_len(m),
                     function(j) .estimate_theta(X[j, ], mean(X[j, ])),
                     numeric(1))
  }

  rho_sum <- numeric(m)
  for (t in seq_len(n_thin)) {
    if (model == "poisson") {
      A <- matrix(stats::rbinom(length(X), size = as.integer(X), prob = 0.5),
                  nrow = m)
    } else {
      A <- matrix(0, nrow = m, ncol = n)
      for (j in seq_len(m)) {
        p <- stats::rbeta(n, r_gene[j] / 2, r_gene[j] / 2)   # one p per cell
        A[j, ] <- stats::rbinom(n, size = as.integer(X[j, ]), prob = p)
      }
    }
    B <- X - A

    aCol <- pmax(colSums(A), 1)             # per-cell library size of each half
    bCol <- pmax(colSums(B), 1)
    rA <- sweep(A, 2, aCol, "/")            # within-cell normalisation
    rB <- sweep(B, 2, bCol, "/")

    rho_sum <- rho_sum + .row_corr(rA, rB)
  }
  .rank_result(rho_sum / n_thin, rownames(X), m)
}
