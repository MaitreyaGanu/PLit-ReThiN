#' PLit: Parametric Length Information Test
#'
#' Ranks genes by a two-part minimum-description-length (MDL) code-length gap
#' between a fitted parametric null and the empirical (saturated) distribution
#' of each gene's counts. Before the complexity penalty, the score equals
#' \code{n * KL(empirical || fitted null)}.
#'
#' The penalty uses \code{d0} null parameters: \code{d0 = 1} for the Poisson
#' (penalty \code{(V_j - 2)/2 * log n}) and \code{d0 = 2} for the Negative
#' Binomial (penalty \code{(V_j - 3)/2 * log n}), where \code{V_j} is the number
#' of distinct observed count values for gene \code{j}.
#'
#' @param X Raw count matrix; genes in rows and cells in columns by default.
#' @param model Null model: \code{"poisson"} (\code{d0 = 1}) or \code{"nb"}
#'   (Negative Binomial, \code{d0 = 2}).
#' @param genes_are_rows Logical; set \code{FALSE} if your matrix is cells x genes.
#' @return A \code{data.frame} with columns \code{feature}, \code{score} and
#'   \code{rank}, ordered best-first (highest score = rank 1).
#' @examples
#' \dontrun{
#'   res <- PLit(X, model = "poisson")
#'   top <- head(res$feature, 500)
#' }
#' @export
PLit <- function(X, model = c("poisson", "nb"), genes_are_rows = TRUE) {
  model <- match.arg(model)
  X  <- .prepare_counts(X, genes_are_rows)
  m  <- nrow(X); n <- ncol(X)
  d0 <- if (model == "poisson") 1L else 2L
  logn <- log(n)

  score <- numeric(m)
  for (j in seq_len(m)) {
    x  <- X[j, ]
    ct <- as.numeric(table(x))          # frequencies of distinct count values
    Vj <- length(ct)

    if (Vj <= d0) {                      # too few distinct values to score
      score[j] <- -Inf
      next
    }

    L1 <- sum(ct * log(ct / n))          # empirical (saturated) log-likelihood

    if (model == "poisson") {
      lam <- mean(x)
      L0  <- sum(stats::dpois(x, lambda = lam, log = TRUE))
    } else {
      mu <- mean(x)
      r  <- .estimate_theta(x, mu)
      L0 <- sum(stats::dnbinom(x, size = r, mu = mu, log = TRUE))
    }

    penalty  <- ((Vj - 1) - d0) / 2 * logn
    score[j] <- L1 - L0 - penalty        # = n * KL(empirical || null) - penalty
  }
  .rank_result(score, rownames(X), m)
}
