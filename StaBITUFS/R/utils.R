# ---- internal helpers (not exported) ----

# Coerce input to a genes x cells numeric count matrix and validate.
.prepare_counts <- function(X, genes_are_rows) {
  X <- as.matrix(X)
  if (!genes_are_rows) X <- t(X)
  storage.mode(X) <- "double"
  if (anyNA(X)) stop("X contains missing values.")
  if (any(X < 0)) stop("X must contain non-negative counts.")
  if (any(abs(X - round(X)) > 1e-8)) {
    warning("X does not look like integer counts; values were rounded.")
    X <- round(X)
  }
  X
}

# MLE of the Negative Binomial dispersion r for one gene, by a one-dimensional
# numerical search over log r on the profile log-likelihood (mu fixed at the
# sample mean). Returns a finite cap (~ Poisson limit) for genes with no
# estimable overdispersion.
.estimate_theta <- function(x, mu) {
  if (mu <= 0) return(1e6)
  n  <- length(x)
  sx <- sum(x)
  neg_profile_ll <- function(logr) {
    r  <- exp(logr)
    ll <- sum(lgamma(x + r)) - n * lgamma(r) +
          n * r * log(r / (r + mu)) + sx * log(mu / (r + mu))
    -ll                                   # optimize() minimises, so negate
  }
  opt <- optimize(neg_profile_ll, interval = c(log(1e-3), log(1e6)))
  min(exp(opt$minimum), 1e6)
}

# Per-row Pearson correlation between two matrices (row j vs row j across cols).
.row_corr <- function(A, B) {
  Am  <- A - rowMeans(A)
  Bm  <- B - rowMeans(B)
  num <- rowSums(Am * Bm)
  den <- sqrt(rowSums(Am^2) * rowSums(Bm^2))
  out <- num / den
  out[!is.finite(out)] <- 0        # zero-variance halves -> correlation 0
  out
}

# Assemble a best-first ranked data.frame of scores.
.rank_result <- function(score, feats, m) {
  if (is.null(feats)) feats <- paste0("feature_", seq_len(m))
  ord <- order(score, decreasing = TRUE)
  data.frame(
    feature = feats[ord],
    score   = score[ord],
    rank    = seq_len(m),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}
