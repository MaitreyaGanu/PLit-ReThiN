# -----------------------------------------------------------------------------
# 1.  Packages
# -----------------------------------------------------------------------------
cat("Loading packages...\n")
suppressPackageStartupMessages({
  library(scRNAseq); library(SingleCellExperiment)
  library(scran); library(scater); library(scry)
  library(Seurat); library(M3Drop)
  library(mclust); library(aricode)
  library(pbapply); library(progress)
  library(ggplot2); library(dplyr); library(tidyr)
})
pboptions(type = "timer")

# -----------------------------------------------------------------------------
# 2.  Data
# -----------------------------------------------------------------------------
cat("Loading Baron Human Pancreas dataset...\n")
sce         <- BaronPancreasData("human")
counts_mat  <- as.matrix(counts(sce))                         # dense (all downstream ops densify anyway)
counts_mat  <- counts_mat[rowSums(counts_mat > 0) >= 10, ]
true_labels <- factor(sce$label)
n_true_k    <- nlevels(true_labels)
cat(sprintf("Dataset: %d genes x %d cells, %d cell types\n\n",
            nrow(counts_mat), ncol(counts_mat), n_true_k))

# Size factors computed ONCE on the full gene set and reused for every K-subset,
# so normalisation never depends on which genes a method selected (fairness control).
sce_full <- logNormCounts(SingleCellExperiment(assays = list(counts = counts_mat)))
full_sf  <- sizeFactors(sce_full)

# -----------------------------------------------------------------------------
# 3.  Experimental constants
# -----------------------------------------------------------------------------
TOP_K         <- c(100L, 200L, 500L, 1000L)
BS_SEEDS      <- 1:5        # bootstrap seeds -> error bars (applied to EVERY method)
B_ROUNDS      <- 20L        # bootstrap rounds per seed
KMEANS_SEEDS  <- 1:30       # k-means repetitions to average out clustering noise
KMEANS_NSTART <- 25L
# -----------------------------------------------------------------------------
# 4.  Core scorers   (counts -> named per-gene score; higher = more important)
#     None of them sees the labels.  No bootstrap inside -- the wrapper adds it.
# -----------------------------------------------------------------------------
row_cor <- function(A, B) {                                   # vectorised row-wise Pearson
  Ac <- A - rowMeans(A); Bc <- B - rowMeans(B)
  num <- rowSums(Ac * Bc); den <- sqrt(rowSums(Ac^2) * rowSums(Bc^2))
  ifelse(den == 0, 0, num / den)
}

# ---- NB helpers shared by PLit_NB and ReThiN_NB ----------------------------
# Mean-dispersion NB: Var(X) = mu + mu^2/r ; r -> Inf recovers Poisson.
# Profile log-likelihood in r only (mu held fixed at the sample mean), dropping
# the r-independent -lgamma(x+1) term; maximised on the log(r) scale for
# numerical stability across the wide dynamic range of r.
nb_profile_loglik <- function(log_r, x, mu) {
  r <- exp(log_r); n <- length(x)
  sum(lgamma(x + r)) - n * lgamma(r) + n * r * log(r / (r + mu)) + sum(x) * log(mu / (r + mu))
}
fit_nb_r <- function(x, mu) {                                 # 1-D ML estimate of r_hat_j at fixed mu_hat_j
  if (mu <= 0) return(1e6)                                    # degenerate (all-zero) gene -> Poisson limit
  exp(optimize(nb_profile_loglik, interval = c(log(1e-3), log(1e6)),
               x = x, mu = mu, maximum = TRUE)$maximum)
}

score_ReThiN_NB <- function(counts, n_thin = 5L) {
  m <- nrow(counts); n <- ncol(counts); cv <- as.vector(counts)
  mu    <- rowMeans(counts)
  r_hat <- vapply(seq_len(m), function(j) fit_nb_r(counts[j, ], mu[j]), numeric(1))  # gene-wise NB size r_j
  nz <- which(cv > 0); xnz <- cv[nz]                          # thin only nonzeros (BetaBin(0,.,.)=0)
  r_nz <- rep(r_hat, n)[nz]                                    # broadcast gene-level r to nonzero cells
  rt <- matrix(0, m, n_thin)
  for (t in seq_len(n_thin)) {
    p_nz <- rbeta(length(xnz), r_nz / 2, r_nz / 2)             # A|X ~ Beta-Binomial(X, r/2, r/2)
    Av <- cv; Av[nz] <- rbinom(length(xnz), xnz, p_nz); Bv <- cv - Av
    Am <- matrix(Av, m); Bm <- matrix(Bv, m)
    sA <- colSums(Am); sA[sA == 0] <- 1; sB <- colSums(Bm); sB[sB == 0] <- 1
    rt[, t] <- row_cor(sweep(Am, 2, sA, "/"), sweep(Bm, 2, sB, "/"))
  }
  setNames(rowMeans(rt), rownames(counts))
}

score_PLit_NB <- function(counts) {
  n <- ncol(counts); mu <- rowMeans(counts)
  tab <- apply(counts, 1L, function(x) {                              # empirical log-lik + #values
    ct <- tabulate(x + 1L); ct <- ct[ct > 0L]
    c(sum(ct * log(ct / n)), length(ct))                       # tilde L1_j, V_j (same empirical code as Poisson PLit)
  })
  r_hat <- vapply(seq_len(nrow(counts)), function(j) fit_nb_r(counts[j, ], mu[j]), numeric(1))
  L0 <- vapply(seq_len(nrow(counts)), function(j) {
    x <- counts[j, ]; m <- mu[j]
    if (m == 0) return(NA_real_)
    r <- r_hat[j]
    sum(lgamma(x + r) - lgamma(r) - lgamma(x + 1L)) +
      n * r * log(r / (r + m)) + sum(x) * log(m / (r + m))     # NB log-likelihood
  }, numeric(1))
  Sj <- tab[1L, ] - L0 - ((tab[2L, ] - 3L) / 2) * log(n)        # d0 = 2 (mu, r) -> penalty (V_j-3)/2 * log n
  Sj[mu == 0 | tab[2L, ] < 3] <- -Inf
  setNames(Sj, rownames(counts))
}

score_pearson <- function(counts) {                           # Lause analytic Pearson residuals
  n <- ncol(counts)
  mu <- outer(rowSums(counts), colSums(counts)) / sum(counts)
  z  <- (counts - mu) / sqrt(mu + mu^2 / 100)
  cl <- sqrt(n); z[z > cl] <- cl; z[z < -cl] <- -cl
  setNames(apply(z, 1, var), rownames(counts))
}

score_scran <- function(counts, sf = NULL) {                  
  s <- SingleCellExperiment(assays = list(counts = counts))
  if (!is.null(sf)) sizeFactors(s) <- sf
  s <- logNormCounts(s); d <- modelGeneVar(s)
  setNames(d$bio, rownames(d))[rownames(counts)]
}

score_deviance <- function(counts)                            # scry binomial deviance
  setNames(as.numeric(devianceFeatureSelection(counts)), rownames(counts))

score_seurat <- function(counts) {                            # Seurat VST
  options(Seurat.object.assay.version = "v3")
  so <- suppressWarnings(suppressMessages(CreateSeuratObject(counts = counts)))
  so <- suppressMessages(FindVariableFeatures(so, selection.method = "vst",
                                              nfeatures = nrow(counts), verbose = FALSE))
  hv <- HVFInfo(so)
  setNames(hv$variance.standardized[match(gsub("_", "-", rownames(counts)), rownames(hv))],
           rownames(counts))
}

score_m3drop <- function(counts) {                            # M3Drop dropout-based
  m3n <- suppressMessages(M3DropConvertData(counts, is.counts = TRUE))
  m3f <- suppressMessages(M3DropFeatureSelection(m3n, mt_method = "fdr",
                                                 mt_threshold = 1, suppress.plot = TRUE))
  sc <- setNames(numeric(nrow(counts)), rownames(counts))
  hit <- intersect(m3f$Gene, names(sc))
  sc[hit] <- -log10(pmax(m3f$p.value[match(hit, m3f$Gene)], 1e-300)); sc
}

score_random <- function(counts) setNames(runif(nrow(counts)), rownames(counts))

SCORERS <- list(
  ReThiN_NB         = score_ReThiN_NB,
  PLit_NB           = score_PLit_NB,
  scran_HVG         = score_scran,
  scry_Deviance     = score_deviance,
  Seurat_VST        = score_seurat,
  Pearson_Residuals = score_pearson,
  M3Drop            = score_m3drop,
  Random_Baseline   = score_random
)

bootstrap_rank <- function(score_fn, counts, B, seed, sf, needs_sf = FALSE, label = "") {
  set.seed(seed); m <- nrow(counts); n <- ncol(counts); acc <- numeric(m)
  boot_idx <- lapply(seq_len(B), function(b) sample(n, round(0.8*n))) 
  pb <- progress_bar$new(
    format = paste0("    ", label, " [:bar] :current/:total | :percent | ETA: :eta"),
    total = B, clear = FALSE, width = 80, force = TRUE)
  for (b in seq_len(B)) {
    idx <- boot_idx[[b]]; bc <- counts[, idx, drop = FALSE]; colnames(bc) <-  paste0("cell", seq_len(ncol(bc))) 
    sc <- tryCatch(
      if (needs_sf) score_fn(bc, sf[idx]) else score_fn(bc),        
      error = function(e) setNames(rep(NA_real_, m), rownames(counts)))
    acc <- acc + rank(-sc[rownames(counts)], ties.method = "average", na.last = TRUE)
    pb$tick()
  }
  rownames(counts)[order(acc)]
}

# -----------------------------------------------------------------------------
# 5.  Run every selector through the SAME bootstrap (runtime logged per call)
#     NOTE: PLit_NB and ReThiN_NB each fit a 1-D NB dispersion MLE per gene per
#     call (via optimize()), so this run is markedly slower than the Poisson
#     version -- this is inherent to the NB instance, not a wrapper change.
# -----------------------------------------------------------------------------
methods_ranked <- list(); runtime_log <- list()
for (mname in names(SCORERS)) {
  cat(sprintf("-- %s  (%d seeds x B=%d) --\n", mname, length(BS_SEEDS), B_ROUNDS))
  for (s in BS_SEEDS) {
    cat(sprintf("  Seed %d/%d:\n", s, max(BS_SEEDS)))
    t0 <- proc.time()[["elapsed"]]
    methods_ranked[[sprintf("%s_s%d", mname, s)]] <-
      bootstrap_rank(SCORERS[[mname]], counts_mat, B_ROUNDS, s,
                     sf = full_sf, needs_sf = (mname == "scran_HVG"), label = mname)
    runtime_log[[length(runtime_log) + 1L]] <- data.frame(
      Method = mname, Seed = s, Runtime_s = proc.time()[["elapsed"]] - t0,
      stringsAsFactors = FALSE)
  }
}
runtime_df <- do.call(rbind, runtime_log)
cat(sprintf("\nAll selectors done. %d (method x seed) rankings.\n\n", length(methods_ranked)))

# -----------------------------------------------------------------------------
# 6.  Evaluation  --  PCA -> k-means (true K) -> ARI + NMI
# -----------------------------------------------------------------------------
cat("-- Evaluation: PCA -> k-means --\n")
grid <- expand.grid(Method_key = names(methods_ranked), K = TOP_K, stringsAsFactors = FALSE)
cat(sprintf("  %d evaluation jobs x %d k-means seeds\n", nrow(grid), length(KMEANS_SEEDS)))

raw_results <- pblapply(seq_len(nrow(grid)), function(i) {
  mkey <- grid$Method_key[i]; k <- grid$K[i]
  method  <- sub("_s[0-9]+$", "", mkey)
  bs_seed <- as.integer(sub(".*_s", "", mkey))
  top_genes <- head(na.omit(methods_ranked[[mkey]]), k)
  if (length(top_genes) < 5L) return(NULL)
  
  sub <- SingleCellExperiment(assays = list(counts = counts_mat[top_genes, ]))
  sizeFactors(sub) <- full_sf                                  # injected (not re-estimated)
  sub <- logNormCounts(sub)
  sub <- suppressWarnings(runPCA(sub, ncomponents = min(15L, length(top_genes) - 1L),
                                 ntop = length(top_genes)))
  pc <- reducedDim(sub, "PCA")
  
  do.call(rbind, lapply(KMEANS_SEEDS, function(ks) {
    set.seed(ks)
    km <- kmeans(pc, centers = n_true_k, nstart = KMEANS_NSTART)
    cl <- factor(km$cluster)
    data.frame(Method = method, BS_seed = bs_seed, K = k, KM_seed = ks,
               ARI = mclust::adjustedRandIndex(true_labels, cl),
               NMI = aricode::NMI(true_labels, cl), stringsAsFactors = FALSE)
  }))
})
raw_df <- do.call(rbind, Filter(Negate(is.null), raw_results))
write.csv(raw_df, "raw_results_BaronPancreas_NB.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 7.  Aggregation  (Step 1: average over k-means seeds; Step 2: mean +/- SD over
#     bootstrap seeds -> EVERY method now has a proper SD)
# -----------------------------------------------------------------------------
per_bs <- raw_df |>
  dplyr::group_by(Method, BS_seed, K) |>
  dplyr::summarise(ARI_bs = mean(ARI), NMI_bs = mean(NMI), .groups = "drop")

summary_df <- per_bs |>
  dplyr::group_by(Method, K) |>
  dplyr::summarise(ARI = round(mean(ARI_bs), 4), ARI_sd = round(sd(ARI_bs), 4),
                   NMI = round(mean(NMI_bs), 4), NMI_sd = round(sd(NMI_bs), 4),
                   .groups = "drop") |>
  dplyr::arrange(K, dplyr::desc(ARI))

cat("\n================== BENCHMARK SUMMARY (NB instance) ==================\n")
print(as.data.frame(summary_df), row.names = FALSE)
write.csv(summary_df, "benchmark_summary_BaronPancreas_NB.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 8.  Runtime summary  (mean +/- SD over bootstrap seeds, for every method)
# -----------------------------------------------------------------------------
runtime_summary <- runtime_df |>
  dplyr::group_by(Method) |>
  dplyr::summarise(Runtime_mean_s = round(mean(Runtime_s), 2),
                   Runtime_sd_s   = round(sd(Runtime_s), 2), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(Runtime_mean_s))
cat("\n================== RUNTIME (seconds) ==================\n")
print(as.data.frame(runtime_summary), row.names = FALSE)
write.csv(runtime_summary, "runtime_BaronPancreas_NB.csv", row.names = FALSE)
      
# -----------------------------------------------------------------------------
# 9.  Figures
# -----------------------------------------------------------------------------
p_ari <- ggplot(summary_df, aes(K, ARI, colour = Method, group = Method)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = ARI - ARI_sd, ymax = ARI + ARI_sd), width = 25, alpha = 0.5) +
  scale_x_continuous(breaks = TOP_K) + theme_bw(base_size = 13) +
  labs(title = "Feature selection - ARI (Baron Human Pancreas, NB instance)",
       subtitle = "Error bars: bootstrap-seed SD (all methods bootstrapped identically)",
       x = "Features selected (K)", y = "ARI (mean +/- SD)", colour = NULL)

p_nmi <- ggplot(summary_df, aes(K, NMI, colour = Method, group = Method)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = NMI - NMI_sd, ymax = NMI + NMI_sd), width = 25, alpha = 0.5) +
  scale_x_continuous(breaks = TOP_K) + theme_bw(base_size = 13) +
  labs(title = "Feature selection - NMI (Baron Human Pancreas, NB instance)",
       subtitle = "Error bars: bootstrap-seed SD (all methods bootstrapped identically)",
       x = "Features selected (K)", y = "NMI (mean +/- SD)", colour = NULL)

p_rt <- runtime_summary |>
  dplyr::mutate(Method = reorder(Method, Runtime_mean_s)) |>
  ggplot(aes(Method, Runtime_mean_s)) +
  geom_col(fill = "steelblue", width = 0.6) +
  geom_errorbar(aes(ymin = pmax(Runtime_mean_s - Runtime_sd_s, 0),
                    ymax = Runtime_mean_s + Runtime_sd_s), width = 0.3) +
  coord_flip() + theme_bw(base_size = 13) +
  labs(title = "Wall-clock runtime - Baron Human Pancreas (NB instance)", x = NULL, y = "Time (seconds)")

ggsave("fig_ARI_BaronPancreas_NB.pdf",     p_ari, width = 8, height = 5)
ggsave("fig_NMI_BaronPancreas_NB.pdf",     p_nmi, width = 8, height = 5)
ggsave("fig_Runtime_BaronPancreas_NB.pdf", p_rt,  width = 6, height = 5)

# -----------------------------------------------------------------------------
# 10.  Session info  (reproducibility)
# -----------------------------------------------------------------------------
cat("\n================== SESSION INFO ==================\n")
si <- sessionInfo(); print(si)
sink("session_info_BaronPancreas_NB.txt"); print(si); sink()
