# -----------------------------------------------------------------------------
# 1.  Packages
# -----------------------------------------------------------------------------

cat("Loading packages...\n")
suppressPackageStartupMessages({
  library(DuoClustering2018) 
  library(SingleCellExperiment)
  library(scran); library(scater); library(scry)
  library(Seurat); library(M3Drop)
  library(mclust); library(aricode)
  library(pbapply); library(progress)
  library(ggplot2); library(dplyr); library(tidyr)
})
pboptions(type = "timer")
options(ExperimentHub.ask = FALSE)

# -----------------------------------------------------------------------------
# 2.  Data  ( 8 purified populations mixed in known proportions)
# -----------------------------------------------------------------------------
cat("Loading Zhengmix8eq...\n")
sce         <- suppressMessages(sce_full_Zhengmix8eq())
counts_mat  <- as.matrix(counts(sce))                         # dense (downstream ops densify anyway)
true_labels <- factor(colData(sce)$phenoid)                
n_true_k    <- nlevels(true_labels)
print(table(true_labels))

counts_mat  <- counts_mat[rowSums(counts_mat > 0) >= 10, ]
cat(sprintf("\nDataset: %d genes x %d cells, %d cell types\n\n",
            nrow(counts_mat), ncol(counts_mat), n_true_k))

# Size factors computed ONCE on the full gene set, reused for every K-subset, so
# normalisation never depends on which genes a method selected (fairness control).
sce_full <- logNormCounts(SingleCellExperiment(assays = list(counts = counts_mat)))
full_sf  <- sizeFactors(sce_full)

# -----------------------------------------------------------------------------
# 3.  Experimental constants
# -----------------------------------------------------------------------------
TOP_K         <- c(100L, 200L, 500L, 1000L)
BS_SEEDS      <- 1:5        # bootstrap seeds -> error bars (applied to EVERY method)
B_ROUNDS      <- 20L        # bootstrap rounds per seed
KMEANS_SEEDS  <- 1:30       
KMEANS_NSTART <- 25L

# -----------------------------------------------------------------------------
# 4.  Core scorers   (counts -> named per-gene score; higher = more important)
# -----------------------------------------------------------------------------
row_cor <- function(A, B) {
  Ac <- A - rowMeans(A); Bc <- B - rowMeans(B)
  num <- rowSums(Ac * Bc); den <- sqrt(rowSums(Ac^2) * rowSums(Bc^2))
  ifelse(den == 0, 0, num / den)
}

score_ReThiN <- function(counts, n_thin = 5L) {       
  m <- nrow(counts); n <- ncol(counts); cv <- as.vector(counts)
  nz <- which(cv > 0); xnz <- cv[nz]
  rt <- matrix(0, m, n_thin)
  for (t in seq_len(n_thin)) {
    Av <- cv; Av[nz] <- rbinom(length(xnz), xnz, 0.5); Bv <- cv - Av
    Am <- matrix(Av, m); Bm <- matrix(Bv, m)
    sA <- colSums(Am); sA[sA == 0] <- 1; sB <- colSums(Bm); sB[sB == 0] <- 1
    rt[, t] <- row_cor(sweep(Am, 2, sA, "/"), sweep(Bm, 2, sB, "/"))
  }
  setNames(rowMeans(rt), rownames(counts))
}

score_PLit <- function(counts) {               
  n <- ncol(counts); lam <- rowMeans(counts)
  log_lam <- ifelse(lam > 0, log(lam), 0)
  L0 <- n * lam * log_lam - n * lam - rowSums(lgamma(counts + 1L))   # Poisson log-likelihood
  tab <- apply(counts, 1L, function(x) {
    ct <- tabulate(x + 1L); ct <- ct[ct > 0L]; c(sum(ct * log(ct / n)), length(ct))
  })
  Sj <- tab[1L, ] - L0 - ((tab[2L, ] - 2L) / 2) * log(n)            # = L0_code - L1_code
  Sj[lam == 0] <- -Inf
  setNames(Sj, rownames(counts))
}

score_pearson <- function(counts) {                           # Lause analytic Pearson residuals
  n <- ncol(counts)
  mu <- outer(rowSums(counts), colSums(counts)) / sum(counts)
  z  <- (counts - mu) / sqrt(mu + mu^2 / 100)
  cl <- sqrt(n); z[z > cl] <- cl; z[z < -cl] <- -cl
  setNames(apply(z, 1, var), rownames(counts))
}

score_scran <- function(counts, sf = NULL) {                  # use the injected (full) size factors, not re-estimated ones
  s <- SingleCellExperiment(assays = list(counts = counts))
  if (!is.null(sf)) sizeFactors(s) <- sf
  s <- logNormCounts(s); d <- modelGeneVar(s)
  setNames(d$bio, rownames(d))[rownames(counts)]
}

score_deviance <- function(counts)
  setNames(as.numeric(devianceFeatureSelection(counts)), rownames(counts))

score_seurat <- function(counts) {
  options(Seurat.object.assay.version = "v3")
  so <- suppressWarnings(suppressMessages(CreateSeuratObject(counts = counts)))
  so <- suppressMessages(FindVariableFeatures(so, selection.method = "vst",
                                              nfeatures = nrow(counts), verbose = FALSE))
  hv <- HVFInfo(so)
  setNames(hv$variance.standardized[match(gsub("_", "-", rownames(counts)), rownames(hv))],
           rownames(counts))
}

score_m3drop <- function(counts) {
  m3n <- suppressMessages(M3DropConvertData(counts, is.counts = TRUE))
  m3f <- suppressMessages(M3DropFeatureSelection(m3n, mt_method = "fdr",
                                                 mt_threshold = 1, suppress.plot = TRUE))
  sc <- setNames(numeric(nrow(counts)), rownames(counts))
  hit <- intersect(m3f$Gene, names(sc))
  sc[hit] <- -log10(pmax(m3f$p.value[match(hit, m3f$Gene)], 1e-300)); sc
}

score_random <- function(counts) setNames(runif(nrow(counts)), rownames(counts))

SCORERS <- list(
  ReThiN = score_ReThiN, PLit = score_PLit, scran_HVG = score_scran,
  scry_Deviance = score_deviance, Seurat_VST = score_seurat,
  Pearson_Residuals = score_pearson, M3Drop = score_m3drop, Random_Baseline = score_random
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
# 5.  Run every selector (subsampling aggregation)
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
  sizeFactors(sub) <- full_sf
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
write.csv(raw_df, "raw_results_ZhengPBMC.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 7.  Aggregation  (avg over k-means seeds; then mean +/- SD over bootstrap seeds)
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
cat("\n================== BENCHMARK SUMMARY ==================\n")
print(as.data.frame(summary_df), row.names = FALSE)
write.csv(summary_df, "benchmark_summary_ZhengPBMC.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 8.  Runtime summary
# -----------------------------------------------------------------------------
runtime_summary <- runtime_df |>
  dplyr::group_by(Method) |>
  dplyr::summarise(Runtime_mean_s = round(mean(Runtime_s), 2),
                   Runtime_sd_s   = round(sd(Runtime_s), 2), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(Runtime_mean_s))
cat("\n================== RUNTIME (seconds) ==================\n")
print(as.data.frame(runtime_summary), row.names = FALSE)
write.csv(runtime_summary, "runtime_ZhengPBMC.csv", row.names = FALSE)
      
# -----------------------------------------------------------------------------
# 9.  Figures
# -----------------------------------------------------------------------------
p_ari <- ggplot(summary_df, aes(K, ARI, colour = Method, group = Method)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = ARI - ARI_sd, ymax = ARI + ARI_sd), width = 25, alpha = 0.5) +
  scale_x_continuous(breaks = TOP_K) + theme_bw(base_size = 13) +
  labs(title = "Feature selection - ARI (Zheng 2017 FACS PBMC)",
       subtitle = "Error bars: bootstrap-seed SD (all methods bootstrapped identically)",
       x = "Features selected (K)", y = "ARI (mean +/- SD)", colour = NULL)
p_nmi <- ggplot(summary_df, aes(K, NMI, colour = Method, group = Method)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = NMI - NMI_sd, ymax = NMI + NMI_sd), width = 25, alpha = 0.5) +
  scale_x_continuous(breaks = TOP_K) + theme_bw(base_size = 13) +
  labs(title = "Feature selection - NMI (Zheng 2017 FACS PBMC)",
       subtitle = "Error bars: bootstrap-seed SD (all methods bootstrapped identically)",
       x = "Features selected (K)", y = "NMI (mean +/- SD)", colour = NULL)
p_rt <- runtime_summary |>
  dplyr::mutate(Method = reorder(Method, Runtime_mean_s)) |>
  ggplot(aes(Method, Runtime_mean_s)) +
  geom_col(fill = "steelblue", width = 0.6) +
  geom_errorbar(aes(ymin = pmax(Runtime_mean_s - Runtime_sd_s, 0),
                    ymax = Runtime_mean_s + Runtime_sd_s), width = 0.3) +
  coord_flip() + theme_bw(base_size = 13) +
  labs(title = "Wall-clock runtime - Zheng 2017 FACS PBMC", x = NULL, y = "Time (seconds)")
ggsave("fig_ARI_ZhengPBMC.pdf",     p_ari, width = 8, height = 5)
ggsave("fig_NMI_ZhengPBMC.pdf",     p_nmi, width = 8, height = 5)
ggsave("fig_Runtime_ZhengPBMC.pdf", p_rt,  width = 6, height = 5)

# -----------------------------------------------------------------------------
# 10.  Session info
# -----------------------------------------------------------------------------
cat("\n================== SESSION INFO ==================\n")
si <- sessionInfo(); print(si)
sink("session_info_ZhengPBMC.txt"); print(si); sink()
