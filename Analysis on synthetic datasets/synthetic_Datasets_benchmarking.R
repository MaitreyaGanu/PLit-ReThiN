suppressPackageStartupMessages({ library(StaBITUFS); library(aricode); library(pheatmap) })

OUT <- path.expand("~/Desktop/synthetic_benchmark")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

G          <- 10000        # genes
N_TYPES    <- 20           # cell types
PER_TYPE   <- 60           # cells per type  -> n = 1200 cells
N_IMP      <- 1000         # important genes per dataset (fixed for fair comparison)
KS         <- c(100, 200, 500, 1000)   # feature budgets (as in the real study)
N_THIN     <- 5            # ReThiN thinning repetitions
N_KM_SEEDS <- 30          
R_OVER     <- 1.5          # NB dispersion for the "nb" noise datasets (smaller = more overdispersed)
INSTANCES  <- c("PLit_Pois", "PLit_NB", "ReThiN_Pois", "ReThiN_NB")

USE_WRAPPER <- FALSE       # TRUE = full 5 seeds x 20 subsampling rounds (hours). FALSE = score once.
WRAP_SEEDS  <- 5
WRAP_ROUNDS <- 20

## ---- one synthetic dataset (20 types, a chosen signal + noise model) --------
gen <- function(signal, noise, seed = 1) {
  set.seed(seed)
  n    <- N_TYPES * PER_TYPE
  type <- rep(seq_len(N_TYPES), each = PER_TYPE)
  is_imp <- c(rep(TRUE, N_IMP), rep(FALSE, G - N_IMP))
  
  Mu   <- matrix(runif(G, 1, 5), G, n)   # null genes: flat baseline (row-constant)
  keep <- matrix(1L, G, n)               # dropout mask (1 = kept); only used for "dropout"
  for (g in which(is_imp)) {
    if (signal == "marker") {            # ON in one type, low elsewhere (strong)
      m <- sample(N_TYPES, 1); Mu[g, ] <- ifelse(type == m, runif(1, 15, 30), runif(1, 0.5, 1.5))
    } else if (signal == "subtle") {     # faint shift in one type (low SNR)
      m <- sample(N_TYPES, 1); b <- runif(1, 4, 9); Mu[g, ] <- ifelse(type == m, b * 1.4, b)
    } else if (signal == "multi") {      # elevated across a handful of types (broad)
      ms <- sample(N_TYPES, sample(3:6, 1)); Mu[g, ] <- ifelse(type %in% ms, runif(1, 8, 16), runif(1, 1, 3))
    } else if (signal == "dropout") {    # expressed in one type, dropped (zeros) elsewhere
      m <- sample(N_TYPES, 1); Mu[g, ] <- runif(1, 6, 12)
      keep[g, ] <- rbinom(n, 1, ifelse(type == m, 0.9, 0.15))
    }
  }
  drw <- if (noise == "poisson") function(mu) rpois(length(mu), mu)
  else                    function(mu) rnbinom(length(mu), size = R_OVER, mu = mu)
  X <- matrix(drw(as.vector(Mu)), G, n) * keep
  storage.mode(X) <- "integer"; rownames(X) <- paste0("g", seq_len(G))
  
  ok <- rowSums(X > 0) >= 10             # same gene filter as the real protocol
  X <- X[ok, ]; is_imp <- is_imp[ok]
  list(name = paste0(signal, "_", noise), signal = signal, noise = noise,
       X = X, type = type, important = rownames(X)[is_imp], n_imp = sum(is_imp))
}

## ---- feature selection: score once, or the full subsampling wrapper --------
score_one <- function(X, m) switch(m,
                                   PLit_Pois   = PLit(X, "poisson"),
                                   PLit_NB     = PLit(X, "nb"),
                                   ReThiN_Pois = ReThiN(X, "poisson", n_thin = N_THIN),
                                   ReThiN_NB   = ReThiN(X, "nb",      n_thin = N_THIN))

consensus_ranking <- function(X, m) {
  if (!USE_WRAPPER) return(score_one(X, m)$feature)
  acc <- setNames(numeric(nrow(X)), rownames(X)); reps <- 0
  for (s in seq_len(WRAP_SEEDS)) { set.seed(s)
    for (r in seq_len(WRAP_ROUNDS)) {
      cells <- sample(ncol(X), floor(0.8 * ncol(X)))
      res <- score_one(X[, cells, drop = FALSE], m)
      rk  <- setNames(res$rank, res$feature)
      acc <- acc + rk[names(acc)]; reps <- reps + 1
    }
  }
  names(sort(acc / reps))          # ascending mean rank = best first
}

## ---- clustering evaluation on a selected gene set (uses aricode) -----------
## Matches the real protocol's clustering step: k-means with 25 internal restarts,
## AVERAGED over N_KM_SEEDS independent seeds (so the ARI/NMI are reproducible and
## do not depend on a single random start).
cluster_eval <- function(X, genes, type, sf, n_seeds = N_KM_SEEDS) {
  sub  <- X[genes, , drop = FALSE]
  logn <- log1p(sweep(sub, 2, sf, "/") * median(sf))     # library-size normalise + log
  npc  <- min(15, nrow(logn) - 1, ncol(logn) - 1)
  pc   <- prcomp(t(logn), rank. = npc)$x
  a <- numeric(n_seeds); m <- numeric(n_seeds)
  for (s in seq_len(n_seeds)) {
    set.seed(s)
    km   <- kmeans(pc, centers = length(unique(type)), nstart = 25, iter.max = 50)$cluster
    a[s] <- ARI(type, km); m[s] <- NMI(type, km)
  }
  c(ARI = mean(a), NMI = mean(m))
}

## ---- build the 8 datasets (noise x signal grid) ----------------------------
grid <- expand.grid(signal = c("marker", "subtle", "multi", "dropout"),
                    noise  = c("poisson", "nb"), stringsAsFactors = FALSE)
datasets <- Map(function(sig, noi, sd) gen(sig, noi, seed = sd),
                grid$signal, grid$noise, seq_len(nrow(grid)))

## ---- run the benchmark -----------------------------------------------------
metric_rows <- list(); recov_rows <- list()
for (d in datasets) {
  message(sprintf("== %s : %d important of %d genes ==", d$name, d$n_imp, nrow(d$X)))
  sf     <- pmax(colSums(d$X), 1)
  ranked <- setNames(lapply(INSTANCES, function(m) consensus_ranking(d$X, m)), INSTANCES)
  
  for (m in INSTANCES) {
    r <- ranked[[m]]
    recov_rows[[length(recov_rows) + 1]] <- data.frame(
      dataset = d$name, signal = d$signal, noise = d$noise, instance = m,
      recovered_pct = round(100 * mean(d$important %in% r[seq_len(d$n_imp)]), 1))
    for (K in KS) {
      ce <- cluster_eval(d$X, r[seq_len(K)], d$type, sf)
      metric_rows[[length(metric_rows) + 1]] <- data.frame(
        dataset = d$name, signal = d$signal, noise = d$noise, instance = m, K = K,
        ARI = round(ce["ARI"], 3), NMI = round(ce["NMI"], 3))
    }
  }
}
metric_df <- do.call(rbind, metric_rows); rownames(metric_df) <- NULL
recov_df  <- do.call(rbind, recov_rows)
write.csv(metric_df, file.path(OUT, "clustering_metrics.csv"), row.names = FALSE)
write.csv(recov_df,  file.path(OUT, "recovery.csv"),           row.names = FALSE)

## ---- scoreboard: mean ARI over K, winner per dataset, heatmap --------------
agg  <- aggregate(cbind(ARI, NMI) ~ dataset + instance, metric_df, mean)
write.csv(agg, file.path(OUT, "mean_metrics.csv"), row.names = FALSE)
wide <- reshape(agg[, c("dataset", "instance", "ARI")],
                idvar = "dataset", timevar = "instance", direction = "wide")
rownames(wide) <- wide$dataset; wide$dataset <- NULL
colnames(wide) <- sub("ARI.", "", colnames(wide), fixed = TRUE)
wide <- as.matrix(wide)
winner <- data.frame(dataset = rownames(wide),
                     best_instance = colnames(wide)[max.col(wide, ties.method = "first")],
                     best_ARI = round(apply(wide, 1, max), 3))
write.csv(winner, file.path(OUT, "winner_per_dataset.csv"), row.names = FALSE)
pheatmap(wide, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE,
         number_format = "%.3f", main = "Mean ARI (dataset x instance) - who wins where",
         filename = file.path(OUT, "scoreboard_ARI.pdf"), width = 7, height = 5.5)

cat("\n=== WINNER PER DATASET (by mean ARI over K) ===\n"); print(winner, row.names = FALSE)
cat("\n=== RECOVERY (% of important genes) ===\n"); print(recov_df, row.names = FALSE)
cat("\nAll outputs saved to:", OUT, "\n")
