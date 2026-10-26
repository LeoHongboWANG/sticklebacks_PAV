#!/usr/bin/env Rscript
# 正式跑：4 条链并行（每核一条）+ 收敛诊断 + 合并后验
# Usage: Rscript scripts/phylogenetic_model.R data results/phylogenetic_model [--validate-only]
suppressPackageStartupMessages({
  library(MCMCglmm); library(ape); library(coda); library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args) >= 1) args[1] else "data", mustWork = TRUE)
outdir <- if (length(args) >= 2) args[2] else "results/phylogenetic_model"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
outdir <- normalizePath(outdir, mustWork = TRUE)
setwd(root)

d  <- read.csv("analysis_outputs/excov_depth_robustness/sample_losttxt_with_coverage.csv",
               stringsAsFactors = FALSE)
tr <- read.tree("p_dis_2.nwk")

d <- d[d$SampleID %in% tr$tip.label, ]
d$animal         <- d$SampleID
d$Ecotype        <- factor(d$EcotypeSimple, levels = c("Marine", "Freshwater"))  # Marine = 参照
d$Log10NePiPlus1 <- as.numeric(scale(log10(d$NePi + 1)))
d$Lost           <- as.integer(round(d$LostTxt))
cat(sprintf("n = %d   freshwater = %d   marine = %d\n",
            nrow(d), sum(d$Ecotype == "Freshwater"), sum(d$Ecotype == "Marine")))

if (!is.ultrametric(tr)) tr <- chronos(tr, quiet = TRUE)
tr <- multi2di(tr, random = FALSE)
tr$edge.length[tr$edge.length <= 0] <- 1e-8
Ainv <- inverseA(tr, nodes = "TIPS", scale = TRUE)$Ainv

prior <- list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1, nu = 0.002)))

nitt <- 1300000; burn <- 300000; thin <- 500; nch <- 4

if ("--validate-only" %in% args) {
  cat(sprintf("Validation passed: %d samples; %d tree tips; NePi scaled; Marine reference\n", nrow(d), length(tr$tip.label)))
  quit(status = 0)
}
cat(sprintf("launching %d chains in parallel (nitt=%d, burnin=%d, thin=%d)\n",
            nch, nitt, burn, thin))
t0 <- Sys.time()

chains <- mclapply(seq_len(nch), function(i) {
  set.seed(1000 + i)
  MCMCglmm(Lost ~ Log10NePiPlus1 + Ecotype,
           random = ~animal, ginverse = list(animal = Ainv),
           family = "poisson", data = d, prior = prior,
           nitt = nitt, burnin = burn, thin = thin, verbose = FALSE)
}, mc.cores = nch, mc.set.seed = FALSE)

cat(sprintf("all chains done in %.1f h\n",
            as.numeric(difftime(Sys.time(), t0, units = "hours"))))

ok <- vapply(chains, function(x) inherits(x, "MCMCglmm"), logical(1))
if (!all(ok)) { cat("FAILED chains:", which(!ok), "\n"); print(chains[!ok]); quit(status = 1) }

sols   <- lapply(chains, function(m) as.mcmc(m$Sol))
ml     <- mcmc.list(sols)
pooled <- do.call(rbind, lapply(sols, as.matrix))

cat("\n===== Gelman-Rubin (potential scale reduction factor) =====\n")
print(gelman.diag(ml, autoburnin = FALSE, multivariate = FALSE))
cat("\n===== Effective sample size (pooled) =====\n")
print(round(effectiveSize(as.mcmc(pooled)), 1))

cat("\n===== Pooled posterior =====\n")
for (v in colnames(pooled)) {
  x  <- pooled[, v]; ci <- quantile(x, c(.025, .975))
  pm <- 2 * min(mean(x > 0), mean(x < 0))
  cat(sprintf("%-22s beta = %+8.4f  [%+8.4f, %+8.4f]  IRR = %6.3f [%.3f, %.3f]  pMCMC = %s\n",
              v, mean(x), ci[1], ci[2], exp(mean(x)), exp(ci[1]), exp(ci[2]),
              ifelse(pm < 1/nrow(pooled), paste0("<", signif(1/nrow(pooled), 2)),
                     signif(max(pm, 1/nrow(pooled)), 3))))
}

cat("\n===== Per-chain posterior means (range across chains) =====\n")
pc <- sapply(sols, function(s) colMeans(as.matrix(s)))
for (v in rownames(pc))
  cat(sprintf("%-22s %+8.4f .. %+8.4f\n", v, min(pc[v, ]), max(pc[v, ])))

saveRDS(chains, file.path(outdir, "mcmcglmm_full_chains.rds"))
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
