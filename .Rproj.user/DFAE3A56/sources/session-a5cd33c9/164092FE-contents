# ============================================================
# Longitudinal Multivariate Dataset Generator (Grouped Correlation Structure)
#
# N subjects, T time points, P variables, split into 3 blocks:
#   - Block 1 (size floor(P/3)): constant pairwise correlations
#         (each pair has its own user-specified fixed correlation)
#   - Block 2 (size floor(P/3)): cyclic pairwise correlations
#         (all pairs in this block share one sinusoidal cycle:
#          r(t) = base + amp * sin(2*pi*t/period + phase))
#   - Block 3 (remainder, P - 2*floor(P/3)): uncorrelated
#         (zero correlation with everything, including each other)
#
# Required packages: install.packages(c("dplyr"))
# ============================================================

library(tidyverse)
theme_set(theme_minimal())

# ---------- Nearest valid correlation matrix ----------
# Used as a safety net: if a block's specified correlations don't
# combine into a positive semi-definite matrix, project to the
# nearest valid correlation matrix (clip negative eigenvalues,
# rescale to unit diagonal).
nearest_correlation_matrix <- function(R) {
  eig <- eigen(R, symmetric = TRUE)
  vals <- pmax(eig$values, 1e-8)
  M <- eig$vectors %*% diag(vals) %*% t(eig$vectors)
  d <- sqrt(diag(M))
  C <- M / outer(d, d)
  diag(C) <- 1
  C
}

is_valid_corr <- function(R) {
  ch <- tryCatch(chol(R), error = function(e) NULL)
  !is.null(ch)
}

# ---------- Partition P variables into the three blocks ----------
# Block 1: constant-correlation block (size floor(P/3))
# Block 2: cyclic-correlation block   (size floor(P/3))
# Block 3: uncorrelated block         (remainder, including any leftover
#          from P not being divisible by 3)
partition_variables <- function(P) {
  block_size <- floor(P / 3)
  const_idx  <- if (block_size > 0) 1:block_size else integer(0)
  cyclic_idx <- if (block_size > 0) (block_size + 1):(2 * block_size) else integer(0)
  uncorr_idx <- setdiff(1:P, c(const_idx, cyclic_idx))
  list(const_idx = const_idx, cyclic_idx = cyclic_idx, uncorr_idx = uncorr_idx)
}

# ---------- Build the full P x P correlation matrix at time t ----------
# const_corr_matrix: a matrix of size length(const_idx) x length(const_idx)
#                     giving the user-specified fixed correlation for every
#                     pair in the constant block (diagonal = 1)
# cyclic_params: list(base, amp, period, phase) shared by ALL pairs in the
#                cyclic block
build_corr_matrix <- function(t, P, blocks, const_corr_matrix, cyclic_params) {
  R <- diag(P)

  # --- Constant block: pairwise correlations fixed at user-specified values ---
  ci <- blocks$const_idx
  if (length(ci) >= 2) {
    R[ci, ci] <- const_corr_matrix
  }

  # --- Cyclic block: all pairs share one sinusoidal correlation over time ---
  cy <- blocks$cyclic_idx
  if (length(cy) >= 2) {
    freq <- 2 * pi / cyclic_params$period
    r <- cyclic_params$base + cyclic_params$amp * sin(freq * t + cyclic_params$phase)
    r <- max(-0.999, min(0.999, r))
    R[cy, cy] <- r
    diag(R)[cy] <- 1
  }

  # --- correct correlation matrix: Precompute correlation matrix + Cholesky factor at every time point ---
  correction <- 0
  if (!is_valid_corr(R)) {
    correction <- 1
    # Correct block-by-block so uncorrelated / cross-block zeros aren't
    # disturbed by a global projection. Only the offending block(s) need it,
    # but a full-matrix nearest-PSD projection is simplest and still leaves
    # cross-block entries at (near) zero if the blocks are already fine.
    R <- nearest_correlation_matrix(R)
  }
  L <- t(chol(R))

  return(list(correction=correction, cor_mat = R, chol = L))
}

# ---------- Main dataset generator ----------
# Arguments:
#   n_subjects        - number of subjects (N)
#   n_time            - number of time points (T)
#   n_vars            - number of variables (P)
#   seed              - random seed for reproducibility
#   noise_sd          - global multiplier on noise scale
#   var_params        - list of length P; each element has
#                        intercept, trendSlope, seasonAmp, seasonPeriod, varSd
#   const_corr_values - named vector or matrix giving the fixed correlation
#                        for each pair in the constant block. Easiest way to
#                        supply this: a single number (applied to every pair
#                        in the block) or a full symmetric matrix.
#   cyclic_params     - list(base, amp, period, phase) for the cyclic block
#
# Returns: list(dataset, corr_by_time, n_corrections, blocks)
generate_dataset <- function(n_subjects, n_time, n_vars, seed, noise_sd,
                             var_params, const_corr_values, cyclic_params) {

  set.seed(seed)
  P <- n_vars
  var_names <- paste0("V", 1:P)
  blocks <- partition_variables(P)

  # --- Build the constant block's correlation matrix from user input ---
  n_const <- length(blocks$const_idx)
  if (n_const >= 2) {
    if (is.matrix(const_corr_values)) {
      stopifnot(nrow(const_corr_values) == n_const, ncol(const_corr_values) == n_const)
      const_corr_matrix <- const_corr_values
      diag(const_corr_matrix) <- 1
    } else if (length(const_corr_values) == 1) {
      # single value applied to every pair in the block
      const_corr_matrix <- matrix(const_corr_values, n_const, n_const)
      diag(const_corr_matrix) <- 1
    } else {
      stop("const_corr_values must be a single number or an n_const x n_const matrix, ",
           "where n_const = floor(P/3) = ", n_const)
    }
  } else {
    const_corr_matrix <- matrix(numeric(0), 0, 0)
  }

  # --- Precompute correlation matrix + Cholesky factor at every time point ---
  corr_by_time <- vector("list", n_time)
  chol_by_time <- vector("list", n_time)
  n_corrections <- 0

  for (t in 1:n_time) {
    R <- build_corr_matrix(t, P, blocks, const_corr_matrix, cyclic_params)

    n_corrections <- n_corrections+R$correction
    corr_by_time[[t]] <- R$cor_mat
    chol_by_time[[t]] <- R$chol
  }

  # --- Simulate subjects ---
  rows <- vector("list", n_subjects * n_time)
  i <- 1
  for (subj in 1:n_subjects) {
    for (t in 1:n_time) {
      z <- rnorm(P)
      L <- chol_by_time[[t]]
      correlated <- as.numeric(L %*% z)

      values <- numeric(P)
      for (v in 1:P) {
        vp <- var_params[[v]]
        trend  <- vp$intercept + vp$trendSlope * t
        season <- vp$seasonAmp * sin(2 * pi * t / vp$seasonPeriod)
        values[v] <- trend + season + correlated[v] * vp$varSd * noise_sd
      }

      row <- c(list(subject = subj, time = t), as.list(round(values, 4)))
      names(row)[3:(2 + P)] <- var_names
      rows[[i]] <- row
      i <- i + 1
    }
  }

  dataset <- bind_rows(rows)
  list(dataset = dataset, corr_by_time = corr_by_time,
       n_corrections = n_corrections, blocks = blocks, var_names = var_names)
}

# ---------- Main dataset generator ----------

P <- 9   # variables -> 3 constant, 3 cyclic, 3 uncorrelated
N <- 100  # subjects
T <- 100 # time points

# Default trend/seasonality parameters for each variable
# (feel free to customize per-variable; this just uses simple defaults)
var_params <- lapply(1:P, function(i) {
  list(
    intercept    = 0,
    trendSlope   = 0.02 * i * ifelse(i %% 2 == 0, 1, -1),
    seasonAmp    = 0,
    seasonPeriod = 1,
    varSd        = 1
  )
})

# Constant block: floor(9/3) = 3 variables -> 3 pairs.
# Supply a full 3x3 matrix so each pair gets its own value:
#   V1-V2 = 0.6, V1-V3 = 0.3, V2-V3 = 0.5
const_corr_matrix <- matrix(c(
  1.0, 0.6, 0.3,
  0.6, 1.0, 0.5,
  0.3, 0.5, 1.0
), nrow = 3, byrow = TRUE)

# Cyclic block: floor(9/3) = 3 variables, all pairs share one sinusoidal cycle
cyclic_params <- list(base = 0.2, amp = 0.5, period = 25, phase = 0)

result <- generate_dataset(
  n_subjects        = N,
  n_time            = T,
  n_vars            = P,
  seed              = 42,
  noise_sd          = 1,
  var_params        = var_params,
  const_corr_values = const_corr_matrix,
  cyclic_params     = cyclic_params
)

# Inspect
head(result$dataset)
dim(result$dataset)
result$corr_by_time

cat("Variable blocks:\n")
cat("  Constant-correlation block: ", paste(result$var_names[result$blocks$const_idx], collapse = ", "), "\n")
cat("  Cyclic-correlation block:   ", paste(result$var_names[result$blocks$cyclic_idx], collapse = ", "), "\n")
cat("  Uncorrelated block:         ", paste(result$var_names[result$blocks$uncorr_idx], collapse = ", "), "\n")
cat("Time points needing PSD correction:", result$n_corrections, "out of", T, "\n")


# graphics
result$dataset %>% pivot_longer(V1:V9) %>%
  filter(subject==1) %>%
  ggplot(aes(x=time, y = value, col=name, group=name)) +
  geom_line()+
  labs(x="Time", y = " ", col = "variable")

## constant cor
result$dataset %>%
  select(subject, time, V1, V2, V3) %>%
  group_by(time) %>%
  summarise(cor12 = cor(V1, V2, method = "spearman"),
            cor23 = cor(V3, V2, method = "spearman"),
            cor13 = cor(V1, V3, method = "spearman")) %>%
  pivot_longer(cor12:cor13) %>%
  ggplot(aes(x=time, y = value, col=name, group=name)) +
  geom_line()+
  ylim(0, 1)+
  labs(x="Time", y = "Cor", col = "Pair")

result$corr_by_time[[1]]

## cyclic cor
result$dataset %>%
  select(subject, time, V4, V5, V6) %>%
  group_by(time) %>%
  summarise(cor45 = cor(V4, V5, method = "spearman"),
            cor46 = cor(V4, V6, method = "spearman"),
            cor56 = cor(V5, V6, method = "spearman")) %>%
  pivot_longer(cor45:cor56) %>%
  ggplot(aes(x=time, y = value, col=name, group=name)) +
  geom_line()+
  ylim(-1, 1)+
  labs(x="Time", y = "Cor", col = "Pair")


# Save to CSV
saveRDS(result$dataset, file ="Data_P9_T100_N100.rds")
