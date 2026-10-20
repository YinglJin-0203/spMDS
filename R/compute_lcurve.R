#' Construct an L-curve by fitting spMDS across a range of lambda values
#'
#' @param delta_list list of T dissimilarity matrices
#' @param phi_list list of T spline basis vectors
#' @param P, K, D dimensions (nodes, basis functions, embedding dim)
#' @param lambda_seq numeric vector of lambda values to try (any order;
#'   internally sorted descending so warm-starts go from most- to
#'   least-regularized, which tends to be more numerically stable)
#' @param rho ADMM penalty parameter (held fixed across lambda)
#' @param admm_iter, max_iter, tol passed through to spMDS at each lambda
#' @param warm_start if TRUE, use previous lambda's fitted B as the next
#'   B_init, which substantially speeds up the sweep
#' @param verbose_fit passed to spMDS's own verbose flag (per-iteration
#'   printing) -- usually FALSE here since you already get progress
#'   printed per lambda
#' @return a data.frame with columns lambda, stress, penalty, log_stress,
#'   log_penalty, plus a list column-like attribute "fits" holding each
#'   full spMDS result object (attr(result, "fits"))
compute_lcurve <- function(delta_list, phi_list, P, K, D,
                           lambda_seq, rho = 1,
                           admm_iter = 15, max_iter = 1000, tol = 1e-5,
                           warm_start = TRUE, verbose_fit = FALSE) {

  lambda_seq <- sort(lambda_seq, decreasing = TRUE)
  n_lambda <- length(lambda_seq)

  stress_vec <- numeric(n_lambda)
  penalty_vec <- numeric(n_lambda)
  movement_vec <- numeric(n_lambda)
  time_vec <- numeric(n_lambda)

  fits <- vector("list", n_lambda)

  # Precompute Phi_list once -- shared across all lambda fits
  Phi_list <- lapply(phi_list, function(phi_t) build_Phi_t(phi_t, P))

  B_warm <- NULL

  for (i in seq_len(n_lambda)) {
    lam <- lambda_seq[i]
    cat(sprintf("\n--- Fitting lambda = %.5g (%d/%d) ---\n", lam, i, n_lambda))

    t_elapsed <- system.time({
      fit_i <- spMDS(
        delta_list = delta_list,
        phi_list = phi_list,
        P = P, K = K, D = D,
        lambda = lam, rho = rho,
        admm_iter = admm_iter, max_iter = max_iter,
        tol = tol, B_init = if (warm_start) B_warm else NULL,
        verbose = verbose_fit
      )
    })

    sp <- compute_stress_and_penalty(delta_list, fit_i$C_list, fit_i$B)
    stress_vec[i] <- sp$stress
    penalty_vec[i] <- sp$penalty
    movement_vec[i] <- sp$movement
    time_vec[i] <- t_elapsed[3]
    fits[[i]] <- fit_i

    if (warm_start) B_warm <- fit_i$B

    cat(sprintf("  stress = %.4f, penalty (||B||_1) = %.4f, time = %.1fs\n",
                sp$stress, sp$penalty, t_elapsed["elapsed"]))
  }

  lcurve_df <- data.frame(
    lambda = lambda_seq,
    stress = stress_vec,
    penalty = penalty_vec,
    movement = movement_vec,
    time = time_vec
  )

  attr(lcurve_df, "fits") <- fits
  return(lcurve_df)
}
