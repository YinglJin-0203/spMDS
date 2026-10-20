
# --------------- Matrix calculations --------------------------

## -------------------------------------------------------------
## 1. Utility: build Phi_t = t(phi_t) %x% I_P   (P x PK)
## -------------------------------------------------------------
build_Phi_t <- function(phi_t, P) {
  # phi_t: K x 1 vector of spline basis values at time t
  # returns: P x (P*K) matrix
  K <- length(phi_t)
  I_P <- Diagonal(P)
  Phi_t <- kronecker(I_P, t(matrix(phi_t, nrow = K, ncol = 1)))
  Phi_t <- as.matrix(Phi_t)
  return(Phi_t)  # P x (P*K)
}


# --------------- Guttman transformation -----------------------
## -------------------------------------------------------------
## 4. Compute Guttman transform weight matrix W_t
##
##    W_ij,t = - delta_ij,t / d_ij,t   (i != j)
##    W_ii,t = - sum_{j != i} W_ij,t
##
##    Handles d_ij,t == 0 (or very small) by setting the
##    corresponding weight contribution to zero, which is the
##    standard convention in SMACOF-type majorization to avoid
##    division by zero.
## -------------------------------------------------------------

compute_Wt <- function(delta_t, d_t, tol = 1e-10) {
  P <- nrow(delta_t)
  W <- matrix(0, P, P)
  for (i in 1:P) {
    for (j in 1:P) {
      if (i != j) {
        if (d_t[i, j] > tol) {
          W[i, j] <- -delta_t[i, j] / d_t[i, j]
        } else {
          W[i, j] <- 0
        }
      }
    }
  }
  for (i in 1:P) {
    W[i, i] <- -sum(W[i, -i])
  }
  return(W)
}

# --------------- B-update --------------------------
## -------------------------------------------------------------
## 5. B-update (ridge-type closed form from ADMM subproblem)
##
##    B^(m+1) = (2*sum_t Phi_t' Phi_t + rho*I)^{-1} *
##              (2*sum_t Phi_t' W_t Phi_t B^(m-1) - Y + rho*Z)
##
##    Precompute the (PK x PK) matrix sum_t Phi_t'Phi_t and its
##    factorization once per outer MM iteration (Phi_t doesn't
##    change across m, so this can even be cached across the
##    whole algorithm).
## -------------------------------------------------------------
update_B <- function(Phi_list, Wt_list, B_prev, Y, Z, rho,
                     PtP_sum = NULL, PK = NULL, P = NULL) {
  T_len <- length(Phi_list)
  if (is.null(PK)) PK <- ncol(Phi_list[[1]])
  if(is.null(P)) P <- nrow(Phi_list[[1]])

  if (is.null(PtP_sum)) {
    PtP_sum <- Reduce(`+`, lapply(Phi_list, function(Phi_t) t(Phi_t) %*% Phi_t))
  }

  # RHS term: 2 * sum_t Phi_t' W_t Phi_t B_prev
  rhs_term <- Reduce(`+`, lapply(seq_len(T_len), function(t) {
    Phi_t <- Phi_list[[t]]
    Wt <- Wt_list[[t]]
    t(Phi_t) %*% (Wt %*% (Phi_t %*% B_prev))
  }))

  LHS <- 2 * P * PtP_sum + rho * diag(PK)
  RHS <- 2 * rhs_term - Y + rho * Z

  B_new <- solve(LHS, RHS)
  return(B_new)
}


# --------------- Z-update -------------------------------------


## -------------------------------------------------------------
## For temporal penalization: build the full (PK x PK) "sum of difference" matrix directly
##      from Phi_list (no separate phi_t vectors needed).
##
##    G = sum_{t=2}^T D_t' D_t,   D_t := Phi_t - Phi_{t-1}
##
##    This follows directly from linearity of the Kronecker
##    product: since Phi_t = phi_t' %x% I_P, we have
##    Phi_t - Phi_{t-1} = (phi_t' - phi_{t-1}') %x% I_P = D_t.
##    So D_t is obtained by simple matrix subtraction of the
##    Phi_t's you already build via build_Phi_t
##
##    Phi_list: list of T matrices (each P x PK), as produced by
##              build_Phi_t()
## -------------------------------------------------------------
build_Dsum <- function(Phi_list) {
  T_len <- length(Phi_list)
  PK <- ncol(Phi_list[[1]])
  G <- matrix(0, PK, PK)
  for (t in 2:T_len) {
    D_t <- Phi_list[[t]] - Phi_list[[t - 1]]
    G <- G + t(D_t) %*% D_t
  }
  return(G)  # PK x PK
}

## -------------------------------------------------------------
## For L1 penalization: soft_threshold proximal
## -------------------------------------------------------------

soft_threshold <- function(x, kappa) {
  sign(x) * pmax(abs(x) - kappa, 0)
}

# -------------------------------------------------------------
## 6. Z-update: now supports two penalty types
##
##    penalty_type = "l1" (default, original behavior):
##      Z^(k+1) = S_{lambda/rho}(B^(k+1) + Y/rho)      [soft-threshold]
##
##    penalty_type = "temporal":
##      Z^(k+1) = (2*lambda*G + rho*I)^{-1} (Y + rho*B)  [ridge solve]
##      requires G = build_Gsum(Phi_list), precomputed once
##      outside the ADMM loop.
##
##    Arguments:
##      B, Y, rho, lambda : as before
##      penalty_type      : "l1" or "temporal"
##      G                 : (PK x PK) matrix from build_Gsum();
##                          REQUIRED if penalty_type = "temporal",
##                          ignored otherwise
## -------------------------------------------------------------

update_Z <- function(B, Y, rho, lambda, penalty_type = c("l1", "temporal"), Dsum = NULL) {
  penalty_type <- match.arg(penalty_type)

  if (penalty_type == "l1") {
    return(soft_threshold(B + Y / rho, lambda / rho))
  }

  # penalty_type == "temporal"
  if (is.null(Dsum)) {
    stop("penalty_type = 'temporal' requires build_Dsum(Phi_list)")
  }
  PK <- nrow(Dsum)
  LHS <- 2 * lambda * Dsum + rho * diag(PK)
  RHS <- Y + rho * B
  Z_new <- solve(LHS, RHS)
  return(Z_new)
}


## -------------------------------------------------------------
## 7. Y-update: dual ascent step
##
##    Note: the document writes Y^(m+1) = Y^m + |B - Z|, but this
##    is very likely a typo for the standard ADMM dual update
##    Y^(m+1) = Y^m + rho*(B - Z), which is what makes the
##    augmented Lagrangian's stationarity conditions consistent.
##    We implement the standard form here and flag the assumption
##    clearly, since using |B - Z| (absolute value) is not a valid
##    dual update in classical ADMM.
## -------------------------------------------------------------
update_Y <- function(Y, B, Z, rho) {
  Y + rho * (B - Z)
}

## -------------------------------------------------------------
## 8. Objective function evaluation: now supports two penalty types
##
##    penalty_type = "l1":
##      o(B) = stress + lambda * ||B||_1
##
##    penalty_type = "temporal":
##      o(B) = stress + lambda * sum_{t=2}^T ||C_t - C_{t-1}||_F^2
##      (equivalent to lambda * sum_t ||D_t B||_F^2, evaluated via
##      the fitted embeddings directly -- no need for G here)
## -------------------------------------------------------------
compute_objective <- function(delta_list, Phi_list, B, lambda,
                              penalty_type = c("l1", "temporal")) {
  penalty_type <- match.arg(penalty_type)

  T_len <- length(delta_list)
  # stress
  stress_total <- 0
  C_list <- vector("list", T_len)
  for (t in 1:T_len) {
    Ct <- Phi_list[[t]] %*% B
    C_list[[t]] <- Ct
    d_t <- as.matrix(dist(Ct))
    delta_t <- delta_list[[t]]
    idx_upper <- upper.tri(delta_t, diag = TRUE)
    stress_total <- stress_total + sum((delta_t[idx_upper] - d_t[idx_upper])^2)
  }

  # penalty
  if (penalty_type == "l1") {
    penalty <- sum(abs(B))
  } else {
    penalty <- 0
    if (T_len >= 2) {
      for (t in 2:T_len) {
        penalty <- penalty + sum((C_list[[t]] - C_list[[t - 1]])^2)
      }
    }
  }

  # movement
  movement <-  sum(sapply(2:T_len, function(t) sum((C_list[[t]] - C_list[[t-1]])^2)))

  return(c(stress = stress_total,
           penalty = penalty,
           movement = movement,
           objective = stress_total + lambda * penalty))
}

## -------------------------------------------------------------
## 2. Compute embeddings C_t = Phi_t %*% B for all t
## -------------------------------------------------------------
compute_Ct_list <- function(Phi_list, B) {
  lapply(Phi_list, function(Phi_t) Phi_t %*% B)  # each P x D
}


## -------------------------------------------------------------
## 4. Compute Guttman transform weight matrix W_t
##
##    W_ij,t = - delta_ij,t / d_ij,t   (i != j)
##    W_ii,t = - sum_{j != i} W_ij,t
##
##    Handles d_ij,t == 0 (or very small) by setting the
##    corresponding weight contribution to zero, which is the
##    standard convention in SMACOF-type majorization to avoid
##    division by zero.
## -------------------------------------------------------------
# compute_Wt <- function(delta_t, d_t, tol = 1e-10) {
#   P <- nrow(delta_t)
#   W <- matrix(0, P, P)
#   off_diag <- (row(W) != col(W))
#   ratio <- matrix(0, P, P)
#   valid <- off_diag & (d_t > tol)
#   ratio[valid] <- -delta_t[valid] / d_t[valid]
#   W[off_diag] <- ratio[off_diag]
#   diag(W) <- -rowSums(W) + diag(W)  # since diag currently 0, this just negates row sums
#   diag(W) <- -rowSums(W - diag(diag(W)))  # ensure diagonal = -sum of off-diag row entries
#   return(W)
# }

# =============================================================
## L-curve construction over a range of lambda values
## =============================================================
## For each lambda, fit spMDS and record:
##   - stress: sum_t sum_{i<=j} (delta_ij_t - d_ij_t)^2   (data fidelity)
##   - penalty: ||B||_1                                    (regularization)
## The L-curve plots log(stress) vs log(penalty); the "corner" of the
## resulting L-shape is a common heuristic for choosing lambda.
## =============================================================

#' Split the objective into its stress and penalty components
#'
#' @param delta_list list of T dissimilarity matrices
#' @param Ct_list list of T precomputed Phi_t matrices (P x PK)
#' @param B fitted coefficient matrix (PK x D)
#' @return list with stress and penalty (unweighted, i.e. penalty = ||B||_1,
#'         not lambda*||B||_1 -- this is what the L-curve needs)
# compute_stress_and_penalty <- function(delta_list, C_list, B) {
#   T_len <- length(delta_list)
#   stress_total <- 0
#   for (t in 1:T_len) {
#     Ct <- C_list[[t]]
#     d_t <- as.matrix(dist(Ct))
#     delta_t <- delta_list[[t]]
#     idx_upper <- upper.tri(delta_t, diag = TRUE)
#     stress_total <- stress_total + sum((delta_t[idx_upper] - d_t[idx_upper])^2)
#   }
#   penalty <- sum(abs(B))
#   movement <-  sum(sapply(2:nT, function(t) sum((C_list[[t]] - C_list[[t-1]])^2)))
#   list(stress = stress_total, penalty = penalty, movement = movement)
# }
