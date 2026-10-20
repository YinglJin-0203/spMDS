## =============================================================
## Example usage with simulated data
## =============================================================
library(splines)
library(Matrix)
library(igraph)
library(tidyverse)
library(magick)
source("R/Helpers.R")
source("R/spMDS.R")

set.seed(123)

# data
df <- read_rds("Manuscripts/Simulation/Data/Data_P9_T100_N100.rds")
head(df)

# dissimilarity
delta_list <- df %>% select(-subject) %>% group_by(time) %>%
  group_map(~cor(.x), method = "spearman")
delta_list <- lapply(delta_list, function(x){1-abs(x)})

# splines basis
nT <- length(delta_list)
K <- 30
P <- dim(delta_list[[1]])[1]
knots <- seq(0, nT, length.out = K - 2)[-c(1, K - 2)]
basis_obj <- bs(1:nT, degree = 3, knots = knots, intercept = TRUE)
phi_list <- lapply(1:nT, function(t) as.numeric(basis_obj[t, ]))


#### Test functions ####
Phit <- build_Phi_t(phi_list[[2]], P)
dim(Phit)
View(Phit)

Phi_list <- lapply(phi_list, function(phi_t) build_Phi_t(phi_t, P))
# Precompute sum_t Phi_t' Phi_t once (does not depend on B or m)
PtP_sum <- Reduce(`+`, lapply(Phi_list, function(Phi_t) t(Phi_t) %*% Phi_t))

# Initialize B, Z, Y
PK <-P*K
D <- 2
B <- matrix(rnorm(PK * D, sd = 0.1), PK, D)
Z <- B
Y <- matrix(0, PK, D)

lambda = 0.1

compute_objective(delta_list, Phi_list, B, lambda)

Ct_list <- compute_Ct_list(Phi_list, B)
T_len <- length(Ct_list)
Wt_list <- vector("list", T_len)
for (t in 1:T_len) {
  d_t <-   as.matrix(dist(Ct_list[[t]]))
  Wt_list[[t]] <- compute_Wt(delta_list[[t]], d_t)
}


rho <- 0.1
for (inner in 1:admm_iter) {
  B <- update_B(Phi_list, Wt_list, B, Y, Z, rho, PtP_sum, PK)
  Z <- update_Z(B, Y, rho, lambda, Dsum)
  Y <- update_Y(Y, B, Z, rho)
}

spMDS_fit <- spMDS(delta_list, phi_list, P, K, 2, 0.1, 1, max_iter = 500)

#### Fit model ####

# dissimilarity
delta_list <- df %>% select(-subject) %>% group_by(time) %>%
  group_map(~cor(.x), method = "spearman")
delta_list <- lapply(delta_list, function(x){1-abs(x)})

# splines basis
nT <- length(delta_list)
K <- 30
P <- dim(delta_list[[1]])[1]
knots <- seq(0, nT, length.out = K - 2)[-c(1, K - 2)]
basis_obj <- bs(1:nT, degree = 3, knots = knots, intercept = TRUE)
phi_list <- lapply(1:nT, function(t) as.numeric(basis_obj[t, ]))


# hyper parameters
lambda_vec <- c(0.01, 0.1, 1)
rho_vec <- c(0.01, 0.1, 1)

df_combo <- expand.grid(lambda_vec = lambda_vec, rho_vec = rho_vec)
df_combo[ , c('stress', "movement", "time")] <- NA

# loop over rows
for (i in 1:nrow(df_combo)) {
  lambda <- df_combo$lambda_vec[i]
  rho <- df_combo$rho_vec[i]
  time_i <- system.time({
    spMDS_fit <- spMDS(
      delta_list = delta_list,
      phi_list = phi_list,
      penalty_type = "l1",
      P = P, K = K, D = 2,
      lambda = lambda, rho = rho,
      admm_iter = 15, max_iter = 1000,
      tol = 1e-5, verbose = TRUE
    )
  })

  df_combo[i ,"time"] <- time_i[3]
  df_combo[i , "stress"] <- tail(spMDS_fit$stress, 1)
  df_combo[i , "movement"] <- tail(spMDS_fit$movement, 1)
}


#### Overall inspection ####
spMDS_fit$C_list

## Plot objective trajectory
plot(spMDS_fit$objective_history, type = "b", pch = 20,
     xlab = "MM Iteration", ylab = "Objective value",
     main = "Convergence of MM-ADMM Splines MDS")

plot(spMDS_fit$stress[1:286])
plot(spMDS_fit$movement[1:286])
# network plot
inet <- make_empty_graph(n=P, directed = F)
net_list <- list()

for (t in seq_len(nT)) {
  layout_mat <- cbind(spMDS_fit$C_list[[t]])
  png(filename = "Manuscripts/Simulation/Figures/snapshot.png")
  plot(inet, layout = layout_mat,
       vertex.size = 30,
       vertex.label.cex = 2,
       vertex.label.color = "black",
       main = paste0("t =", t, ", lambda = ", 0.5, ", rho = ", 1))
  dev.off()

  net_list[[t]] <- image_read("Manuscripts/Simulation/Figures/snapshot.png")
}


## animate at 2 frames per second
net_gif <- image_animate(image_join(net_list), fps = 2)

## view animated image
net_gif

## save to disk
image_write(image = net_gif,
            path = "Manuscripts/Simulation/Figures/GIF_N100_P9_T100_lambda0.5_rho0.1.gif")


# stress and regularization
## stress
stress_total <- 0
for (t in 1:nT) {
  Ct <- spMDS_fit$C_list[[t]]
  d_t <- as.matrix(dist(Ct))
  delta_t <- delta_list[[t]]
  idx_upper <- upper.tri(delta_t, diag = TRUE)
  stress_total <- stress_total + sum((delta_t[idx_upper] - d_t[idx_upper])^2)
}
stress_total

# penalization
sum(abs(spMDS_fit$B))

# movement
movement <- sum(sapply(2:nT, function(t) sum((spMDS_fit$C_list[[t]] - spMDS_fit$C_list[[t-1]])^2)))
movement

#### Detail Inspection ####
# compare actual corrlation, observed correlation, distance

vars <- setdiff(colnames(df), c("time", "subject"))

# get all unique pairs
pairs <- combn(vars, 2, simplify = FALSE)

# compute correlation for each pair at each time point
cor_df <- map_dfr(pairs, function(pair) {
  df %>%
    group_by(time) %>%
    summarise(
      var1        = pair[1],
      var2        = pair[2],
      correlation = cor(.data[[pair[1]]], .data[[pair[2]]],
                        use = "pairwise.complete.obs"),
      .groups = "drop"
    )
}) %>%
  select(time, var1, var2, correlation) %>%  # reorder columns
  arrange(time, var1, var2)                  # sort by time then variable pair

head(cor_df)

# compute euclidean distance for each pair at each time point
dist_df <- map_dfr(seq_along(spMDS_fit$C_list), function(t_idx) {
  coords <- spMDS_fit$C_list[[t_idx]]
  rownames(coords) <- paste0("V", 1:P)

  map_dfr(pairs, function(pair) {
    coord1 <- coords[pair[1], ]
    coord2 <- coords[pair[2], ]
    dist   <- sqrt(sum((coord1 - coord2)^2))

    data.frame(
      time     = t_idx-1,   # your vector of unique time points
      var1     = pair[1],
      var2     = pair[2],
      distance = dist
    )
  })
})

# join with the correlation data frame
cor_df <- cor_df %>%
  left_join(dist_df, by = c("time", "var1", "var2"))

head(cor_df)


cor_df %>%
  ggplot()+
  # geom_line(aes(x=time, y=correlation))+
  geom_line(aes(x=time, y=distance))+
  facet_grid(rows = vars(var1), cols = vars(var2))

