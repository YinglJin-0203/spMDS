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
source("R/compute_lcurve.R")
source("R/plot_lcurve.R")

set.seed(123)








#### Example ####
df <- read_rds("Manuscripts/Simulation/Data/Data_P9_T100_N100.rds")
head(df)

#### fit model ####

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

# lambda search
lambda_seq <- 10^seq(-3, 1, length.out = 10)  # log-spaced grid
lcurve_df <- compute_lcurve(
  delta_list = delta_list, phi_list = phi_list,
  P = P, K = K, D = 2,
  lambda_seq = lambda_seq, rho = 1,
  admm_iter = 15, max_iter = 1000, tol = 1e-5,
  warm_start = TRUE, verbose_fit = FALSE
)

# result

plot_lcurve(lcurve_df)

lcurve_df %>%
  mutate(obj = stress+lambda*penalty) %>%
  ggplot()+
  geom_line(aes(x=lambda, y = stress, col = "Stress"))+
  geom_line(aes(x=lambda, y = penalty, col = "regularization"))+
  geom_point(aes(x=lambda, y = stress, col = "Stress"))+
  geom_point(aes(x=lambda, y = penalty, col = "regularization"))+
  scale_x_log10()

sum(lcurve_df$time)

