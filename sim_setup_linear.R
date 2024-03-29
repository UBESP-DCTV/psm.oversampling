# Oversampling: simulation setup with linear prob model ----------------
library(tidyverse)
library(doParallel)
library(MASS)
library(furrr)
library(tictoc)
options(mc.cores = future::availableCores() - 1)

set.seed(456987)

source(here::here("sim_functions_linear.R"))

# 1) Super populations -------------------------------------------------
# 1A) Define inputs ----------------------------------------------------
n <- 1000000L
mus <- rep(0, 6L)
sds <- rep(1, 6L)
corr_matrix <- matrix(
  c(
    1, 0.3, 0.4, 0.1, 0.2, 0.5,
    0.3, 1, 0.2, 0.1, 0.3, 0.2,
    0.4, 0.2, 1, 0.2, 0.2, 0.2,
    0.1, 0.1, 0.2, 1, 0.5, 0.4,
    0.2, 0.3, 0.2, 0.5, 1, 0.2,
    0.5, 0.2, 0.2, 0.4, 0.2, 1
  ),
  nrow = 6L, ncol = 6L, byrow = TRUE
)

beta_list <- list(
  "weak" = c(
    log(1.25), log(1.5), log(1.25), log(1.5), log(1.25), log(1.5)
  ),
  "strong" = c(
    log(1.5), log(1.75), log(1.5), log(1.75), log(1.5), log(1.75)
  )
)

intercept_treat <- seq(from = -3, to = 3, by = 0.05)
p_treat <- c(0.3, 0.5, 0.7)
delta_treat <- 0.01
alpha <- c(0.25, -0.01, 0.02, 0.003, -0.005, -0.001)
intercept_outcome_grid <- seq(from = 0.1, to = 0.4, by = 0.01)
marginal_rd <- -0.15
p_outcome <- 0.2

# 1B) Simulate super-population ----------------------------------------
scenario_grid <- tidyr::expand_grid(
  beta = beta_list, p_treat = p_treat
)

plan(multisession)

tic()

scenarios <- furrr::future_map2_dfr(
  .x = scenario_grid[["beta"]], .y = scenario_grid[["p_treat"]],
  ~ sim_pop(
    n = n, mus = mus, sds = sds, corr_matrix = corr_matrix,
    beta = .x, intercept_treat = intercept_treat, p_treat = .y,
    delta_treat = delta_treat, alpha = alpha,
    outcome_grid = intercept_outcome_grid, p_outcome = p_outcome,
    marginal_rd = marginal_rd
  ), .progress = TRUE, .options = furrr::furrr_options(seed = TRUE)
) %>%
  mutate(
    treat_assign = rep(c("weak", "strong"), each = nrow(.)/2),
    p_treat = rep(p_treat, 2L)
  )

toc()

# 2) Sample from the superpopulation -----------------------------------
n_samples <- 10000L
sample_size <- c(100L, 250L, 500L, 1000L)

scen_data <- tidyr::expand_grid(
  scenarios,
  ss = sample_size
)

plan(multisession)

tic()

sim_samples <- furrr::future_map2(
  .x = scen_data[["population"]], .y = scen_data[["ss"]],
  ~ sample_dfs(
    population = .x, sample_size = .y, n_samples = n_samples
  ), .progress = TRUE, .options = furrr::furrr_options(seed = TRUE)
)

toc()

sim_data <- scen_data %>%
  dplyr::mutate(samples = sim_samples) %>%
  dplyr::select(-population)

# 3) Save the data -----------------------------------------------------
save(
  scenarios, sim_data,
  file = here::here(
    glue::glue("simulation_data_pop_{n}_n_{n_samples}_linear.rda")
  )
)

