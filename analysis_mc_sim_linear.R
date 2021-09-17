# Analysis: Monte carlo simulations ------------------------------------
library(tidyverse)
library(Matching)
library(cobalt)
library(furrr)
library(tictoc)
library(progressr)
options(mc.cores = future::availableCores() - 1)

set.seed(456987)

# 1) Load the data -----------------------------------------------------
n <- 1000000
n_sim <- 10000

load(here::here(
  glue::glue("simulation_data_pop_1000000_n_{n_sim}_linear.rda")
))

source(here::here("analysis_functions.R"))

# 2) Define inputs for the analysis ------------------------------------
covariates <- c("x1", "x2", "x3", "x4", "x5", "x6")
treatment <- "treatment"
outcome <- "y"
method <- "nearest"
ratio <- 1:5
order <- "largest"
replacement <- c(FALSE, TRUE)
rep_over_grid <- tidyr::expand_grid(
  replacement = replacement,
  ratio = ratio
)

plan(multisession)

tic()

simulations_results <- furrr::future_map(
  .x = sim_data$samples,
  ~ multi_scenario_analysis_linear(
      mc_samples = .x,
      att_true = 0.15,
      covariates = covariates,
      outcome = outcome,
      treatment = treatment,
      rep_over_grid = rep_over_grid
    ),
  .options = furrr::furrr_options(seed = TRUE), .progress = TRUE
)

toc()

sim_res <- sim_data %>%
  dplyr::mutate(
    results = simulations_results
  ) %>%
  dplyr::select(-samples)

# Save everything ------------------------------------------------------

save(
  sim_res, tab_parms,
  file = here::here(
    glue::glue(
      "simulation_results_pop_1000000_n_{n_sim}_linear.rda"
    )
  )
)


