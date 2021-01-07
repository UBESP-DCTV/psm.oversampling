# Analysis: Monte carlo simulations ------------------------------------
library(tidyverse)
library(Matching)
library(cobalt)
library(furrr)
library(tictoc)
options(mc.cores = future::availableCores() - 1)

set.seed(456987)

# 1) Load the data -----------------------------------------------------
n <- 1000000
n_sim <- 1000

load(here::here(
  glue::glue("simulation_data_pop_1000000_n_{n_sim}.rda")
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

simulation_results <- furrr::future_map2(
  .x = sim_data$samples, .y = sim_data$att_true,
  ~ multi_scenario_analysis(
    mc_samples = .x,
    att_true = .y,
    covariates = covariates,
    treatment = treatment,
    outcome = outcome,
    rep_over_grid = rep_over_grid
  ), .progress = TRUE
)

toc()

sim_res <- sim_data %>%
  dplyr::mutate(
    results = simulation_results
  )

# Save everything ------------------------------------------------------
save(
  simulation_results, sim_res,
  file = here::here(
    glue::glue(
      "simulation_results_pop_1000000_n_{n_sim}.rda"
    )
  )
)


