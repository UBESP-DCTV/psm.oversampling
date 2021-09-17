# Functions for the analysis -------------------------------------------
# Estimate PS ----------------------------------------------------------
ps_estimate <- function(data, covariates, treatment) {

  assertive::assert_is_data.frame(data)
  assertive::assert_is_character(covariates)
  assertive::assert_is_character(treatment)

  # Construct the formula for PS
  ps_form <- paste(
    treatment,
    paste0(covariates, collapse = " + "),
    sep = " ~ "
  ) %>%
    stats::as.formula()

  # PS model
  ps_fit <- stats::glm(
    ps_form, data = data, family = stats::binomial("logit")
  )

  # Add PS (on the logit scale) to the data
  data %>%
    dplyr::mutate(ps_logit = stats::predict(ps_fit))

}

# Matching -------------------------------------------------------------
matching_fun <- function(
  ps_data, outcome, treatment, covariates, ratio, replacement
) {

  assertive::assert_is_data.frame(ps_data)
  assertive::assert_is_character(covariates)
  assertive::assert_is_character(treatment)
  assertive::assert_is_character(outcome)
  assertive::assert_is_a_number(ratio)
  assertive::assert_is_logical(replacement)

  # Define the inputs for Matching
  y <- ps_data[[outcome]]
  tr <- ps_data[[treatment]]
  ps <- ps_data[["ps_logit"]]

  # # Define the caliper
  # caliper <- 0.2 * sd(ps)

  # Perform the matching
  Matching::Match(
    Y = y,
    Tr = tr,
    X = ps,
    estimand ="ATT",
    M = ratio,
    caliper = 0.2,
    replace = replacement,
    ties = TRUE,
    Var.calc = 1
  )

}

# Estimate treatment effect --------------------------------------------
est_treat_effect <- function(matching_obj) {

  if(class(matching_obj) != "Match") {
    usethis::ui_stop(
      "'matching_obj' must be of class 'Match'"
    )
  }

  # Create an if else condition to discard Matching that were not
  # performed given the absence of valid control matches

  if (!is.na(matching_obj)[1]) {

  rd <- abs(matching_obj$est[1, 1])
  se_ai <- matching_obj$se
  se_std <- matching_obj$se.standard
  rd_lower_ai <- rd - 1.96 * se_ai
  rd_upper_ai <- rd + 1.96 * se_ai
  rd_lower_std <- rd - 1.96 * se_std
  rd_upper_std <- rd + 1.96 * se_std

  # Store results into a tibble
  res <- tibble::tibble(
    rd_est = rd,
    rd_lower_ai = ifelse(purrr::is_empty(rd_lower_ai), NA, rd_lower_ai),
    rd_upper_ai = ifelse(purrr::is_empty(rd_upper_ai), NA, rd_upper_ai),
    rd_lower_std = rd_lower_std,
    rd_upper_std = rd_upper_std
  )

  } else {

    res <- tibble::tibble(
      rd_est = NA_real_,
      rd_lower_ai = NA_real_,
      rd_upper_ai = NA_real_,
      rd_lower_std = NA_real_,
      rd_upper_std = NA_real_
    )

  }

  res

}

# Run a single analysis ------------------------------------------------
single_analysis_linear <- function(
  data, covariates, treatment, outcome, ratio, replacement
) {

  assertive::assert_is_data.frame(data)
  assertive::assert_is_character(covariates)
  assertive::assert_is_character(treatment)
  assertive::assert_is_character(outcome)
  assertive::assert_is_a_number(ratio)
  assertive::assert_is_logical(replacement)

  # Step 1: estimate PS
  ps_est <- ps_estimate(data, covariates, treatment)

  # Step 2: perform matching
  mtch <- matching_fun(
    ps_data = ps_est, covariates = covariates,
    treatment = treatment, outcome = outcome,
    ratio = ratio, replacement = replacement
  )

  # Step 3: estimate treatment effect
  tr_eff <- est_treat_effect(matching_obj = mtch)

  # Store results into a tibble
  tibble::tibble(tr_eff)

}

# Run the analysis on all the datasets ---------------------------------
mc_metrics_linear <- function(mc_results, att_true) {

  assertive::assert_is_data.frame(mc_results)
  assertive::assert_is_a_number(att_true)

  # bias and relative bias
  bias <- mean(mc_results$rd_est, na.rm = TRUE) - att_true
  rel_bias <- abs(100 * (bias/att_true))

  # RMSE
  rmse <- sqrt(
    sum(mc_results$rd_est - att_true, na.rm = TRUE)^2/nrow(mc_results)
  )

  # 95% nominal coverage
  incl_att_ai <- att_true >= mc_results$rd_lower_ai &
    att_true <= mc_results$rd_upper_ai

  incl_att_std <- att_true >= mc_results$rd_lower_std &
    att_true <= mc_results$rd_upper_std

  nc_ai <- sum(incl_att_ai, na.rm = TRUE)/length(incl_att_ai)
  nc_std <- sum(incl_att_std, na.rm = TRUE)/length(incl_att_std)

  # Store results into a tibble
  tibble::tibble(
    rel_bias = rel_bias,
    rmse = rmse,
    nc_ai = nc_ai, nc_std = nc_std
  )


}

# Single analysis for a scenario ---------------------------------------
scenario_analysis_linear <- function(
  mc_samples, att_true, covariates, treatment, outcome, ratio,
  replacement
) {

  assertive::assert_is_list(mc_samples)
  assertive::assert_is_a_number(att_true)
  assertive::assert_is_character(covariates)
  assertive::assert_is_character(treatment)
  assertive::assert_is_character(outcome)
  assertive::assert_is_a_number(ratio)
  assertive::assert_is_logical(replacement)

  # Get MC results
  mc_results <- mc_analysis_linear(
    mc_samples = mc_samples,
    covariates = covariates,
    treatment = treatment,
    outcome = outcome,
    ratio = ratio,
    replacement = replacement
  )

  # Get MC performances metrics
  mc_metrics_linear(
    mc_results = mc_results, att_true = att_true
  ) %>%
    tibble::add_column(ratio = ratio, .before = 1L) %>%
    tibble::add_column(replace = replacement, .before = 1L)
}

# Multiple analysis for scenario ---------------------------------------
multi_scenario_analysis_linear <- function(
  mc_samples, att_true, covariates, treatment, outcome, rep_over_grid
) {

  assertive::assert_is_list(mc_samples)
  assertive::assert_is_a_number(att_true)
  assertive::assert_is_character(covariates)
  assertive::assert_is_character(treatment)
  assertive::assert_is_character(outcome)
  assertive::assert_is_data.frame(rep_over_grid)

  purrr::map2_dfr(
    .x = rep_over_grid[["replacement"]], .y = rep_over_grid[["ratio"]],
    ~ scenario_analysis_linear(
      mc_samples = mc_samples,
      att_true = att_true,
      covariates = covariates,
      treatment = treatment,
      outcome = outcome,
      ratio = .y,
      replacement = .x
    )
  )


}
