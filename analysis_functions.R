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
  ps <- ps_data$ps_logit

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

# Evaluate balance -----------------------------------------------------
balance_diagnostics <- function(
  matching_obj, ps_data, treatment, covariates
) {

  assertive::assert_is_data.frame(ps_data)
  if(class(matching_obj) != "Match") {
    usethis::ui_stop(
      "'matching_obj' must be of class 'Match'"
    )
  }
  assertive::assert_is_character(treatment)
  assertive::assert_is_character(covariates)

  # Treatment formula
  treat_formula <- paste(
    treatment,
    paste0(c(covariates, "ps_logit"), collapse = " + "),
    sep = " ~ "
  ) %>%
    as.formula()

  # Create an if else condition to discard Matching that were not
  # performed given the absence of valid control matches

  if (!is.na(matching_obj)[1]) {

  # One-dimensional measures
  bal_tab <- cobalt::bal.tab(
    x = matching_obj,
    formula = treat_formula,
    data = ps_data,
    continuous = "std", binary = "std",
    s.d.denom = "pooled", abs = TRUE, quick = FALSE
  )

  smds <- bal_tab$Balance %>%
    tibble::as_tibble(rownames = "covs") %>%
    dplyr::slice(-nrow(.)) %>%
    dplyr::select(covs, Diff.Adj)

  # Overall balance measure
  asmds <- mean(smds$Diff.Adj, na.rm = TRUE)
  ps_smd <- bal_tab$Balance %>%
    tibble::as_tibble(rownames = "covs") %>%
    dplyr::filter(covs == "ps_logit") %>%
    .[["Diff.Adj"]]
  ps_ovc <- bal_tab$Balance %>%
    tibble::as_tibble(rownames = "covs") %>%
    dplyr::filter(covs == "ps_logit") %>%
    .[["OVL.Adj"]]
  prop_matched_treated <- bal_tab$Observations[2, 2]/
    bal_tab$Observations[1, 2]

  # Final tibble with all the measures
  res <- tibble::tibble(
    single_smds = list(smds),
    average_smd = asmds,
    ps_smd = ps_smd,
    ps_ovc = ps_ovc,
    p_matched_treated = prop_matched_treated
  )

  } else {

    res <- tibble::tibble(
      single_smds = list(tibble::tibble(
        covs = covariates,
        Diff.Adj = rep(NA_real_, length(covariates))
      )),
      average_smd = NA_real_,
      ps_smd = NA_real_,
      ps_ovc = NA_real_,
      p_matched_treated = NA_real_
    )

  }

  res

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
single_analysis <- function(
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

  # Step 3: balance evaluation
  balance <- balance_diagnostics(
    matching_obj = mtch, ps_data = ps_est, treatment = treatment,
    covariates = covariates
  )

  # Step 4: estimate treatment effect
  tr_eff <- est_treat_effect(matching_obj = mtch)

  # Store results into a tibble
  tibble::tibble(balance, tr_eff)

}

# Run the analysis on all the datasets ---------------------------------
mc_analysis <- function(
  mc_samples, covariates, treatment, outcome, ratio, replacement
) {

  assertive::assert_is_list(mc_samples)
  assertive::assert_is_character(covariates)
  assertive::assert_is_character(treatment)
  assertive::assert_is_character(outcome)
  assertive::assert_is_a_number(ratio)
  assertive::assert_is_logical(replacement)

  purrr::map_dfr(
    .x = mc_samples,
    ~ {

      single_analysis(
        data = .x,
        covariates = covariates,
        treatment = treatment,
        outcome = outcome,
        ratio = ratio,
        replacement = replacement
      )

    }
  )

}

# Computes performances metrics ----------------------------------------
mc_metrics <- function(mc_results, att_true) {

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

  # Balance metrics ----------------------------------------------------
  # Average single smds
  single_smds <- purrr::map_dfr(
    .x = mc_results$single_smds,
    ~ .x
  ) %>%
    dplyr::group_by(covs) %>%
    dplyr::summarise(asmd = mean(Diff.Adj, na.rm = TRUE)) %>%
    dplyr::ungroup()

  # Average ASMD, PS-SMD, OVC, AUC and Proportion of matched treated
  asmd <- mean(mc_results$average_smd, na.rm = TRUE)
  ps_smd <- mean(mc_results$ps_smd, na.rm = TRUE)
  ps_ovc <- mean(mc_results$ps_ovc, na.rm = TRUE)
  p_mt <- mean(mc_results$p_matched_treated, na.rm = TRUE)

  # Store results into a tibble
  tibble::tibble(
    rel_bias = rel_bias,
    rmse = rmse,
    nc_ai = nc_ai, nc_std = nc_std,
    single_smds = list(single_smds),
    asmd = asmd, ps_smd = ps_smd, ps_ovc = ps_ovc, p_mt = p_mt
  )


}

# Single analysis for a scenario ---------------------------------------
scenario_analysis <- function(
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
  mc_results <- mc_analysis(
    mc_samples = mc_samples,
    covariates = covariates,
    treatment = treatment,
    outcome = outcome,
    ratio = ratio,
    replacement = replacement
  )

  # Get MC performances metrics
  mc_metrics(
    mc_results = mc_results, att_true = att_true
  ) %>%
    tibble::add_column(ratio = ratio, .before = 1L) %>%
    tibble::add_column(replace = replacement, .before = 1L)
}

# Multiple analysis for scenario ---------------------------------------
multi_scenario_analysis <- function(
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
    ~ scenario_analysis(
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

# Estimate treatment effect for case study -----------------------------
est_treat_effect_cs <- function(matching_obj) {

  if(class(matching_obj) != "Match") {
    usethis::ui_stop(
      "'matching_obj' must be of class 'Match'"
    )
  }

  # Create an if else condition to discard Matching that were not
  # performed given the absence of valid control matches

  if (!is.na(matching_obj)[1]) {

    rd <- matching_obj$est[1, 1]
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

# Run a single analysis for the case study -----------------------------
single_analysis_cs <- function(
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

  prop_contr_rep <- sum(table(mtch$index.control) > 1)/
    length(table(mtch$index.control))

  # Step 4: balance evaluation
  balance <- balance_diagnostics(
    matching_obj = mtch, ps_data = ps_est, treatment = treatment,
    covariates = covariates
  )

  # Step 5: estimate treatment effect
  tr_eff <- est_treat_effect_cs(matching_obj = mtch)

  # Store results into a tibble
  tibble::tibble(balance, p_contr_rep = prop_contr_rep, tr_eff)

}

