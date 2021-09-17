# Misc functions for simulations ---------------------------------------

# Generate baseline covariates -----------------------------------------
gen_baseline <- function(n, mus, sds, corr_matrix) {

  assertive::assert_is_a_number(n)
  assertive::assert_is_numeric(mus)
  assertive::assert_is_numeric(sds)
  assertive::assert_is_symmetric_matrix(corr_matrix)

  cov_matrix <- Matrix::nearPD(corr_matrix)$mat * (sds %*% t(sds))

  mv_norm <- MASS::mvrnorm(n = n, mu = mus, Sigma = cov_matrix)

  mv_norm[, 1] <- ifelse(mv_norm[, 1] > 0.7, 1, 0)
  mv_norm[, 2] <- ifelse(mv_norm[, 2] > 0.5, 1, 0)
  mv_norm[, 3] <- ifelse(mv_norm[, 3] > 0.9, 1, 0)

  colnames(mv_norm) <- glue::glue("x{1:ncol(mv_norm)}")

  mv_norm

}

# Proportion associated to an intercept --------------------------------
prop_treat <- function(x_matrix, beta, intercept) {

  assertive::assert_is_matrix(x_matrix)
  assertive::assert_is_numeric(beta)
  assertive::assert_is_a_number(intercept)

  coefs <- c(intercept, beta)
  des_matrix <- cbind(rep(1, nrow(x_matrix)), x_matrix)

  eta <- des_matrix %*% coefs

  treat <- stats::rbinom(
    n = nrow(x_matrix), size = 1, prob = stats::plogis(eta)
  )

  p_treat <- mean(treat)

  tibble::tibble(
    intercept = intercept, treat = list(treat), p_treat = p_treat
  )

}

# Simulate the treatment -----------------------------------------------
gen_treat <- function(x_matrix, beta, int_grid, prop, delta) {

  assertive::assert_is_matrix(x_matrix)
  assertive::assert_is_numeric(beta)
  assertive::assert_is_numeric(int_grid)
  assertive::assert_is_a_number(prop)
  assertive::assert_is_a_number(delta)

  res <- purrr::map_dfr(
    .x = int_grid,
    ~ prop_treat(
      x_matrix = x_matrix, beta = beta, intercept = .x
    )
  ) %>%
    dplyr::filter(p_treat > prop - delta & p_treat < prop + delta) %>%
    dplyr::slice(nrow(.))

  treat <- unlist(res[["treat"]])

  intercept <- res[["intercept"]]

  list("treat" = treat, "intercept" = intercept)

}

# Compute potential outcomes -------------------------------------------
potential_outcomes <- function(
  x_matrix, treatment, alpha, intercept_outcome, alpha_treatment
) {

  assertive::assert_is_matrix(x_matrix)
  assertive::assert_is_list(treatment)
  assertive::assert_is_numeric(alpha)
  assertive::assert_is_a_number(intercept_outcome)
  assertive::assert_is_a_number(alpha_treatment)

  # Design matrix
  des_mat <- cbind(rep(1, n), treatment$treat, x_matrix)
  des_mat_att <- des_mat[des_mat[, 2] == 1, ]
  coef_1 <- c(intercept_outcome, alpha_treatment, alpha)
  coef_2 <- c(intercept_outcome, 0, alpha)

  # 1) ATT -------------------------------------------------------------
  p1 <- mean(des_mat_att %*% coef_1)
  p0 <- mean(des_mat_att %*% coef_2)

  att_true <- abs(p1 - p0)

  # Outcome
  pp <- des_mat %*% coef_1
  y <- stats::rbinom(n, 1, pp)

  # Results in a final tibble
  tibble::tibble(
    intercept_outcome = intercept_outcome,
    coef_treat = alpha_treatment,
    prop_outcome = mean(y),
    att = att_true
  )

}

# Choose intercept and the beta of treatment of the outcome's model ----
intercept_outcome <- function(
  x_matrix, treatment, alpha, outcome_grid,
  p_outcome, marginal_rd
) {

  assertive::assert_is_matrix(x_matrix)
  assertive::assert_is_list(treatment)
  assertive::assert_is_numeric(alpha)
  assertive::assert_is_numeric(outcome_grid)
  assertive::assert_is_numeric(p_outcome)
  assertive::assert_is_numeric(marginal_rd)

  purrr::map_dfr(
    .x = outcome_grid,
    ~ {

      potential_outcomes(
        x_matrix = x_matrix, treatment = treatment,
        alpha = alpha,
        intercept_outcome = .x, alpha_treatment = marginal_rd
      )

    }
  ) %>%
    # Remove scenarios where P < 0 or P > 1
    dplyr::filter(!is.na(prop_outcome)) %>%
    # Filter by proportion of events
    dplyr::mutate(diff_p = abs(prop_outcome - p_outcome)) %>%
    dplyr::filter(diff_p == min(diff_p)) %>%
    dplyr::select(-diff_p)

}

# Simulate the superpopulation -----------------------------------------
sim_pop <- function(
  n, mus, sds, corr_matrix,
  beta, intercept_treat, p_treat, delta_treat,
  alpha, outcome_grid, p_outcome, marginal_rd
) {

  assertive::assert_is_a_number(n)
  assertive::assert_is_numeric(mus)
  assertive::assert_is_numeric(sds)
  assertive::assert_is_symmetric_matrix(corr_matrix)
  assertive::assert_is_numeric(beta)
  assertive::assert_is_numeric(intercept_treat)
  assertive::assert_is_a_number(p_treat)
  assertive::assert_is_a_number(delta_treat)
  assertive::assert_is_numeric(alpha)
  assertive::assert_is_numeric(outcome_grid)
  assertive::assert_is_numeric(p_outcome)
  assertive::assert_is_numeric(marginal_rd)

  # Generate baseline matrix
  x_matrix <- gen_baseline(
    n = n, mus = mus, sds = sds, corr_matrix = corr_matrix
  )

  # Simulate treatment mechanism
  treatment <- gen_treat(
    x_matrix = x_matrix, beta = beta, int_grid = intercept_treat,
    prop = p_treat, delta = delta_treat
  )

  # Choose intercept and treatment's coefficient of the outcome model
  coef_out <- intercept_outcome(
    x_matrix = x_matrix, treatment = treatment, alpha = alpha,
    outcome_grid = outcome_grid, p_outcome = p_outcome,
    marginal_rd = marginal_rd
  )

  # Simulate the outcome
  des_matrix <- cbind(rep(1, n), treatment$treat, x_matrix)
  alphas <- c(
    coef_out$intercept_outcome, coef_out$coef_treat,
    alpha
  )

  pp <- des_matrix %*% alphas
  y <- stats::rbinom(n, 1, pp)

  dd <- tibble::as_tibble(cbind(y, treatment$treat, x_matrix)) %>%
    dplyr::rename(treatment = V2)

  # Store results into a tibble
  tibble::tibble(
    int_treat = treatment$intercept,
    int_outcome = coef_out$intercept_outcome,
    population = list(dd),
    att_true = coef_out$att
  )
}

# Sample from the superpopulation --------------------------------------
sample_dfs <- function(population, sample_size, n_samples) {

  assertive::assert_is_data.frame(population)
  assertive::assert_is_a_number(sample_size)
  assertive::assert_is_a_number(n_samples)

  # Get the population
  purrr::rerun(
    .n = n_samples, dplyr::sample_n(population, sample_size)
  )

}
