# ==============================================================================
# joint_simulation.R
# Project: The Diversification Question — Copula-Based Aggregation of
#          Operational Risk Losses (NexaCore Financial Technologies)
#
# Purpose: (1) Monte Carlo simulation of the joint annual total loss across
#              all 5 categories: the ICT <-> Execution/Process pair is
#              simulated through each candidate copula (Sklar's theorem —
#              copula-drawn uniforms mapped through each category's own
#              empirical marginal annual-loss distribution), while the
#              other 3 categories are simulated independently and summed in.
#          (2) Computes diversified 99.9% VaR and Expected Shortfall (ES) on
#              the joint distribution, and compares against the additive
#              ("silo") sum of standalone category VaR/ES.
#          (3) Repeats (1) across the 5 candidate copula families (AIC-best
#              parameter per family) to show how much the diversification
#              benefit estimate moves with the dependence assumption.
#
# DESIGN NOTE: a single rCopula() draw per family feeds BOTH the headline
# "primary" report and the sensitivity table -- the primary report is the
# row of the sensitivity table matching the AIC-selected family (t, df=2),
# not a separately re-simulated scenario. An earlier version of this script
# drew the primary scenario independently from the sensitivity loop, which
# produced two slightly different Monte Carlo realizations of the identical
# t(df=2) model (~0.8% apart) reported in two different places downstream.
# This version uses one Monte Carlo draw per family as the single source of
# truth for both.
#
# DECISION (INSTRUCTIONS.md Sec 12, 2026-06-23): the primary dependence
# assumption is t(df=2), the best AIC fit among the five families this
# project set out to compare. A VineCopula cross-check found a Tawn
# extreme-value copula fits marginally better; that finding is documented
# in the rendered report as a scope limitation (Sec 3.3.3) and is NOT
# included in the sensitivity sweep below, which is restricted to the five
# originally-specified families to keep the comparison apples-to-apples.
#
# Output:  data/standalone_var_es.csv  (per-category standalone VaR/ES)
#          data/joint_sim.csv          (50,000-draw joint simulation, primary copula)
#          data/var_comparison.csv     (diversified vs. additive VaR/ES, by copula family)
#
# Seed:    42 (project standard)
# ==============================================================================

library(tidyverse)
library(copula)

set.seed(42)

# ------------------------------------------------------------------------------
# 0. Load pre-computed inputs
# ------------------------------------------------------------------------------
marginal_fits <- read.csv("data/marginal_fits.csv", stringsAsFactors = FALSE)
copula_aic    <- read.csv("data/copula_aic.csv",    stringsAsFactors = FALSE)

n_sim     <- 50000
var_level <- 0.999

independent_categories <- c("Internal Fraud", "External Fraud", "Clients & Business Practices")

# ------------------------------------------------------------------------------
# 1. Marginal annual loss simulator (compound Poisson-lognormal)
# ------------------------------------------------------------------------------
# Vectorized: draws n_sim annual event counts, then a single pooled vector
# of severities, and uses rowsum() to aggregate severities back to their
# originating iteration. Avoids an explicit per-iteration R loop, which
# would be materially slower at n_sim = 50,000 for the higher-frequency
# categories (Execution/Process: ~870 events/year-equivalent pool).

simulate_annual_losses <- function(n_sim, lambda_monthly, sev_mu, sev_sigma) {
  lambda_annual <- lambda_monthly * 12
  n_events      <- rpois(n_sim, lambda_annual)
  total_events  <- sum(n_events)

  losses <- numeric(n_sim)
  if (total_events == 0) return(losses)

  severities <- rlnorm(total_events, meanlog = sev_mu, sdlog = sev_sigma)
  idx        <- rep(seq_len(n_sim), times = n_events)

  grouped <- rowsum(severities, group = idx)
  losses[as.integer(rownames(grouped))] <- grouped[, 1]
  losses
}

get_marginal_params <- function(cat_name) {
  row <- marginal_fits |> dplyr::filter(category == cat_name)
  list(lambda_monthly = row$freq_lambda, sev_mu = row$sev_mu, sev_sigma = row$sev_sigma)
}

# ------------------------------------------------------------------------------
# 2. Simulate the 3 independent categories
# ------------------------------------------------------------------------------
independent_losses <- purrr::map(independent_categories, function(cat_name) {
  p <- get_marginal_params(cat_name)
  simulate_annual_losses(n_sim, p$lambda_monthly, p$sev_mu, p$sev_sigma)
})
names(independent_losses) <- independent_categories

# ------------------------------------------------------------------------------
# 3. Marginal annual loss samples for the dependent pair
# ------------------------------------------------------------------------------
# These samples serve double duty: (a) standalone marginal distributions
# for the additive/silo VaR-ES calculation, and (b) the empirical
# inverse-CDF target for the copula-uniform-to-loss mapping in Step 6.
# Reusing the same samples across both uses, and across every copula
# family tested in Step 6, isolates the effect of the dependence
# assumption -- only the joint structure changes; the marginals do not.

ict_params  <- get_marginal_params("ICT & Business Disruption")
exec_params <- get_marginal_params("Execution, Delivery & Process Mgmt")

ict_marginal_sample  <- simulate_annual_losses(n_sim, ict_params$lambda_monthly,  ict_params$sev_mu,  ict_params$sev_sigma)
exec_marginal_sample <- simulate_annual_losses(n_sim, exec_params$lambda_monthly, exec_params$sev_mu, exec_params$sev_sigma)

# ------------------------------------------------------------------------------
# 4. Helpers
# ------------------------------------------------------------------------------
# build_copula_object: constructs a copula object from a (family label,
# parameter) pair as written by fit_copulas.R, e.g. "Gaussian", "t (df=2)",
# "Clayton", "Gumbel", "Frank". Every candidate has exactly one free
# parameter (t's df is fixed, not estimated -- see fit_copulas.R Sec 3).

build_copula_object <- function(family_label, param_value) {
  family_token <- stringr::str_extract(family_label, "^[A-Za-z]+")
  switch(
    family_token,
    "Gaussian" = copula::normalCopula(param = param_value, dim = 2),
    "t"        = copula::tCopula(param = param_value, dim = 2,
                                  df = as.numeric(stringr::str_extract(family_label, "\\d+")),
                                  df.fixed = TRUE),
    "Clayton"  = copula::claytonCopula(param = param_value, dim = 2),
    "Gumbel"   = copula::gumbelCopula(param = param_value, dim = 2),
    "Frank"    = copula::frankCopula(param = param_value, dim = 2)
  )
}

# empirical_quantile_map: maps a vector of copula-drawn uniforms to the
# corresponding quantile of an empirical sample of equal length, by rank.
# Exact (no interpolation), and fast at this scale.

empirical_quantile_map <- function(u, sample_draws) {
  sorted <- sort(sample_draws)
  n      <- length(sorted)
  idx    <- pmin(pmax(ceiling(u * n), 1), n)
  sorted[idx]
}

standalone_var <- function(x, p = var_level) unname(quantile(x, p, type = 7))
standalone_es  <- function(x, p = var_level) {
  v <- standalone_var(x, p)
  mean(x[x >= v])
}

# ------------------------------------------------------------------------------
# 5. Standalone (silo) VaR / ES per category, and the additive sum
# ------------------------------------------------------------------------------
all_category_samples <- list(
  `Internal Fraud`                     = independent_losses[["Internal Fraud"]],
  `External Fraud`                     = independent_losses[["External Fraud"]],
  `ICT & Business Disruption`          = ict_marginal_sample,
  `Execution, Delivery & Process Mgmt` = exec_marginal_sample,
  `Clients & Business Practices`       = independent_losses[["Clients & Business Practices"]]
)

standalone_var_es <- tibble::tibble(
  category = names(all_category_samples),
  var_999  = purrr::map_dbl(all_category_samples, standalone_var),
  es_999   = purrr::map_dbl(all_category_samples, standalone_es)
)

additive_var <- sum(standalone_var_es$var_999)
additive_es  <- sum(standalone_var_es$es_999)

write.csv(standalone_var_es, "data/standalone_var_es.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 6. Joint simulation across the 5 representative copula families
# ------------------------------------------------------------------------------
# One Monte Carlo draw per family. The Tawn finding from fit_copulas.R is
# deliberately excluded here -- see DECISION note above.

representative_families <- copula_aic |>
  dplyr::mutate(family_group = stringr::str_extract(family, "^[A-Za-z]+")) |>
  dplyr::group_by(family_group) |>
  dplyr::slice_min(aic, n = 1) |>
  dplyr::ungroup() |>
  dplyr::arrange(aic)

joint_total_by_family <- list()

var_comparison <- purrr::map_dfr(seq_len(nrow(representative_families)), function(i) {

  row         <- representative_families[i, ]
  param_value <- as.numeric(row$params)
  cop         <- build_copula_object(row$family, param_value)

  u           <- copula::rCopula(n_sim, cop)
  ict_loss_i  <- empirical_quantile_map(u[, 1], ict_marginal_sample)
  exec_loss_i <- empirical_quantile_map(u[, 2], exec_marginal_sample)

  joint_total_i <- ict_loss_i + exec_loss_i +
    independent_losses[["Internal Fraud"]] +
    independent_losses[["External Fraud"]] +
    independent_losses[["Clients & Business Practices"]]

  joint_total_by_family[[row$family]] <<- list(ict = ict_loss_i, exec = exec_loss_i, total = joint_total_i)

  v <- standalone_var(joint_total_i)
  e <- standalone_es(joint_total_i)

  tibble::tibble(
    family          = row$family,
    aic             = row$aic,
    diversified_var = v,
    diversified_es  = e,
    additive_var    = additive_var,
    additive_es     = additive_es,
    benefit_var     = additive_var - v,
    benefit_es      = additive_es  - e,
    benefit_var_pct = (additive_var - v) / additive_var,
    benefit_es_pct  = (additive_es  - e) / additive_es
  )
}) |>
  dplyr::arrange(aic)

write.csv(var_comparison, "data/var_comparison.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 7. Primary ("headline") scenario -- read directly off the table above
# ------------------------------------------------------------------------------
primary_label  <- copula_aic$family[copula_aic$rank == 1]
primary_param  <- as.numeric(copula_aic$params[copula_aic$rank == 1])
primary_result <- var_comparison |> dplyr::filter(family == primary_label)
primary_draws  <- joint_total_by_family[[primary_label]]

diversified_var_primary <- primary_result$diversified_var
diversified_es_primary  <- primary_result$diversified_es
benefit_var_primary     <- primary_result$benefit_var
benefit_es_primary      <- primary_result$benefit_es

write.csv(
  tibble::tibble(
    iter           = seq_len(n_sim),
    internal_fraud = independent_losses[["Internal Fraud"]],
    external_fraud = independent_losses[["External Fraud"]],
    ict            = primary_draws$ict,
    exec           = primary_draws$exec,
    conduct        = independent_losses[["Clients & Business Practices"]],
    joint_total    = primary_draws$total
  ),
  "data/joint_sim.csv", row.names = FALSE
)

# ------------------------------------------------------------------------------
# 8. Console diagnostics
# ------------------------------------------------------------------------------
cat("=== joint_simulation.R: diagnostic summary ===\n\n")

cat("Monte Carlo iterations:", n_sim, "| VaR/ES level:", var_level, "\n\n")

cat("Standalone VaR/ES by category:\n")
print(standalone_var_es)

cat("\nAdditive (silo) VaR_999: ", formatC(additive_var, format = "f", big.mark = ",", digits = 0), "\n", sep = "")
cat("Additive (silo) ES_999:  ",  formatC(additive_es,  format = "f", big.mark = ",", digits = 0), "\n\n", sep = "")

cat("Primary copula family:", primary_label, "| param:", primary_param, "\n")
cat("Diversified VaR_999:  ", formatC(diversified_var_primary, format = "f", big.mark = ",", digits = 0), "\n", sep = "")
cat("Diversified ES_999:   ", formatC(diversified_es_primary,  format = "f", big.mark = ",", digits = 0), "\n\n", sep = "")

cat("Diversification benefit (VaR): ", formatC(benefit_var_primary, format = "f", big.mark = ",", digits = 0),
    " (", round(100 * benefit_var_primary / additive_var, 1), "%)\n", sep = "")
cat("Diversification benefit (ES):  ", formatC(benefit_es_primary, format = "f", big.mark = ",", digits = 0),
    " (", round(100 * benefit_es_primary / additive_es, 1), "%)\n\n", sep = "")

cat("Sensitivity across copula families (marginals fixed, one draw per family):\n")
print(var_comparison)

cat("\nSpread across families -- VaR benefit %: ",
    round(100 * min(var_comparison$benefit_var_pct), 1), "% to ",
    round(100 * max(var_comparison$benefit_var_pct), 1), "%\n", sep = "")
cat("Spread across families -- ES benefit %:  ",
    round(100 * min(var_comparison$benefit_es_pct), 1), "% to ",
    round(100 * max(var_comparison$benefit_es_pct), 1), "%\n\n", sep = "")

cat("Done.\n")
