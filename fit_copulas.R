# ==============================================================================
# fit_copulas.R
# Project: The Diversification Question — Copula-Based Aggregation of
#          Operational Risk Losses (NexaCore Financial Technologies)
#
# Purpose: (1) Fits marginal frequency (Poisson) and severity (lognormal)
#              distributions per category from event-level data.
#          (2) Fits five candidate copula families to the ONE pair with
#              genuine dependence — ICT & Business Disruption <-> Execution,
#              Delivery & Process Mgmt — and selects the best fit by AIC.
#
# DECISION (INSTRUCTIONS.md Sec 12, 2026-06-23): the copula stage is
# restricted to this bivariate pair rather than a 5-dimensional joint fit
# across all categories. The other three categories (Internal Fraud,
# External Fraud, Clients & Business Practices) are independent by
# construction in generate_losses.R. A single-parameter 5-D exchangeable
# copula would average the real ICT-Execution dependence against three
# pairs with no true dependence, diluting the effect this project exists to
# demonstrate, and at n=10 annual observations risks fitting sampling noise
# from the unrelated pairs as if it were structure. The three independent
# categories are simulated directly from their own marginal fits in
# joint_simulation.R and added into the joint annual total.
#
# Output:  data/marginal_fits.csv   (Poisson/lognormal params, all 5 categories)
#          data/pseudo_obs.csv      (rank-transformed pair data, for plotting)
#          data/copula_aic.csv      (AIC comparison, `copula` package fits)
#          data/vine_check.csv      (VineCopula::BiCopSelect cross-check)
#          models/best_copula_fit.rds (fitted copula object, for joint_simulation.R)
#
# Seed:    42 (project standard)
# ==============================================================================

library(tidyverse)
library(copula)
library(fitdistrplus)
library(VineCopula)

set.seed(42)

# ------------------------------------------------------------------------------
# 0. Load pre-computed loss data
# ------------------------------------------------------------------------------
loss_events             <- read.csv("data/loss_events.csv",             stringsAsFactors = FALSE)
monthly_category_losses <- read.csv("data/monthly_category_losses.csv", stringsAsFactors = FALSE)
annual_category_losses  <- read.csv("data/annual_category_losses.csv",  stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# 1. Marginal frequency-severity fitting (fitdistrplus)
# ------------------------------------------------------------------------------
# Frequency: Poisson MLE on monthly event counts, per category.
# Severity:  Lognormal MLE on individual event severities, per category.
#
# Note on the Execution/Process category specifically: its pooled monthly
# event-count series includes the contagion-month bursts, so the fitted
# Poisson lambda here reflects the REALIZED (contagion-inflated) marginal
# process, not the "baseline-only" lambda used to generate it. This is the
# correct marginal to fit -- it is what a risk team observing only the
# historical loss series would estimate -- and is noted here so the
# distinction is explicit rather than silently absorbed.

category_names <- unique(annual_category_losses$category)

fit_marginal <- function(cat_name) {

  freq_data <- monthly_category_losses |> dplyr::filter(category == cat_name) |> dplyr::pull(event_count)
  sev_data  <- loss_events             |> dplyr::filter(category == cat_name) |> dplyr::pull(severity)

  freq_fit <- fitdistrplus::fitdist(freq_data, "pois")
  sev_fit  <- fitdistrplus::fitdist(sev_data,  "lnorm")

  tibble::tibble(
    category       = cat_name,
    n_events       = length(sev_data),
    freq_lambda    = unname(freq_fit$estimate["lambda"]),
    freq_lambda_se = unname(freq_fit$sd["lambda"]),
    sev_mu         = unname(sev_fit$estimate["meanlog"]),
    sev_mu_se      = unname(sev_fit$sd["meanlog"]),
    sev_sigma      = unname(sev_fit$estimate["sdlog"]),
    sev_sigma_se   = unname(sev_fit$sd["sdlog"])
  )
}

marginal_fits <- purrr::map_dfr(category_names, fit_marginal)

write.csv(marginal_fits, "data/marginal_fits.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 2. Dependent pair: pivot to wide annual totals, compute pseudo-observations
# ------------------------------------------------------------------------------
dependent_pair <- c("ICT & Business Disruption", "Execution, Delivery & Process Mgmt")

pair_wide <- annual_category_losses |>
  dplyr::filter(category %in% dependent_pair) |>
  dplyr::select(year, category, total_loss) |>
  tidyr::pivot_wider(names_from = category, values_from = total_loss) |>
  dplyr::arrange(year) |>
  dplyr::rename(
    ict  = `ICT & Business Disruption`,
    exec = `Execution, Delivery & Process Mgmt`
  )

# Empirical dependence measures, descriptive (for document narrative —
# these are reported alongside, not in place of, the copula parameter)
kendall_tau  <- cor(pair_wide$ict, pair_wide$exec, method = "kendall")
spearman_rho <- cor(pair_wide$ict, pair_wide$exec, method = "spearman")

# Pseudo-observations: rank-based, scaled to (0,1) -- the standard input
# for maximum pseudo-likelihood copula fitting (Genest et al.)
pseudo_obs <- copula::pobs(as.matrix(pair_wide[, c("ict", "exec")]))
colnames(pseudo_obs) <- c("ict", "exec")

pseudo_obs_df <- tibble::tibble(
  year     = pair_wide$year,
  ict_raw  = pair_wide$ict,
  exec_raw = pair_wide$exec,
  ict_u    = pseudo_obs[, "ict"],
  exec_u   = pseudo_obs[, "exec"]
)
write.csv(pseudo_obs_df, "data/pseudo_obs.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 3. Fit candidate copula families via maximum pseudo-likelihood
# ------------------------------------------------------------------------------
# t-copula degrees of freedom: estimating df by MLE jointly with rho at
# n=10 is unstable (one extra continuous parameter with very little data
# to constrain it). df is therefore selected from a small fixed grid by
# AIC rather than estimated -- a documented simplification, not an
# oversight (each grid value is fit as its own 1-parameter candidate).

fit_one_copula <- function(cop, label, k_params) {
  result <- tryCatch({
    fit <- copula::fitCopula(cop, data = pseudo_obs, method = "mpl")
    ll  <- fit@loglik
    aic <- -2 * ll + 2 * k_params
    list(
      family    = label,
      params    = paste(round(coef(fit), 4), collapse = ", "),
      loglik    = ll,
      aic       = aic,
      converged = TRUE
    )
  }, error = function(e) {
    list(family = label, params = NA_character_, loglik = NA_real_, aic = NA_real_, converged = FALSE)
  })
  tibble::as_tibble(result)
}

fit_t_copula_grid <- function(df_grid = c(2, 4, 6, 10, 20)) {
  purrr::map_dfr(df_grid, function(df_val) {
    cop <- copula::tCopula(dim = 2, df = df_val, df.fixed = TRUE)
    fit_one_copula(cop, paste0("t (df=", df_val, ")"), k_params = 1)
  })
}

candidates <- dplyr::bind_rows(
  fit_one_copula(copula::normalCopula(dim = 2),  "Gaussian", k_params = 1),
  fit_t_copula_grid(),
  fit_one_copula(copula::claytonCopula(dim = 2), "Clayton",  k_params = 1),
  fit_one_copula(copula::gumbelCopula(dim = 2),  "Gumbel",   k_params = 1),
  fit_one_copula(copula::frankCopula(dim = 2),   "Frank",    k_params = 1)
)

copula_aic <- candidates |>
  dplyr::filter(converged) |>
  dplyr::arrange(aic) |>
  dplyr::mutate(rank = dplyr::row_number()) |>
  dplyr::select(rank, family, params, loglik, aic)

best_family <- copula_aic$family[1]

write.csv(copula_aic, "data/copula_aic.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 4. VineCopula cross-check (independent family-selection diagnostic)
# ------------------------------------------------------------------------------
# BiCopSelect() runs its own AIC-based search across a broader family set
# (including rotated Clayton/Gumbel) as a confirmatory check against the
# `copula` package result above -- a diagnostic, not the primary selection.

vine_check <- VineCopula::BiCopSelect(
  pseudo_obs[, "ict"], pseudo_obs[, "exec"],
  familyset = NA, selectioncrit = "AIC"
)

vine_check_summary <- tibble::tibble(
  vine_family_name = vine_check$familyname,
  vine_family_code = vine_check$family,
  vine_par         = vine_check$par,
  vine_par2        = vine_check$par2,
  vine_aic         = vine_check$AIC
)

write.csv(vine_check_summary, "data/vine_check.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 5. Save best-fit copula object for joint_simulation.R
# ------------------------------------------------------------------------------
dir.create("models", showWarnings = FALSE)

family_token <- stringr::str_extract(best_family, "^[A-Za-z]+")

best_cop_object <- switch(
  family_token,
  "Gaussian" = copula::normalCopula(dim = 2),
  "t"        = copula::tCopula(dim = 2, df = as.numeric(stringr::str_extract(best_family, "\\d+")), df.fixed = TRUE),
  "Clayton"  = copula::claytonCopula(dim = 2),
  "Gumbel"   = copula::gumbelCopula(dim = 2),
  "Frank"    = copula::frankCopula(dim = 2)
)

best_fit <- copula::fitCopula(best_cop_object, data = pseudo_obs, method = "mpl")
saveRDS(best_fit, "models/best_copula_fit.rds")

# ------------------------------------------------------------------------------
# 6. Console diagnostics
# ------------------------------------------------------------------------------
cat("=== fit_copulas.R: diagnostic summary ===\n\n")

cat("Marginal fits (frequency-severity, all 5 categories):\n")
print(marginal_fits)

cat("\nDependent pair -- annual totals (n =", nrow(pair_wide), "years):\n")
print(pair_wide)

cat("\nEmpirical Kendall's tau:  ", round(kendall_tau, 4), "\n", sep = "")
cat("Empirical Spearman's rho: ", round(spearman_rho, 4), "\n\n", sep = "")

cat("Copula AIC comparison (`copula` package, ranked):\n")
print(copula_aic)

cat("\nSelected family:", best_family, "\n\n")

cat("VineCopula cross-check (BiCopSelect):\n")
print(vine_check_summary)

cat("\nConverged candidates:", nrow(copula_aic), "of", nrow(candidates), "fitted\n")

cat("\nDone.\n")
