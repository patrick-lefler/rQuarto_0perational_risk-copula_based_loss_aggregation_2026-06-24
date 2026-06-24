# ==============================================================================
# generate_losses.R
# Project: The Diversification Question — Copula-Based Aggregation of
#          Operational Risk Losses (NexaCore Financial Technologies)
#
# Purpose: Generates 10 years (120 months) of synthetic operational loss
#          event data across 5 Basel-style event-type categories using a
#          compound Poisson-lognormal frequency-severity framework.
#
#          A deliberate, documented dependence structure is embedded between
#          ICT & Business Disruption and Execution, Delivery & Process
#          Management: months containing an ICT disruption event trigger an
#          elevated frequency AND severity of process-failure losses in the
#          same month, simulating operational contagion from a systems
#          outage into downstream processing errors. This is the dependence
#          signal the copula-fitting stage (fit_copulas.R) is designed to
#          recover and quantify.
#
#          All other category pairs are simulated independently — the
#          contagion link is intentionally isolated to one pair so that the
#          copula-selection stage has a clear, attributable dependence
#          source to detect rather than diffuse cross-category noise.
#
# Output:  data/loss_events.csv             (event-level detail)
#          data/monthly_category_losses.csv (monthly category aggregates)
#          data/annual_category_losses.csv  (annual category aggregates —
#                                             primary input for copula fitting
#                                             in fit_copulas.R)
#
# Seed:    42 (project standard, set.seed(42))
#
# Decision log ref: INSTRUCTIONS.md §12 — severity-tail depth intentionally
# simplified (lognormal only, no GPD/POT) to keep scope on the aggregation
# question. See "The Tail You Cannot Fit" for tail-fitting methodology.
# ==============================================================================

library(tidyverse)

set.seed(42)

# ------------------------------------------------------------------------------
# 1. Category parameters
# ------------------------------------------------------------------------------
# Frequency: annual Poisson lambda (events/year), converted to monthly lambda.
# Severity:  lognormal, parameterized by target MEDIAN severity and sigma
#            (mu = log(median); median chosen over mean for interpretability
#            when communicating assumptions to a non-technical reader).
#
# NOTE: All severity figures are illustrative synthetic values, not derived
# from any real NexaCore data or external loss database. Currency unit: USD,
# per the project-wide number-formatting standard (INSTRUCTIONS.md §9).

categories <- tibble::tribble(
  ~category,                             ~lambda_annual, ~severity_median, ~severity_sigma,
  "Internal Fraud",                      8,              50000,            1.2,
  "External Fraud",                      40,             2000,             0.9,
  "ICT & Business Disruption",           15,             80000,            1.0,
  "Execution, Delivery & Process Mgmt",  60,             5000,             0.8,
  "Clients & Business Practices",        10,             150000,           1.3
) |>
  dplyr::mutate(
    lambda_monthly = lambda_annual / 12,
    severity_mu    = log(severity_median)
  )

# ------------------------------------------------------------------------------
# 2. Contagion design (documented misconfiguration)
# ------------------------------------------------------------------------------
# In any month where ICT & Business Disruption generates >= 1 event:
#   - Execution, Delivery & Process Mgmt draws ADDITIONAL events that month
#     from a separate Poisson(contagion_lambda_extra)
#   - ALL Execution/Process severities in that month (baseline + contagion
#     events) are scaled by contagion_severity_multiplier, representing
#     remediation effort and processing-backlog costs that compound during
#     an outage
#
# This creates loss co-movement between the two categories without altering
# their marginal frequency/severity parameters in non-contagion months — the
# dependence is structural and time-localized, not a marginal artifact. That
# is the property the copula-fitting stage is meant to detect and quantify.

contagion_lambda_extra        <- 3
contagion_severity_multiplier <- 1.4

n_years  <- 10
n_months <- n_years * 12

# ------------------------------------------------------------------------------
# 3. Simulate ICT disruption event counts first (drives contagion timing)
# ------------------------------------------------------------------------------
ict_params <- categories |> dplyr::filter(category == "ICT & Business Disruption")

ict_monthly_counts <- rpois(n_months, lambda = ict_params$lambda_monthly)
contagion_months   <- which(ict_monthly_counts > 0)

# ------------------------------------------------------------------------------
# 4. Event-level simulation
# ------------------------------------------------------------------------------
simulate_category_events <- function(cat_row, month_idx, n_events, severity_mult = 1) {
  if (n_events == 0) return(NULL)
  tibble::tibble(
    category    = cat_row$category,
    month_index = month_idx,
    severity    = rlnorm(n_events, meanlog = cat_row$severity_mu, sdlog = cat_row$severity_sigma) * severity_mult
  )
}

event_list <- vector("list", length = n_months * nrow(categories))
counter <- 1

for (m in seq_len(n_months)) {
  for (i in seq_len(nrow(categories))) {

    cat_row <- categories[i, ]

    if (cat_row$category == "ICT & Business Disruption") {

      n_events <- ict_monthly_counts[m]
      events   <- simulate_category_events(cat_row, m, n_events)

    } else if (cat_row$category == "Execution, Delivery & Process Mgmt") {

      n_baseline   <- rpois(1, lambda = cat_row$lambda_monthly)
      is_contagion <- m %in% contagion_months
      n_contagion  <- if (is_contagion) rpois(1, lambda = contagion_lambda_extra) else 0
      sev_mult     <- if (is_contagion) contagion_severity_multiplier else 1

      events <- simulate_category_events(cat_row, m, n_baseline + n_contagion, sev_mult)

    } else {

      n_events <- rpois(1, lambda = cat_row$lambda_monthly)
      events   <- simulate_category_events(cat_row, m, n_events)
    }

    event_list[[counter]] <- events
    counter <- counter + 1
  }
}

loss_events <- dplyr::bind_rows(event_list) |>
  dplyr::arrange(month_index, category) |>
  dplyr::mutate(
    event_id = dplyr::row_number(),
    year     = ((month_index - 1) %/% 12) + 1,
    month    = ((month_index - 1) %% 12) + 1
  ) |>
  dplyr::select(event_id, year, month, month_index, category, severity)

# ------------------------------------------------------------------------------
# 5. Aggregations
# ------------------------------------------------------------------------------
monthly_category_losses <- loss_events |>
  dplyr::group_by(category, month_index) |>
  dplyr::summarise(
    event_count = dplyr::n(),
    total_loss  = sum(severity),
    .groups = "drop"
  ) |>
  # complete category x month_index grid so months with zero losses are
  # represented explicitly rather than silently dropped
  tidyr::complete(
    category, month_index = 1:n_months,
    fill = list(event_count = 0, total_loss = 0)
  ) |>
  dplyr::mutate(
    year  = ((month_index - 1) %/% 12) + 1,
    month = ((month_index - 1) %% 12) + 1
  ) |>
  dplyr::arrange(category, month_index) |>
  dplyr::select(category, year, month, month_index, event_count, total_loss)

annual_category_losses <- monthly_category_losses |>
  dplyr::group_by(category, year) |>
  dplyr::summarise(
    event_count = sum(event_count),
    total_loss  = sum(total_loss),
    .groups = "drop"
  ) |>
  dplyr::arrange(category, year)

# ------------------------------------------------------------------------------
# 6. Write outputs
# ------------------------------------------------------------------------------
dir.create("data", showWarnings = FALSE)

write.csv(loss_events,             "data/loss_events.csv",             row.names = FALSE)
write.csv(monthly_category_losses, "data/monthly_category_losses.csv", row.names = FALSE)
write.csv(annual_category_losses,  "data/annual_category_losses.csv",  row.names = FALSE)

# ------------------------------------------------------------------------------
# 7. Console diagnostics
# ------------------------------------------------------------------------------
# Per project pattern: Quarto's message/warning suppression can mask
# pipeline issues, so this script reports row counts and summary stats
# directly to console when run standalone (Rscript scripts/generate_losses.R)
# or sourced from pipeline_diagnostic.R.

cat("=== generate_losses.R: diagnostic summary ===\n\n")

cat("Event-level rows:        ", nrow(loss_events), "\n")
cat("Monthly aggregate rows:  ", nrow(monthly_category_losses),
    " (expected: ", n_months * nrow(categories), ")\n", sep = "")
cat("Annual aggregate rows:   ", nrow(annual_category_losses),
    " (expected: ", n_years * nrow(categories), ")\n\n", sep = "")

cat("Contagion months triggered:", length(contagion_months), "of", n_months, "total months\n\n")

cat("Event counts by category:\n")
print(
  loss_events |>
    dplyr::count(category, name = "n_events") |>
    dplyr::arrange(dplyr::desc(n_events))
)

cat("\nAnnual total loss summary by category (mean / sd across 10 years):\n")
print(
  annual_category_losses |>
    dplyr::group_by(category) |>
    dplyr::summarise(
      mean_annual_loss = mean(total_loss),
      sd_annual_loss   = sd(total_loss),
      .groups = "drop"
    )
)

cat("\nDone.\n")
