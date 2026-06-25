# The Diversification Question: Copula-Based Aggregation of Operational Risk Losses

Author: Patrick Lefler
Published: 2026-06-25
Rendered: https://patrick-lefler.github.io/rQuarto_0perational_risk-copula_based_loss_aggregation_2026-06-24/

## Project Introduction
> Tests whether operational risk capital built from a silo sum overstates the true number, and whether that answer depends on which dependence model you trust.

## Overview
This project models operational losses across five categories for NexaCore Financial Technologies using a compound Poisson-lognormal framework, then fits a Student-t copula to the one category pair carrying a genuine, deliberately engineered dependence relationship. A 50,000-iteration Monte Carlo simulation of the resulting joint loss distribution produces a diversified 99.9% Value-at-Risk and Expected Shortfall, compared against the additive sum of standalone category risk measures. The headline result is a real diversification benefit of roughly 30%, and a more interesting finding that this benefit barely moves across five different copula families — because most of the loss base is independent by construction, not because any one dependence model is uniquely correct. The intended outcome is a defensible internal economic capital framework, explicitly not a regulatory submission.

## Tech Stack
* **Language:** R
* **Framework:** [Quarto](https://quarto.org/)
* **Primary Libraries:** tidyverse, copula, fitdistrplus, VineCopula, kableExtra, scales, showtext, sysfonts, sessioninfo
* **Deployment/Output:** Self-contained HTML report (`embed-resources: true`)

## Repository Structure
```
copula-operational-risk/
├── data/                       # CSV outputs from pre-computation scripts
├── scripts/
│   ├── generate_losses.R       # Synthetic frequency-severity simulation
│   ├── fit_copulas.R           # Pseudo-observation transform + copula MLE/AIC
│   └── joint_simulation.R      # Joint loss simulation + VaR/ES calculation
├── models/                     # Saved copula fit object (.rds)
├── output/                     # Rendered HTML
├── _brand.yml
├── _quarto.yml
├── INSTRUCTIONS.md
└── index.qmd
```

## Key Findings
> 1. Diversifying across operational loss categories, rather than summing standalone 99.9% risk measures, reduces the capital estimate by 29.5% (VaR) and 33.8% (Expected Shortfall).
> 2. That benefit is strikingly insensitive to the choice of copula family — a roughly one-percentage-point spread across five candidates — because only one of five categories carries any modeled dependence at all. Most of the benefit comes from summing independent risks, not from the specific dependence model chosen.
> 3. A goodness-of-fit cross-check found a better-fitting extreme-value (Tawn) copula outside the five tested families, and one candidate (Clayton) converged to an implausible parameter sign. Both are reported transparently as scope-bounded limitations rather than smoothed over.

## License
This project is licensed under the MIT License. See the LICENSE file for details.

## Contact
Patrick Lefler [https://www.linkedin.com/in/patricklefler/] | [patrick-lefler.github.io] | [https://substack.com/@pflefler]
