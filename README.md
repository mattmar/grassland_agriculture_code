# Grassland–Agriculture Modelling Code

This repository contains the R scripts used to reproduce the statistical analyses for the study
“Surrounding low-intensity agriculture benefits grassland plant diversity across Europe”.

The workflow implements a spatially explicit Bayesian framework using INLA with a barrier-SPDE
to quantify the effects of agricultural intensity, nitrogen inputs, climate, land use, and
spatial structure on grassland plant diversity across Europe.

--------------------------------------------------------------------
Repository structure
--------------------------------------------------------------------

scripts/
├── 01_inla_model-part01.R
├── 02_inla_model-part02.R
├── 03_figures-part01.R
├── 04_figures-part02.R
└── 05_inla_model-absolute-contribution.R
└── 06_inla_model-marginal-contribution.R

--------------------------------------------------------------------
Script overview
--------------------------------------------------------------------

PART 01 — Data screening and filtering
--------------------------------------------------
01_inla_model-part01.R

• Loads vegetation plot data and environmental covariates.
• Harmonises land-use classes and survey years.
• Applies exclusion criteria based on spatial uncertainty.
• Fits a preliminary INLA model without nitrogen effects.
• Identifies and removes observations inconsistent with theory
• Outputs indices used to subset data in downstream analyses.

--------------------------------------------------

PART 02 — Main spatial INLA model (paper results)
--------------------------------------------------
02_inla_model-part02.R

• Fits the full Gamma regression model with log link.
• Includes fixed effects for nitrogen fertiliser (NBud) and
  nitrogen deposition (NDep).
• Includes random effects for:
  – land-cover class (CLC level 2),
  – biogeographical region,
  – survey year,
  – plot area,
  – soil pH,
  – temperature,
  – precipitation,
  – elevation.
• Models spatial autocorrelation using a barrier-SPDE.
• Uses penalised complexity (PC) priors throughout.
• Back-transforms nitrogen coefficients to the original scale
  and reports percentage changes in diversity.

--------------------------------------------------

PART 03 — Main figures
--------------------------------------------------
03_figures-part01.R

• Generates all main manuscript figures.
• Includes:
  – raw data visualisations,
  – model-based response curves for nitrogen inputs,
  – land-use effects,
  – spatial predictions and maps.
• All figures are produced directly from the fitted INLA object
  generated in Part 02.

--------------------------------------------------

PART 04 — Random-effects visualisation
--------------------------------------------------
04_figures-part02.R

• Visualises marginal posterior effects of random terms only.
• Includes plots for:
  – biogeographical regions,
  – year,
  – plot area,
  – climate variables,
  – soil pH,
  – elevation.
• All effects are transformed back to the response scale
  (relative change in diversity).

--------------------------------------------------

PART 05 — Model absolute contribution analyses
--------------------------------------------------
05_inla_model-absolute-contribution.R

• Assesses the contribution of each model component using:
  – absolute (single-term) model comparisons.
• Computes DIC, marginal likelihood differences, and Bayesian R².

PART 06 — Model marginal contribution analyses
--------------------------------------------------------------------
06_inla_model-marginal-contribution.R

• Assesses the contribution of each model component using:
  – marginal (stepwise) log-evidence differences,

REQUIREMENTS
--------------------------------------------------------------------
• R version ≥ 4.3
• Core packages:
  INLA, inlabru, fmesher, sf, terra, sp, ggplot2, mgcv, parallel

--------------------------------------------------------------------
Usage
--------------------------------------------------------------------

Scripts are designed to be run sequentially:
All scripts assume that local paths are set at the top of each file
using the variable `baseDir`.
