# Absolute model contribution (per-variable models vs intercept-only baseline)
## Assumes objects from modelling script exist: euSD.dt2, clc_idx, bioG_idx, AN2, euSD.spde2

###############
pc <- list(prec = list(prior = "pc.prec", param = c(0.1, 0.01)))

###############
## INLA controls

control.fixed <- list(mean.intercept = 0, prec.intercept = 1e-3)

control.inla <- list(
  diagonal = 1e-4,
  h        = 0.1,
  tolerance= 1e-5,
  strategy = "gaussian",
  improved.simplified.laplace = TRUE
)

control.compute <- list(
  dic    = TRUE,
  waic   = FALSE,
  cpo    = FALSE,
  config = TRUE,
  return.marginals.predictor = FALSE
)

control.family <- list(link = "log")

###############
## Absolute contribution models

f1  <- L ~ 1
f2  <- L ~ 1 + NBud + NDep
f3  <- L ~ 1 + f(clcL2, model = "iid", constr = TRUE, hyper = pc)
f4  <- L ~ 1 + f(bioG,  model = "iid", constr = TRUE,
                 hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))),
                 diagonal = 1e-6)
f5  <- L ~ 1 + f(area_g, model = "rw1", scale.model = TRUE, hyper = pc)
f6  <- L ~ 1 + f(year,   model = "rw2", hyper = pc)
f7  <- L ~ 1 + f(pH_g,   model = "rw1", scale.model = TRUE, hyper = pc, diagonal = 1e-6)
f8  <- L ~ 1 + f(temp_g, model = "rw1", scale.model = TRUE, hyper = pc)
f9 <- L ~ 1 + f(prep_g, model = "rw1", scale.model = TRUE, hyper = pc)
f10 <- L ~ 1 + f(elev_g, model = "rw1", scale.model = TRUE, hyper = pc)
f11 <- L ~ 1 + f(nodes,  model = euSpde2)

formulas <- list(f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11)
formula_names <- c("base", "fixed_N", "clcL2", "bioG", "area", "year", "pH", "temp", "prec", "altitude", "space")

###############
## Fit loop

Sys.setenv(
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS      = "1",
  OMP_NUM_THREADS      = "3"
)
inla.setOption(num.threads = "3:1")
set.seed(42)

resultsABS <- lapply(seq_along(formulas), function(i) {
  cat("\n▶ Running model:", formula_names[i], "(", i, "of", length(formulas), ")\n")

  tryCatch({
    mod <- inla(
      formulas[[i]],
      family            = "gamma",
      data              = inla.stack.data(stk2),
      control.family    = control.family,
      control.predictor = list(A = inla.stack.A(stk2), compute = FALSE, link = 1),
      control.fixed     = control.fixed,
      control.inla      = control.inla,
      control.compute   = control.compute
    )

    cat("   ✓ Model", formula_names[i], "fit completed. Rerunning...\n")
    modR <- inla.rerun(mod)
    print(summary(modR))
    modR
  }, error = function(e) {
    cat("Model", formula_names[i], "failed:\n   ", conditionMessage(e), "\n")
    NULL
  })
})

names(resultsABS) <- formula_names

###############
## Log-evidence differences vs base

logZ_base <- resultsABS[["base"]]$mlik[1]

abs_evidence <- sapply(resultsABS, function(mod) {
  round(mod$mlik[1] - logZ_base, 2)
})

abs_evidence_df <- data.frame(
  model = formula_names,
  log_evidence_diff = abs_evidence
)

dput(abs_evidence_df)
