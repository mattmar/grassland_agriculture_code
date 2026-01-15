# Assess marginal contribution of covariates (stepwise log-evidence / DIC / CPO)
## Assumes objects from modelling script exist: euSD.dt2, clc_idx, bioG_idx, AN2, euSD.spde2

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
## Stepwise formulas

m1  <- L ~ 1
m2  <- L ~ 1 + NBud + NDep

m3  <- L ~ 1 + NBud + NDep +
  f(clcL2, model = "iid", constr = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))),
    diagonal = 1e-6)

m4  <- L ~ 1 + NBud + NDep +
  f(clcL2, model = "iid", constr = TRUE, hyper = pc) +
  f(bioG,  model = "iid", constr = TRUE, hyper = pc)

m5  <- L ~ 1 + NBud + NDep +
  f(clcL2, model = "iid", constr = TRUE, hyper = pc) +
  f(bioG,  model = "iid", constr = TRUE, hyper = pc) +
  f(area_g, model = "rw1", scale.model = TRUE, hyper = pc)

m6  <- L ~ 1 + NBud + NDep +
  f(clcL2, model = "iid", constr = TRUE, hyper = pc) +
  f(bioG,  model = "iid", constr = TRUE, hyper = pc) +
  f(area_g, model = "rw1", scale.model = TRUE, hyper = pc) +
  f(year,   model = "rw2", hyper = pc)

m7  <- L ~ 1 + NBud + NDep +
  f(clcL2, model = "iid", constr = TRUE, hyper = pc) +
  f(bioG,  model = "iid", constr = TRUE, hyper = pc) +
  f(area_g, model = "rw1", scale.model = TRUE, hyper = pc) +
  f(year,   model = "rw1", hyper = pc) +
  f(pH_g,   model = "rw2", scale.model = TRUE, hyper = pc, diagonal = 1e-6)

m8  <- L ~ 1 + NBud + NDep +
  f(clcL2, model = "iid", constr = TRUE, hyper = pc) +
  f(bioG,  model = "iid", constr = TRUE, hyper = pc) +
  f(area_g, model = "rw1", scale.model = TRUE, hyper = pc) +
  f(year,   model = "rw2", hyper = pc) +
  f(pH_g,   model = "rw1", scale.model = TRUE, hyper = pc, diagonal = 1e-6) +
  f(temp_g, model = "rw1", scale.model = TRUE, hyper = pc)

m9  <- L ~ 1 + NBud + NDep +
  f(clcL2, model = "iid", constr = TRUE, hyper = pc) +
  f(bioG,  model = "iid", constr = TRUE, hyper = pc) +
  f(area_g, model = "rw1", scale.model = TRUE, hyper = pc) +
  f(year,   model = "rw2", hyper = pc) +
  f(pH_g,   model = "rw1", scale.model = TRUE, hyper = pc, diagonal = 1e-6) +
  f(temp_g, model = "rw1", scale.model = TRUE, hyper = pc) +
  f(prep_g, model = "rw1", scale.model = TRUE, hyper = pc)

m10 <- L ~ 1 + NBud + NDep +
  f(clcL2, model = "iid", constr = TRUE, hyper = pc) +
  f(bioG,  model = "iid", constr = TRUE, hyper = pc) +
  f(area_g, model = "rw1", scale.model = TRUE, hyper = pc) +
  f(year,   model = "rw2", hyper = pc) +
  f(pH_g,   model = "rw1", scale.model = TRUE, hyper = pc, diagonal = 1e-6) +
  f(temp_g, model = "rw1", scale.model = TRUE, hyper = pc) +
  f(prep_g, model = "rw1", scale.model = TRUE, hyper = pc) +
  f(elev_g, model = "rw1", scale.model = TRUE, hyper = pc)

m11 <- L ~ 1 + NBud + NDep +
  f(clcL2, model = "iid", constr = TRUE, hyper = pc) +
  f(bioG,  model = "iid", constr = TRUE, hyper = pc) +
  f(area_g, model = "rw1", scale.model = TRUE, hyper = pc) +
  f(year,   model = "rw2", hyper = pc) +
  f(pH_g,   model = "rw1", scale.model = TRUE, hyper = pc, diagonal = 1e-6) +
  f(temp_g, model = "rw1", scale.model = TRUE, hyper = pc) +
  f(prep_g, model = "rw1", scale.model = TRUE, hyper = pc) +
  f(elev_g, model = "rw1", scale.model = TRUE, hyper = pc) +
  f(nodes, model = euSpde2)

formulas <- list(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11)
formula_names <- c("base","fixed_N","clcL2","bioG","area","year","pH","temp","prec","altitude","space")

###############
## Fit models

Sys.setenv(
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS      = "1",
  OMP_NUM_THREADS      = "3"
)
inla.setOption(num.threads = "3:1")
set.seed(42)

resultsMAR1 <- lapply(1, function(i) {
  message("Running model ", i, " of ", length(formulas))

  tryCatch({
    mod <- inla(
      formulas[[i]],
      family            = "gamma",
      data              = inla.stack.data(stk2),
      control.family    = control.family,
      control.predictor = list(A = inla.stack.A(stk2), compute = FALSE, link = 1),
      control.fixed     = control.fixed,
      control.inla      = control.inla,
      control.compute   = control.compute,
      num.threads       = "3:1"
    )
    modR <- inla.rerun(mod)
    print(summary(modR))
    modR
  }, error = function(e) {
    message("Model ", i, " failed: ", conditionMessage(e))
    NULL
  })
})

###############
## Evidence summary

abs_evidence <- structure(
  list(
    model = c("base","fixed_N","clcL2","bioG","area","year","pH","temp","prec","altitude","space"),
    log_evidence_diff = c(0, 2285.54, 1233.03, 3324.03, 3000.76, -85.84, -43.75, 2599.3, 1519.13, 1519.13, 5748.59)
  ),
  class = "data.frame",
  row.names = c("base","fixed_N","clcL2","bioG","area","year","pH","temp","prec","altitude","space")
)

logZ_reduced <- resultsMAR[[1]]$mlik[1]

mar_evidence <- data.frame(t(sapply(resultsMAR, function(mod) {
  logZ_full <- mod$mlik[1]
  dlogZ <- logZ_full - logZ_reduced
  cbind.data.frame(
    mar_evidence = round(dlogZ, 2),
    dic = mod$dic$dic,
    loo = sum(log(mod$cpo$cpo))
  )
})))

row.names(mar_evidence) <- c("base","fixed","clcL2","bioG","area","year","pH","temp","prec","alt","space")
mar_evidence <- as.data.frame(mar_evidence)

evidence_df <- cbind(abs_evidence, mar_evidence)
evidence_df$dmar <- c(NA, diff(unlist(evidence_df$mar_evidence)))

###############
## Bayesian R2 (median only)

n_samples <- 1000
y_obs <- euSD.dt2$L

bayes_R2_CI <- mclapply(resultsMAR, function(mod) {
  post_samples <- inla.posterior.sample(n_samples, mod)

  idx <- grep("^APredictor", rownames(post_samples[[1]]$latent))

  fitted_samples <- sapply(post_samples, function(s) exp(s$latent[idx]))

  bayes_R2_post <- apply(fitted_samples, 2, function(y_hat) {
    var(y_hat) / (var(y_hat) + var(y_obs - y_hat))
  })

  bayes_R2_q <- quantile(bayes_R2_post, probs = c(0.025, 0.5, 0.975))
  cat("Bayesian R² median:", round(bayes_R2_q[2], 3),
      "\n95% CI:", round(bayes_R2_q[1], 3), "to", round(bayes_R2_q[3], 3), "\n")

  round(bayes_R2_q[2], 3)
}, mc.cores = length(resultsMAR))

evidence_df <- cbind(evidence_df, unlist(bayes_R2_CI))

###############
## Nitrogen coefficients

nitro_coef <- lapply(resultsMAR, function(mod) {
  nitro <- mod$summary.fixed$mean[c(2, 3)]
  names(nitro) <- c("NBud", "NDep")
  nitro
})

evidence_df <- cbind(evidence_df, do.call(rbind, nitro_coef))
