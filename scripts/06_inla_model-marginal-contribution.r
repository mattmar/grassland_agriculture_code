# Sequential contribution of model predictors
# Assumes the objects created in inla_model-part02.R are available:
# euDatDf2, stk2 and euSpde2.

library(parallel)

# -------------------------------------------------------------------------
# Priors
# -------------------------------------------------------------------------

pc <- list(
  prec = list(
    prior = "pc.prec",
    param = c(0.1, 0.01)
  )
)

# -------------------------------------------------------------------------
# INLA controls
# -------------------------------------------------------------------------

control.fixed <- list(
  mean = list(
    intercept = 0,
    NBud = 0,
    NDep = 0,
    default = 0
  ),
  prec = list(
    intercept = 1e-3,
    NBud = 1,
    NDep = 1,
    default = 1
  )
)

control.inla <- list(
  diagonal = 1e-4,
  h = 0.1,
  tolerance = 1e-5,
  strategy = "gaussian",
  improved.simplified.laplace = TRUE
)

control.compute <- list(
  dic = TRUE,
  waic = FALSE,
  cpo = TRUE,
  config = TRUE,
  return.marginals.predictor = FALSE
)

control.family <- list(
  link = "log"
)

# -------------------------------------------------------------------------
# Sequential models
# -------------------------------------------------------------------------

m1 <- L ~ -1 +
  intercept

m2 <- L ~ -1 +
  intercept +
  NBud +
  NDep

m3 <- L ~ -1 +
  intercept +
  NBud +
  NDep +
  f(
    clcL2,
    model = "iid",
    constr = TRUE,
    hyper = list(
      prec = list(
        prior = "pc.prec",
        param = c(1, 0.01)
      )
    ),
    diagonal = 1e-6
  )

m4 <- L ~ -1 +
  intercept +
  NBud +
  NDep +
  f(clcL2, model = "iid", constr = TRUE, hyper = pc) +
  f(bioG, model = "iid", constr = TRUE, hyper = pc)

m5 <- L ~ -1 +
  intercept +
  NBud +
  NDep +
  f(clcL2, model = "iid", constr = TRUE, hyper = pc) +
  f(bioG, model = "iid", constr = TRUE, hyper = pc) +
  f(area_g, model = "rw1", scale.model = TRUE, hyper = pc)

m6 <- L ~ -1 +
  intercept +
  NBud +
  NDep +
  f(clcL2, model = "iid", constr = TRUE, hyper = pc) +
  f(bioG, model = "iid", constr = TRUE, hyper = pc) +
  f(area_g, model = "rw1", scale.model = TRUE, hyper = pc) +
  f(year, model = "rw2", hyper = pc)

m7 <- L ~ -1 +
  intercept +
  NBud +
  NDep +
  f(clcL2, model = "iid", constr = TRUE, hyper = pc) +
  f(bioG, model = "iid", constr = TRUE, hyper = pc) +
  f(area_g, model = "rw1", scale.model = TRUE, hyper = pc) +
  f(year, model = "rw2", hyper = pc) +
  f(
    pH_g,
    model = "rw1",
    scale.model = TRUE,
    hyper = pc,
    diagonal = 1e-6
  )

m8 <- L ~ -1 +
  intercept +
  NBud +
  NDep +
  f(clcL2, model = "iid", constr = TRUE, hyper = pc) +
  f(bioG, model = "iid", constr = TRUE, hyper = pc) +
  f(area_g, model = "rw1", scale.model = TRUE, hyper = pc) +
  f(year, model = "rw2", hyper = pc) +
  f(pH_g, model = "rw1", scale.model = TRUE, hyper = pc, diagonal = 1e-6) +
  f(temp_g, model = "rw1", scale.model = TRUE, hyper = pc)

m9 <- L ~ -1 +
  intercept +
  NBud +
  NDep +
  f(clcL2, model = "iid", constr = TRUE, hyper = pc) +
  f(bioG, model = "iid", constr = TRUE, hyper = pc) +
  f(area_g, model = "rw1", scale.model = TRUE, hyper = pc) +
  f(year, model = "rw2", hyper = pc) +
  f(pH_g, model = "rw1", scale.model = TRUE, hyper = pc, diagonal = 1e-6) +
  f(temp_g, model = "rw1", scale.model = TRUE, hyper = pc) +
  f(prep_g, model = "rw1", scale.model = TRUE, hyper = pc)

m10 <- L ~ -1 +
  intercept +
  NBud +
  NDep +
  f(clcL2, model = "iid", constr = TRUE, hyper = pc) +
  f(bioG, model = "iid", constr = TRUE, hyper = pc) +
  f(area_g, model = "rw1", scale.model = TRUE, hyper = pc) +
  f(year, model = "rw2", hyper = pc) +
  f(pH_g, model = "rw1", scale.model = TRUE, hyper = pc, diagonal = 1e-6) +
  f(temp_g, model = "rw1", scale.model = TRUE, hyper = pc) +
  f(prep_g, model = "rw1", scale.model = TRUE, hyper = pc) +
  f(elev_g, model = "rw1", scale.model = TRUE, hyper = pc)

m11 <- L ~ -1 +
  intercept +
  NBud +
  NDep +
  f(clcL2, model = "iid", constr = TRUE, hyper = pc) +
  f(bioG, model = "iid", constr = TRUE, hyper = pc) +
  f(area_g, model = "rw1", scale.model = TRUE, hyper = pc) +
  f(year, model = "rw2", hyper = pc) +
  f(pH_g, model = "rw1", scale.model = TRUE, hyper = pc, diagonal = 1e-6) +
  f(temp_g, model = "rw1", scale.model = TRUE, hyper = pc) +
  f(prep_g, model = "rw1", scale.model = TRUE, hyper = pc) +
  f(elev_g, model = "rw1", scale.model = TRUE, hyper = pc) +
  f(nodes, model = euSpde2)

formulasMAR <- list(
  base = m1,
  fixed_N = m2,
  clcL2 = m3,
  bioG = m4,
  area = m5,
  year = m6,
  pH = m7,
  temp = m8,
  prec = m9,
  altitude = m10,
  space = m11
)

# -------------------------------------------------------------------------
# Fit models
# -------------------------------------------------------------------------

Sys.setenv(
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OMP_NUM_THREADS = "3"
)

inla.setOption(num.threads = "3:1")
set.seed(42)

resultsMAR <- lapply(
  names(formulasMAR),
  function(modelName) {

    message("Running model: ", modelName)

    tryCatch({

      mod <- inla(
        formulasMAR[[modelName]],
        family = "gamma",
        data = inla.stack.data(stk2),
        control.family = control.family,
        control.predictor = list(
          A = inla.stack.A(stk2),
          compute = FALSE,
          link = 1
        ),
        control.fixed = control.fixed,
        control.inla = control.inla,
        control.compute = control.compute,
        num.threads = "3:1"
      )

      mod <- inla.rerun(mod)

      message("Completed: ", modelName)

      mod

    }, error = function(e) {

      message(
        "Model failed: ",
        modelName,
        " - ",
        conditionMessage(e)
      )

      NULL
    })
  }
)

names(resultsMAR) <- names(formulasMAR)

# -------------------------------------------------------------------------
# Model comparison
# -------------------------------------------------------------------------

logZ <- sapply(
  resultsMAR,
  function(mod) {

    if (is.null(mod)) {
      return(NA_real_)
    }

    as.numeric(mod$mlik[1, 1])
  }
)

dic <- sapply(
  resultsMAR,
  function(mod) {

    if (is.null(mod)) {
      return(NA_real_)
    }

    mod$dic$dic
  }
)

logCPO <- sapply(
  resultsMAR,
  function(mod) {

    if (is.null(mod)) {
      return(NA_real_)
    }

    cpo <- mod$cpo$cpo
    cpo[!is.finite(cpo) | cpo <= 0] <- .Machine$double.eps

    sum(log(cpo))
  }
)

evidenceMAR <- data.frame(
  model = names(resultsMAR),
  log_evidence = logZ,
  delta_log_evidence = c(NA, diff(logZ)),
  DIC = dic,
  delta_DIC = c(NA, diff(dic)),
  log_CPO = logCPO,
  delta_log_CPO = c(NA, diff(logCPO)),
  row.names = NULL
)

# -------------------------------------------------------------------------
# Bayesian R2
# -------------------------------------------------------------------------

nSamples <- 1000
yObs <- euDatDf2$L

bayesR2 <- mclapply(
  resultsMAR,
  function(mod) {

    if (is.null(mod)) {
      return(c(
        R2_low = NA_real_,
        R2_median = NA_real_,
        R2_high = NA_real_
      ))
    }

    postSamples <- inla.posterior.sample(
      nSamples,
      mod
    )

    idx <- grep(
      "^APredictor",
      rownames(postSamples[[1]]$latent)
    )

    fittedSamples <- sapply(
      postSamples,
      function(s) {
        exp(s$latent[idx])
      }
    )

    bayesR2Post <- apply(
      fittedSamples,
      2,
      function(yHat) {
        var(yHat) /
          (var(yHat) + var(yObs - yHat))
      }
    )

    q <- quantile(
      bayesR2Post,
      probs = c(0.025, 0.5, 0.975),
      na.rm = TRUE
    )

    c(
      R2_low = unname(q[1]),
      R2_median = unname(q[2]),
      R2_high = unname(q[3])
    )
  },
  mc.cores = 3
)

bayesR2Mat <- do.call(
  rbind,
  bayesR2
)

evidenceMAR$R2_low <- bayesR2Mat[, "R2_low"]
evidenceMAR$R2_median <- bayesR2Mat[, "R2_median"]
evidenceMAR$R2_high <- bayesR2Mat[, "R2_high"]

# -------------------------------------------------------------------------
# Nitrogen coefficients
# -------------------------------------------------------------------------

nitroCoef <- lapply(
  resultsMAR,
  function(mod) {

    if (is.null(mod)) {
      return(c(
        NBud = NA_real_,
        NDep = NA_real_
      ))
    }

    fixedNames <- rownames(mod$summary.fixed)

    c(
      NBud = if ("NBud" %in% fixedNames) {
        mod$summary.fixed["NBud", "mean"]
      } else {
        NA_real_
      },
      NDep = if ("NDep" %in% fixedNames) {
        mod$summary.fixed["NDep", "mean"]
      } else {
        NA_real_
      }
    )
  }
)

nitroCoef <- do.call(
  rbind,
  nitroCoef
)

evidenceMAR$NBud <- nitroCoef[, "NBud"]
evidenceMAR$NDep <- nitroCoef[, "NDep"]

print(evidenceMAR)
