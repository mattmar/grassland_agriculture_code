# -------------------------------------------------------------------------
# Absolute contribution of individual model predictors
# -------------------------------------------------------------------------
# Each predictor is fitted separately against an intercept-only baseline.
# Contribution is expressed as the difference in log marginal likelihood
# relative to the baseline model.

pc <- list(
  prec = list(
    prior = "pc.prec",
    param = c(0.1, 0.01)
  )
)

# -------------------------------------------------------------------------
# INLA controls
# -------------------------------------------------------------------------

ctrlFixedAbs <- list(
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

ctrlInlaAbs <- list(
  diagonal = 1e-4,
  h = 0.1,
  tolerance = 1e-5,
  strategy = "gaussian",
  int.strategy = "eb",
  improved.simplified.laplace = TRUE
)

ctrlComputeAbs <- list(
  dic = TRUE,
  waic = FALSE,
  cpo = FALSE,
  config = TRUE,
  return.marginals.predictor = FALSE
)

ctrlFamilyAbs <- list(
  link = "log"
)

# -------------------------------------------------------------------------
# Individual predictor models
# -------------------------------------------------------------------------

fBase <- L ~ -1 + intercept

fNBud <- L ~ -1 +
  intercept +
  NBud

fNDep <- L ~ -1 +
  intercept +
  NDep

fClc <- L ~ -1 +
  intercept +
  f(
    clcL2,
    model = "iid",
    constr = TRUE,
    hyper = pc
  )

fBioG <- L ~ -1 +
  intercept +
  f(
  bioG,
  model = "iid",
  constr = TRUE,
  hyper = list(
    prec = list(
      prior = "pc.prec",
      param = c(2, 0.01)
    )
  ),
  diagonal = 1e-6
)

fArea <- L ~ -1 +
  intercept +
  f(
    area_g,
    model = "rw1",
    scale.model = TRUE,
    hyper = pc
  )

fYear <- L ~ -1 +
  intercept +
  f(
    year,
    model = "rw2",
    hyper = pc
  )

fPH <- L ~ -1 +
  intercept +
  f(
    pH_g,
    model = "rw1",
    scale.model = TRUE,
    hyper = pc,
    diagonal = 1e-6
  )

fTemp <- L ~ -1 +
  intercept +
  f(
    temp_g,
    model = "rw1",
    scale.model = TRUE,
    hyper = pc
  )

fPrep <- L ~ -1 +
  intercept +
  f(
    prep_g,
    model = "rw1",
    scale.model = TRUE,
    hyper = pc
  )

fElev <- L ~ -1 +
  intercept +
  f(
    elev_g,
    model = "rw1",
    scale.model = TRUE,
    hyper = pc
  )

fSpace <- L ~ -1 +
  intercept +
  f(
    nodes,
    model = euSpde2
  )

formulasAbs <- list(
  base = fBase,
  N_budget = fNBud,
  N_deposition = fNDep,
  land_cover = fClc,
  biogeographic_region = fBioG,
  plot_area = fArea,
  year = fYear,
  pH = fPH,
  temperature = fTemp,
  precipitation = fPrep,
  elevation = fElev,
  space = fSpace
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

resultsAbs <- lapply(
  names(formulasAbs),
  function(modelName) {

    cat(
      "\nRunning model:",
      modelName,
      "\n"
    )

    tryCatch({

      mod <- inla(
        formulasAbs[[modelName]],
        family = "gamma",
        data = inla.stack.data(stk2),
        control.family = ctrlFamilyAbs,
        control.predictor = list(
          A = inla.stack.A(stk2),
          compute = FALSE,
          link = 1
        ),
        control.fixed = ctrlFixedAbs,
        control.inla = ctrlInlaAbs,
        control.compute = ctrlComputeAbs
      )

      mod <- inla.rerun(mod)

      cat(
        "Completed:",
        modelName,
        "\n"
      )

      mod

    }, error = function(e) {

      cat(
        "Model failed:",
        modelName,
        "\n",
        conditionMessage(e),
        "\n"
      )

      NULL
    })
  }
)

names(resultsAbs) <- names(formulasAbs)

# -------------------------------------------------------------------------
# Difference in log marginal likelihood from intercept-only model
# -------------------------------------------------------------------------

logZBase <- as.numeric(
  resultsAbs[["base"]]$mlik[1, 1]
)

logEvidenceDiff <- sapply(
  resultsAbs,
  function(mod) {

    if (is.null(mod)) {
      return(NA_real_)
    }

    as.numeric(mod$mlik[1, 1]) - logZBase
  }
)

absEvidence <- data.frame(
  predictor = names(logEvidenceDiff),
  log_evidence_diff = round(logEvidenceDiff, 2),
  row.names = NULL
)

absEvidence <- absEvidence[
  order(
    absEvidence$log_evidence_diff,
    decreasing = TRUE,
    na.last = TRUE
  ),
]

print(absEvidence)