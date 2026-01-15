### Main INLA model (paper results)
## Assumes Part 1 has been run and these objects exist:
## baseDir, euDatSp, badIdx, europeSimp (sf), euDatDf (optional)
set.seed(42)

# Free some memory
rm(fitNull)
gc()

## Parallel settings
Sys.setenv(
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS      = "1",
  OMP_NUM_THREADS      = "3"
)
inla.setOption(num.threads = "3:1")

## Data subset (exclude observations flagged in Part 1)
euDatSp2 <- euDatSp[-badIdx, ]

## Model data.frame
euDatDf2 <- as.data.frame(euDatSp2)

## Convert boundary + points to sp objects for INLA
europeSp2 <- as(st_transform(europeSimp, crs = 3035), "Spatial")
euDatSp2  <- as(st_as_sf(euDatDf2, coords = c("x", "y"), crs = 3035, remove = FALSE), "Spatial")

## Mesh boundary segments
bdrySeg2 <- fmesher::fm_as_segm(europeSp2)

## Thin points for mesh construction (all points kept for projector)
nAll  <- nrow(euDatSp2@coords)
nSub  <- min(5000L, nAll)
idxMesh2 <- sample.int(nAll, nSub)

locsMesh2 <- euDatSp2[idxMesh2, ]
locsAll2  <- euDatSp2

## Mesh (parameters unchanged)
euMesh2 <- fmesher::fm_mesh_2d_inla(
  loc      = locsMesh2,
  boundary = bdrySeg2,
  max.edge = c(25e3, 100e3),
  offset   = c(10e3, 200e3),
  cutoff   = 10e3,
  crs      = sp::CRS(SRS_string = "EPSG:3035")
)

## Barrier triangles
landTri2  <- fmesher::fm_contains(europeSp2, euMesh2, type = "centroid")
allTri2   <- seq_len(nrow(euMesh2$graph$tv))
waterTri2 <- setdiff(allTri2, landTri2)

euSpde2 <- inla.barrier.pcmatern(
  mesh              = euMesh2,
  barrier.triangles = waterTri2,
  prior.range       = c(50e3, 0.1),
  prior.sigma       = c(1.0, 0.1)
)

## Projector matrix
A2 <- inla.spde.make.A(euMesh2, loc = locsAll2)

## Indices for IID terms
clcIdx2  <- as.integer(factor(euDatDf2$clcL2100))
bioGIdx2 <- as.integer(factor(euDatDf2$bioGregions))

## Grouped smooth covariates
mDat2 <- data.frame(
  intercept = 1,
  avgT      = scale(euDatDf2$temp),
  cumPrep   = scale(euDatDf2$prec),
  altitude  = scale(euDatDf2$elev),
  pH        = scale(euDatDf2$pH, scale = FALSE),
  NBud      = scale(euDatDf2$NFIELD),
  NDep      = scale(euDatDf2$NDEP),
  year      = as.numeric(as.character(euDatDf2$year)),
  clcL2     = clcIdx2,
  bioG      = bioGIdx2,
  area_g    = inla.group(scale(as.integer(euDatDf2$area), scale = FALSE), n = 20, method = "cut"),
  pH_g      = inla.group(euDatDf2$pH,   n = 25, method = "quantile"),
  prep_g    = inla.group(euDatDf2$prec, n = 25, method = "cut"),
  temp_g    = inla.group(euDatDf2$temp, n = 50, method = "cut"),
  elev_g    = inla.group(euDatDf2$elev, n = 25, method = "quantile")
)

## INLA stack
stk2 <- inla.stack(
  data    = list(L = euDatDf2$L),
  A       = list(A2, 1),
  effects = list(nodes = 1:euSpde2$f$n, mDat2)
)

## Formula 
pc <- list(prec = list(prior = "pc.prec", param = c(0.1, 0.01)))

formMain <- L ~ 1 + NBud + NDep +
  f(clcL2, model = "iid", constr = TRUE, hyper = pc) +
  f(bioG,  model = "iid", constr = TRUE, hyper = pc) +
  f(area_g, model = "rw1", scale.model = TRUE, hyper = pc) +
  f(year,   model = "rw2", hyper = pc) +
  f(pH_g,   model = "rw1", scale.model = TRUE, hyper = pc, diagonal = 1e-6) +
  f(temp_g, model = "rw1", scale.model = TRUE, hyper = pc) +
  f(prep_g, model = "rw1", scale.model = TRUE, hyper = pc) +
  f(elev_g, model = "rw1", scale.model = TRUE, hyper = pc) +
  f(nodes, model = euSpde2)

## Controls
ctrlFixed2 <- list(mean.intercept = 0, prec.intercept = 1e-3)

ctrlInla2 <- list(
  diagonal = 1e-4,
  h        = 0.1,
  tolerance= 1e-5,
  strategy = "simplified.laplace",
  improved.simplified.laplace = TRUE
)

ctrlCompute2 <- list(
  dic    = TRUE,
  waic   = FALSE,
  cpo    = FALSE,
  config = TRUE,
  return.marginals.predictor = FALSE
)

ctrlFamily2 <- list(link = "log")

## Fit
system.time({
  fitMain <- inla(
    formMain,
    family            = "gamma",
    data              = inla.stack.data(stk2),
    control.family    = ctrlFamily2,
    control.predictor = list(A = inla.stack.A(stk2), compute = FALSE, link = 1),
    control.fixed     = ctrlFixed2,
    control.inla      = ctrlInla2,
    control.compute   = ctrlCompute2
  )
  fitMain <- inla.rerun(fitMain)
})

print(summary(fitMain))

## Back-transform N effects to original units (% change in L)
betaNBud <- fitMain$summary.fixed["NBud", "mean"]
sdNBud   <- sd(euDatDf2$NFIELD, na.rm = TRUE)

## % change in L for +10 units of NBud (original scale)
pctNBud10 <- (exp(betaNBud * (1 / sdNBud)) - 1) * 100
message("% change in L for +1 units of NBud: ", round(pctNBud10, 3))

betaNDep <- fitMain$summary.fixed["NDep", "mean"]
sdNDep   <- sd(euDatDf2$NDEP, na.rm = TRUE)

## % change in L for +1 unit of NDep (original scale)
pctNDep1 <- (exp(betaNDep * (1 / sdNDep)) - 1) * 100
message("% change in L for +1 unit of NDep: ", round(pctNDep1, 3))