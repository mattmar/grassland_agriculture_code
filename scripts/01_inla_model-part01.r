### Prepare dataset for spatial SPDE model (INLA)

## Libraries
library(terra)
library(sf)
library(INLA)
library(inlabru)
library(parallel)
library(sp)
library(ggplot2)
library(rmapshaper)
library(raster)

## Parallel settings
Sys.setenv(
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS      = "1",
  OMP_NUM_THREADS      = "1"
)
inla.setOption(num.threads = "1:1")

## Paths
baseDir <- "~/git/grassland_agriculture_code/"
setwd(baseDir)
source(file.path(baseDir, "data", "NNcoords.R"))

## Boundaries
europeSf <- st_read(file.path(baseDir, "data", "nuts0_simple_dissolved_cleaned.gpkg"))

## Data with response + covariates
euDatSp <- readRDS(file.path(baseDir, "./data/evaVS_covariates.RDS"))

## Year handling
euDatSp$year <- factor(euDatSp$year, ordered = TRUE)

## Recode: forest surveys from 2022–2023 moved to 2021
euDatSp[euDatSp$year %in% c(2022, 2023) & euDatSp$clcL2500 == "Forest", ]$year <- 2021

## Harmonise L2 class label
euDatSp$clcL2100[euDatSp$clcL2100 %in% "artificial_veg"] <- "het_agri_areas"

## Cap years at 2021 (later years have low sample size)
euDatSp$year[as.numeric(as.character(euDatSp$year)) > 2021] <- 2021
euDatSp$year <- droplevels(euDatSp$year)

## Data exclusions / transforms used throughout
euDatSp <- euDatSp[-which(euDatSp$uncer > 100), ]    # exclude high spatial-uncertainty observations
euDatSp$pH <- euDatSp$pH / 1000                      # pH scaling (input stored x1000)
euDatSp$NFIELD <- ifelse(euDatSp$NFIELD >= 400, 400, euDatSp$NFIELD)  # cap N budget

## Coerce for INLA mesh + SPDE
euDatDf <- as.data.frame(euDatSp)
euDatPt <- SpatialPoints(geom(euDatSp)[, c("x", "y")], proj4string = CRS("EPSG:3035"))

## Model boundary (simplified)
europeSimp <- ms_simplify(st_as_sf(europeSf), keep = 0.1)

## Ensure consistent projection (LAEA, EPSG:3035)
europeV   <- vect(europeSimp)
euDatV    <- euDatSp

europeLaea <- project(europeV, "EPSG:3035")
euDatLaea  <- project(euDatV, "EPSG:3035")

europeSp <- as(europeLaea, "Spatial")
euDatSp  <- as(euDatLaea, "Spatial")

bdrySeg <- inla.sp2segment(as(europeSimp, "Spatial"))

## Thin points for mesh construction
set.seed(42)
idxMesh   <- sample(seq_len(nrow(euDatSp)), 5000)
locsMesh  <- crop(euDatSp[idxMesh, ], europeSimp)
locsAll   <- euDatSp

## Mesh parameters
autoRange <- diff(range(coordinates(rbind(euDatPt)))) / 10

edgeInner <- autoRange / 4
edgeOuter <- edgeInner * 4
cutoff    <- edgeInner / 4

euMesh <- inla.mesh.2d(
  loc      = locsMesh,
  boundary = bdrySeg,
  max.edge = c(edgeInner, edgeOuter),
  offset   = c(edgeInner, edgeOuter),
  cutoff   = cutoff,
  crs      = CRS("+init=EPSG:3035")
)

## Barrier-SPDE
landTri <- inla.over_sp_mesh(
  europeSp,
  euMesh,
  type = "centroid",
  ignore.CRS = FALSE
)

euSpde <- inla.barrier.pcmatern(
  mesh              = euMesh,
  barrier.triangles = landTri,
  prior.range       = c(autoRange, 0.5),
  prior.sigma       = c(0.1, 0.01)
)

## Projector matrix
A <- inla.spde.make.A(euMesh, loc = locsAll)

## Model frame
mDat <- data.frame(
  intercept = 1,
  avgT      = scale(euDatDf$temp),
  cumPrep   = scale(euDatDf$prec),
  altitude  = scale(euDatDf$elev),
  pH        = scale(euDatDf$pH, scale = FALSE),
  NBud      = scale(euDatDf$NFIELD),
  NDep      = scale(euDatDf$NDEP),
  year      = euDatDf$year,
  clcL2     = as.integer(as.factor(euDatDf$clcL2100)),
  bioG      = euDatDf$bioGregions,
  area      = scale(as.integer(euDatDf$area), scale = FALSE)
)

stkNull <- inla.stack(
  data = list(L = euDatDf$L),
  A = list(A, 1),
  effects = list(nodes = 1:euSpde$f$n, mDat)
)

## Null model: climate + edaphic + land cover + spatial field
formNull <- L ~ 1 +
  f(clcL2, model = "iid", constr = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
  f(bioG, model = "iid", constr = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
  f(area, model = "rw2", scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
  f(year, model = "rw1",
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
  f(inla.group(pH, n = 10, "quantile"), model = "rw1", scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
  f(inla.group(cumPrep, n = 25), model = "rw1", scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
  f(inla.group(avgT, n = 50), model = "rw1", scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
  f(inla.group(altitude, n = 25), model = "rw1", scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
  f(nodes, model = euSpde)

ctrlFixed <- list(
  mean.intercept = 0, prec.intercept = 1e-3,
  mean = 0, prec = 1
)

ctrlInla <- list(
  strategy = "simplified.laplace",
  improved.simplified.laplace = TRUE,
  diagonal = 1e-5
)

ctrlCompute <- list(
  dic    = TRUE,
  waic   = TRUE,
  cpo    = TRUE,
  config = TRUE,
  return.marginals.predictor = FALSE
)

ctrlFamily <- list(
  link  = "log",
  hyper = list(prec = list(prior = "pc.prec", param = c(0.5, 0.5)))
)

## Takes time
system.time({
  fitNull <- inla(
    formNull,
    family            = "gamma",
    data              = inla.stack.data(stkNull),
    control.family    = ctrlFamily,
    control.predictor = list(A = inla.stack.A(stkNull), compute = FALSE, link = 1),
    control.fixed     = ctrlFixed,
    control.inla      = ctrlInla,
    control.compute   = ctrlCompute
  )
})

## Diagnostics used to flag observations inconsistent with the null model
cpo  <- fitNull$cpo$cpo
pit  <- fitNull$cpo$pit
fail <- fitNull$cpo$failure  # 1 = CPO failed

## Replace failed/invalid CPOs so filtering steps are stable
cpo[fail == 1 | !is.finite(cpo)] <- .Machine$double.eps

## Predictive extremes (PIT tails)
pitHigh <- which(pit >= 0.95)
pitLow  <- which(pit <= 0.05)

## Poor fit under the null model
cpoPoor <- which(cpo <= 0.05)

## Nitrogen budget quantiles (based on scaled NBud used in model frame)
nBudHigh <- which(mDat$NBud > quantile(mDat$NBud, p = 0.90, na.rm = TRUE))
nBudLow  <- which(mDat$NBud < quantile(mDat$NBud, p = 0.60, na.rm = TRUE))

## Points removed: (i) extreme PIT tail, (ii) poor CPO, (iii) high/low N budget
badIdxHigh <- Reduce(intersect, list(pitHigh, cpoPoor, nBudHigh))
badIdxLow  <- Reduce(intersect, list(pitLow,  cpoPoor, nBudLow))

## Additional filter: very high N budget + very high diversity in arable land
badIdxArable <- which(euDatDf$NFIELD > 100 & euDatDf$L > 3 & euDatDf$clcL2100 %in% "arable_land")
badIdx <- unique(c(badIdxLow, badIdxHigh, badIdxArable))

## Percent of points that will be removed
paste(round(length(badIdx) / nrow(mDat) * 100,2),"% of the initial dataset")

# Follow with inla_model-part02.r
