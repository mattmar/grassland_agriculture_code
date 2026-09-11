# Part 1: Prepare the dataset for the spatial INLA model
# Run from the repository root.

library(terra)
library(sf)
library(INLA)
library(inlabru)
library(sp)
library(rmapshaper)
library(raster)

# -------------------------------------------------------------------------
# Settings and paths
# -------------------------------------------------------------------------

set.seed(42)

Sys.setenv(
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OMP_NUM_THREADS = "12"
)
inla.setOption(num.threads = "12:1")

baseDir <- normalizePath("/home/matteo/git/grassland_agriculture_code")
dataDir <- file.path(baseDir, "data")

# -------------------------------------------------------------------------
# Load and preprocess data
# -------------------------------------------------------------------------

europeSf <- st_read(
  file.path(dataDir, "nuts0_simple_dissolved_cleaned.gpkg"),
  quiet = TRUE
)

euDatSp <- readRDS(
  file.path(dataDir, "evaVS_covariates.RDS")
)

# Harmonise year and land-use variables
euDatSp$year <- factor(euDatSp$year, ordered = TRUE)

euDatSp[
  euDatSp$year %in% c(2022, 2023) &
    euDatSp$clcL2500 == "Forest",
]$year <- 2021

euDatSp$clcL2100[
  euDatSp$clcL2100 == "artificial_veg"
] <- "het_agri_areas"

euDatSp$year[
  as.numeric(as.character(euDatSp$year)) > 2021
] <- 2021

euDatSp$year <- droplevels(euDatSp$year)

# Apply data transformations used in the analysis
euDatSp <- euDatSp[-which(euDatSp$uncer > 100), ]
euDatSp$pH <- euDatSp$pH / 1000
euDatSp$NFIELD <- pmin(euDatSp$NFIELD, 400)

euDatDf <- as.data.frame(euDatSp)

# -------------------------------------------------------------------------
# Build spatial mesh and barrier SPDE
# -------------------------------------------------------------------------

euDatPt <- SpatialPoints(
  geom(euDatSp)[, c("x", "y")],
  proj4string = CRS("EPSG:3035")
)

europeSimp <- ms_simplify(
  st_transform(st_as_sf(europeSf), crs = 3035),
  keep = 0.1
)

europeLaea <- project(vect(europeSimp), "EPSG:3035")
euDatLaea <- project(euDatSp, "EPSG:3035")

europeSp <- as(europeLaea, "Spatial")
euDatModelSp <- as(euDatLaea, "Spatial")

bdrySeg <- inla.sp2segment(
  as(europeSimp, "Spatial")
)

# Use a random subset of points for mesh construction
idxMesh <- sample(seq_len(nrow(euDatModelSp)), 5000)

locsMesh <- raster::crop(
  euDatModelSp[idxMesh, ],
  as(europeSimp, "Spatial")
)

locsAll <- euDatModelSp

autoRange <- diff(range(coordinates(euDatPt))) / 40
edgeInner <- autoRange / 4
edgeOuter <- edgeInner * 4
cutoff <- edgeInner / 4

euMesh <- inla.mesh.2d(
  loc = locsMesh,
  boundary = bdrySeg,
  max.edge = c(edgeInner, edgeOuter),
  offset = c(edgeInner, edgeOuter),
  cutoff = cutoff,
  crs = CRS("+init=EPSG:3035")
)

landTri <- inla.over_sp_mesh(
  europeSp,
  euMesh,
  type = "centroid",
  ignore.CRS = FALSE
)

euSpde <- inla.barrier.pcmatern(
  mesh = euMesh,
  barrier.triangles = landTri,
  prior.range = c(autoRange, 0.5),
  prior.sigma = c(0.1, 0.01)
)

A <- inla.spde.make.A(
  euMesh,
  loc = locsAll
)

# -------------------------------------------------------------------------
# Baseline model excluding nitrogen predictors
# -------------------------------------------------------------------------

mDat <- data.frame(
  intercept = 1,
  avgT = scale(euDatDf$temp),
  cumPrep = scale(euDatDf$prec),
  altitude = scale(euDatDf$elev),
  pH = scale(euDatDf$pH, scale = FALSE),
  NBud = scale(euDatDf$NFIELD),
  NDep = scale(euDatDf$NDEP),
  year = euDatDf$year,
  clcL2 = as.integer(as.factor(euDatDf$clcL2100)),
  bioG = euDatDf$bioGregions,
  area = scale(as.integer(euDatDf$area), scale = FALSE)
)

stkNull <- inla.stack(
  data = list(L = euDatDf$L),
  A = list(A, 1),
  effects = list(
    nodes = seq_len(euSpde$f$n),
    mDat
  )
)

formNull <- L ~ -1 +
  f(clcL2, model = "iid", constr = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
  f(bioG, model = "iid", constr = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
  f(area, model = "rw1", scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
  f(year, model = "rw1",
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
  f(inla.group(pH, n = 10, method = "quantile"),
    model = "rw1", scale.model = TRUE, diagonal = 1e-6,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
  f(inla.group(cumPrep, n = 10),
    model = "rw1", scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
  f(inla.group(avgT, n = 10),
    model = "rw1", scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
  f(inla.group(altitude, n = 10),
    model = "rw1", scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
  f(nodes, model = euSpde)

ctrlFixed <- list(
  mean.intercept = 0,
  prec.intercept = 1e-3,
  mean = c(0, 0, 0, 0, 0, 0),
  prec = c(1, 1, 1, 1, 1, 1)
)

ctrlInla <- list(
  strategy = "simplified.laplace",
  improved.simplified.laplace = TRUE,
  diagonal = 1e-4,
  h = 0.1
)

ctrlCompute <- list(
  dic = TRUE,
  waic = TRUE,
  cpo = TRUE,
  config = TRUE,
  return.marginals.predictor = FALSE
)

ctrlFamily <- list(
  link = "log",
  hyper = list(
    prec = list(prior = "pc.prec", param = c(0.5, 0.5))
  )
)

system.time({
  fitNull <- inla(
    formNull,
    family = "gamma",
    data = inla.stack.data(stkNull),
    control.family = ctrlFamily,
    control.predictor = list(
      A = inla.stack.A(stkNull),
      compute = FALSE,
      link = 1
    ),
    control.fixed = ctrlFixed,
    control.inla = ctrlInla,
    control.compute = ctrlCompute
  )
})

# -------------------------------------------------------------------------
# Identify observations excluded from the final model
# -------------------------------------------------------------------------

cpo <- fitNull$cpo$cpo
pit <- fitNull$cpo$pit
fail <- fitNull$cpo$failure

cpo[fail == 1 | !is.finite(cpo)] <- .Machine$double.eps

pitHigh <- which(pit >= 0.95)
pitLow <- which(pit <= 0.05)
cpoPoor <- which(cpo <= 0.05)

nBudHigh <- which(
  mDat$NBud > quantile(mDat$NBud, 0.90, na.rm = TRUE)
)

nBudLow <- which(
  mDat$NBud < quantile(mDat$NBud, 0.90, na.rm = TRUE)
)

badIdxHigh <- Reduce(
  intersect,
  list(pitHigh, cpoPoor, nBudHigh)
)

badIdxLow <- Reduce(
  intersect,
  list(pitLow, cpoPoor, nBudLow)
)

badIdxArable <- which(
  euDatDf$NFIELD > 100 &
    euDatDf$L > 10 &
    euDatDf$clcL2100 == "arable_land"
)

badIdx <- unique(c(
  badIdxLow,
  badIdxHigh,
  badIdxArable
))

message(
  round(length(badIdx) / nrow(mDat) * 100, 2),
  "% of observations removed"
)

# Save the filtered source data used by Part 2
euDatSp2 <- euDatSp[-badIdx, ]

saveRDS(
  euDatSp2,
  file.path(dataDir, "evaVS_formodel.RDS")
)
