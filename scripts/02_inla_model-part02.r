# Part 2: Fit the main spatial INLA model
# Run from the repository root. Part 1 must have produced
# ./data/evaVS_formodel.RDS. This script does not depend on
# objects remaining in the Part 1 R session.

library(terra)
library(sf)
library(INLA)
library(inlabru)
library(sp)
library(rmapshaper)

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
# Load model data and boundary
# -------------------------------------------------------------------------

euDatSp2 <- readRDS(
  file.path(dataDir, "evaVS_formodel.RDS")
)

euDatDf2 <- as.data.frame(euDatSp2)

europeSf <- st_read(
  file.path(dataDir, "nuts0_simple_dissolved_cleaned.gpkg"),
  quiet = TRUE
)

europeSimp <- ms_simplify(
  st_as_sf(europeSf),
  keep = 0.1
)

europeSp2 <- as(
  st_transform(europeSimp, crs = 3035),
  "Spatial"
)

euDatDf2$x <- geom(euDatSp2)[,1]
euDatDf2$y <- geom(euDatSp2)[,2]

euDatSp2 <- as(
  st_as_sf(
    euDatDf2,
    coords = c("x", "y"),
    crs = 3035,
    remove = FALSE
  ),
  "Spatial"
)

# -------------------------------------------------------------------------
# Build spatial mesh and barrier SPDE
# -------------------------------------------------------------------------

bdrySeg2 <- fmesher::fm_as_segm(europeSp2)

nAll <- nrow(euDatSp2@coords)
nSub <- min(5000L, nAll)
idxMesh2 <- sample.int(nAll, nSub)

locsMesh2 <- euDatSp2[idxMesh2, ]
locsAll2 <- euDatSp2

euMesh2 <- fmesher::fm_mesh_2d_inla(
  loc = locsMesh2,
  boundary = bdrySeg2,
  max.edge = c(25e3, 100e3),
  offset = c(10e3, 200e3),
  cutoff = 10e3,
  crs = sp::CRS(SRS_string = "EPSG:3035")
)

landTri2 <- fmesher::fm_contains(
  europeSp2,
  euMesh2,
  type = "centroid"
)

allTri2 <- seq_len(nrow(euMesh2$graph$tv))
waterTri2 <- setdiff(allTri2, landTri2)

euSpde2 <- inla.barrier.pcmatern(
  mesh = euMesh2,
  barrier.triangles = waterTri2,
  prior.range = c(50e3, 0.1),
  prior.sigma = c(1.0, 0.1)
)

A2 <- inla.spde.make.A(
  euMesh2,
  loc = locsAll2
)

# -------------------------------------------------------------------------
# Model data
# -------------------------------------------------------------------------

clcIdx2 <- as.integer(factor(euDatDf2$clcL2100))
bioGIdx2 <- as.integer(factor(euDatDf2$bioGregions))

mDat2 <- data.frame(
  intercept = 1,
  NBud = scale(euDatDf2$NFIELD),
  NDep = scale(euDatDf2$NDEP),
  year = as.numeric(as.character(euDatDf2$year)),
  clcL2 = clcIdx2,
  bioG = bioGIdx2,
  area_g = inla.group(
    scale(as.integer(euDatDf2$area), scale = FALSE),
    n = 20,
    method = "cut"
  ),
  pH_g = inla.group(
    euDatDf2$pH,
    n = 25,
    method = "quantile"
  ),
  prep_g = inla.group(
    euDatDf2$prec,
    n = 25,
    method = "cut"
  ),
  temp_g = inla.group(
    euDatDf2$temp,
    n = 50,
    method = "cut"
  ),
  elev_g = inla.group(
    euDatDf2$elev,
    n = 25,
    method = "quantile"
  )
)

stk2 <- inla.stack(
  data = list(L = euDatDf2$L),
  A = list(A2, 1),
  effects = list(
    nodes = seq_len(euSpde2$f$n),
    mDat2
  )
)

# -------------------------------------------------------------------------
# Main model
# -------------------------------------------------------------------------

pc <- list(
  prec = list(
    prior = "pc.prec",
    param = c(0.1, 0.01)
  )
)

formMain <- L ~ -1 +
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

ctrlFixed2 <- list(
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

ctrlInla2 <- list(
  diagonal = 1e-4,
  h = 0.1,
  tolerance = 1e-5,
  strategy = "gaussian",
  int.strategy = "eb",
  improved.simplified.laplace = TRUE
)

ctrlCompute2 <- list(
  dic = TRUE,
  waic = FALSE,
  cpo = FALSE,
  config = TRUE,
  return.marginals.predictor = FALSE
)

ctrlFamily2 <- list(link = "log")

system.time({
  fitMain <- inla(
    formMain,
    family = "gamma",
    data = inla.stack.data(stk2),
    control.family = ctrlFamily2,
    control.predictor = list(
      A = inla.stack.A(stk2),
      compute = FALSE,
      link = 1
    ),
    control.fixed = ctrlFixed2,
    control.inla = ctrlInla2,
    control.compute = ctrlCompute2
  )

  fitMain <- inla.rerun(fitMain)
})

print(summary(fitMain))

# -------------------------------------------------------------------------
# Back-transform nitrogen effects to original units
# -------------------------------------------------------------------------

betaNBud <- fitMain$summary.fixed["NBud", "mean"]
sdNBud <- sd(euDatDf2$NFIELD, na.rm = TRUE)

pctNBud1 <- (
  exp(betaNBud * (1 / sdNBud)) - 1
) * 100

message(
  "% change in L for +1 unit of fertiliser N input: ",
  round(pctNBud1, 3)
)

betaNDep <- fitMain$summary.fixed["NDep", "mean"]
sdNDep <- sd(euDatDf2$NDEP, na.rm = TRUE)

pctNDep1 <- (
  exp(betaNDep * (1 / sdNDep)) - 1
) * 100

message(
  "% change in L for +1 unit of atmospheric N deposition: ",
  round(pctNDep1, 3)
)

# -------------------------------------------------------------------------
# Save outputs used by figure scripts
# -------------------------------------------------------------------------

saveRDS(
  fitMain,
  file.path(dataDir, "fitMain.RDS")
)

write.csv(
  mDat2,
  file.path(dataDir, "mDat2.csv"),
  row.names = FALSE
)
