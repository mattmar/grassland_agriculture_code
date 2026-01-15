### Visualise random effects from INLA model
## Assumes these objects already exist in memory:
## fitMain (or h.inlaf2), euDatDf2, mDat2

## Libraries
library(INLA)
library(ggplot2)
library(ggpubr)

## Model object
bestM <- fitMain

## Helper: transform marginal to response scale (relative change)
toResp <- function(m) INLA::inla.tmarginal(function(x) exp(x) - 1, m)

## Helper: posterior mean + 95% CI from a list of marginals
summMarg <- function(mList) {
  mResp <- lapply(mList, toResp)
  data.frame(
    mean = sapply(mResp, INLA::inla.emarginal, fun = identity),
    lwr  = sapply(mResp, INLA::inla.qmarginal, p = 0.025),
    upr  = sapply(mResp, INLA::inla.qmarginal, p = 0.975)
  )
}

## Helper: standard RW plot
plotRw <- function(df, xVar, xLab, yLab = NULL) {
  ggplot(df, aes(x = .data[[xVar]], y = mean)) +
    geom_point(colour = "grey50", alpha = 0.6, size = 2) +
    geom_smooth(method = "loess", span = 0.9, se = TRUE, colour = "grey30", alpha = 0.25) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
    labs(x = xLab, y = yLab) +
    theme_minimal()
}

####################
## Biogeographical region (iid)
bioGSum <- bestM$summary.random$bioG
bioGSum[, 2:7] <- exp(bioGSum[, 2:7]) - 1

## Map iid index back to original bioGregions codes
bioGFac <- factor(euDatDf2$bioGregions)
bioGLook <- data.frame(
  bioG = seq_along(levels(bioGFac)),
  bioGregions = as.integer(levels(bioGFac))
)

bioGSum <- merge(bioGSum, bioGLook, by.x = "ID", by.y = "bioG", all.x = TRUE)

bioGLabels <- c(
  `1`  = "Alpine",
  `4`  = "Atlantic",
  `5`  = "Black Sea",
  `6`  = "Boreal",
  `7`  = "Continental",
  `9`  = "Mediterranean",
  `11` = "Steppic",
  `12` = "Transition/Other"
)

bioGSum$region <- bioGLabels[as.character(bioGSum$bioGregions)]

## Sample sizes per region (uses the model index, not the original codes)
bioGCounts <- aggregate(n ~ bioG, data = data.frame(bioG = mDat2$bioG, n = 1), sum)
bioGSum <- merge(bioGSum, bioGCounts, by.x="bioGregions", by.y = "bioG", all.x = TRUE)

plotBioG <- data.frame(
  region  = factor(bioGSum$region),
  mean    = bioGSum$`0.5quant`,
  lwr     = bioGSum$`0.025quant`,
  upr     = bioGSum$`0.975quant`,
  samples = bioGSum$n
)

plotBioG <- plotBioG[!is.na(plotBioG$region), ]
plotBioG$region <- factor(plotBioG$region, levels = levels(plotBioG$region)[order(plotBioG$mean)], ordered = TRUE)

ggBioG <- ggplot(plotBioG, aes(x = region, y = mean)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.35, colour = "grey50") +
  geom_point(size = 3, shape = 22, fill = "grey85") +
  geom_text(aes(label = paste0("n=", samples)), y = 0.9, size = 3) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
  labs(x = "Biogeographical region", y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

####################
## Year (rw2)
yearSum <- summMarg(bestM$marginals.random$year)
plotYear <- data.frame(
  year = 1990:2021,
  yearSum
)

ggYear <- plotRw(plotYear, "year", "Year")

####################
## Area (rw1; grouped)
areaSum <- summMarg(bestM$marginals.random$area_g)
areaVals <- tapply(euDatDf2$area, mDat2$area_g, mean)

plotArea <- data.frame(
  area = as.numeric(areaVals),
  areaSum
)

ggArea <- plotRw(plotArea, "area", expression(Survey~Area~(m^2)))

####################
## Temperature (rw1; grouped)
tempSum <- summMarg(bestM$marginals.random$temp_g)
tempVals <- tapply(euDatDf2$temp, mDat2$temp_g, mean)

plotTemp <- data.frame(
  temp = as.numeric(tempVals),
  tempSum
)

ggTemp <- plotRw(plotTemp, "temp", expression(Average~Temperature~("\u00B0"*C)), "Relative change in L")

####################
## Precipitation (rw1; grouped)
prepSum <- summMarg(bestM$marginals.random$prep_g)
prepVals <- tapply(euDatDf2$prec, mDat2$prep_g, mean)

plotPrep <- data.frame(
  prec = as.numeric(prepVals),
  prepSum
)

ggPrep <- plotRw(plotPrep, "prec", expression(Cumulative~Precipitation~(mm%.%year^{-1})))

####################
## Elevation (rw1; grouped)
elevSum <- summMarg(bestM$marginals.random$elev_g)
elevVals <- tapply(euDatDf2$elev, mDat2$elev_g, mean)

plotElev <- data.frame(
  elev = as.numeric(elevVals),
  elevSum
)

ggElev <- plotRw(plotElev, "elev", "Elevation (m)")

####################
## pH (rw1; grouped)
pHSum <- summMarg(bestM$marginals.random$pH_g)
pHVals <- tapply(euDatDf2$pH, mDat2$pH_g, mean)

plotPH <- data.frame(
  pH = as.numeric(pHVals),
  pHSum
)

ggPH <- plotRw(plotPH, "pH", "pH")

####################
## Arrange + save
gArr <- ggpubr::ggarrange(
  ggTemp, ggPrep,
  ggElev, ggPH,
  ggArea, ggYear,
  NULL,   ggBioG, NULL,
  ncol = 3, nrow = 3,
  labels = c("A","B","C","D","E","F","","G","")
)

plot(gArr)