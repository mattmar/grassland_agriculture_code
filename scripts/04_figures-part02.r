### Visualise random effects from INLA model
## Assumes these objects already exist in memory:
## fitMain, euDatDf2, mDat2

library(INLA)
library(ggplot2)
library(sf)
library(terra)
library(scales)
library(gridExtra)
library(viridis)
library(patchwork)
library(reshape2)
library(ggdist)
library(ggbreak)
library(ggpubr)
library(dplyr)
library(grid)
library(tidytext)
library(forcats)

baseDir <- normalizePath("../git/grassland_agriculture_code/")
data_dir <- file.path(baseDir, "data")
noteapp_dir  <- file.path(baseDir, "Figures_notes")

fitMain <- readRDS(file.path(data_dir, "fitMain.RDS"))
euDatDf2 <- as.data.frame(readRDS(file.path(data_dir, "evaVS_formodel.RDS")))
mDat2 <- read.csv(file.path(data_dir, "mDat2.csv"))
euSpde2 <- readRDS(file.path(data_dir,"euSpde2.RDS"))
euSD <- as.data.frame(readRDS(file.path(data_dir, "evaVS_formodel.RDS")))

# -------------------------------------------------------------------------
# Appendix Figure 1: Land use share
# -------------------------------------------------------------------------

# Sorrounding land use share for each type of grasslands
euSD$clcL2 <- as.factor(euSD$clcL2100)
euSD <- euSD[-which(is.na(euSD$clcL2100)), ]
euSD$clcL2 <- droplevels(euSD$clcL2, exclude = "<NA>")

levels(euSD$clcL2) <- clcL2Labs

levels3100 <- c(
  "1.1.1 Continuous urban fabric", "1.1.2 Discontinuous urban fabric",
  "1.2.1 Industrial or commercial units", "1.2.2 Road and rail networks and associated land",
  "1.2.3 Port areas", "1.2.4 Airports",
  "1.3.1 Mineral extraction sites", "1.3.2 Dump sites", "1.3.3 Construction sites",
  "1.4.1 Green urban areas", "1.4.2 Sport and leisure facilities",
  "2.1.1 Non-irrigated arable land", "2.1.2 Permanently irrigated land", "2.1.3 Rice fields",
  "2.2.1 Vineyards", "2.2.2 Fruit trees and berry plantations", "2.2.3 Olive groves",
  "2.3.1 Pastures",
  "2.4.1 Annual crops associated with permanent crops", "2.4.2 Complex cultivation patterns",
  "2.4.3 Agri land, with large areas of natural vegetation", "2.4.4 Agro-forestry areas",
  "3.1.1 Broad-leaved forest", "3.1.2 Coniferous forest", "3.1.3 Mixed forest",
  "3.2.1 Natural grasslands", "3.2.2 Moors and heathland", "3.2.3 Sclerophyllous vegetation",
  "3.2.4 Transitional woodland-shrub",
  "3.3.1 Beaches, dunes, sands", "3.3.2 Bare rocks", "3.3.3 Sparsely vegetated areas",
  "3.3.4 Burnt areas", "3.3.5 Glaciers and perpetual snow",
  "4.1.1 Inland marshes", "4.1.2 Peat bogs",
  "4.2.1 Salt marshes", "4.2.2 Salines", "4.2.3 Intertidal flats",
  "5.1.1 Water courses", "5.1.2 Water bodies",
  "5.2.1 Coastal lagoons", "5.2.2 Estuaries", "5.2.3 Sea and ocean",
  "No Data", "No Data", "No Data"
)

levels3100 <- data.frame(
  id       = seq_along(levels3100),
  clcL3100n = gsub("^[0-9.]+ ", "", levels3100)
)

euSD$clcL2 <- factor(euSD$clcL2, levels = clcL2Labs, ordered = TRUE)
euSDm <- merge(euSD, levels3100, by.y = "id", by.x = "clcL3100", all.x = TRUE)

## Count plots in each original CLC class
plot_counts <- euSDm %>%
  count(clcL2, clcL3100n, name = "n")

## Order aggregated groups
group_order <- c(
  "Arable land",
  "Forests",
  "Heterogeneous\nagricultural areas",
  "Natural grasslands",
  "Open spaces",
  "Pastures",
  "Shrubland",
  "Urban fabric",
  "Waters/Wetlands"
)
levels(euSDm$clcL2) <- group_order

plot_counts$clcL2 <- factor(
  plot_counts$clcL2,
  levels = group_order
)

## Sort groups, then sort classes by abundance within each group
plot_counts <- plot_counts[
  order(plot_counts$clcL2, plot_counts$n),
]

## Set CLC class order
## Largest bar will appear at the top of each facet
plot_counts$clcL3100n <- factor(
  plot_counts$clcL3100n,
  levels = unique(plot_counts$clcL3100n)
)

cols <- c(
  "Arable land" = "#7B3294",
  "Forests" = "#9970AB",
  "Heterogeneous\nagricultural areas" = "#C2A5CF",
  "Natural grasslands" = "#DECBE4",
  "Open spaces" = "#F7F7F7",
  "Pastures" = "#D9F0D3",
  "Shrubland" = "#A6DBA0",
  "Urban fabric" = "#5AAE61",
  "Waters/Wetlands" = "#1B7837"
)

g_bar <- ggplot(
  plot_counts,
  aes(
    x = n,
    y = clcL3100n,
    fill = clcL2
  )
) +
  geom_col(
    width = 0.8,
    colour = "grey45",
    linewidth = 0.3
  ) +
  facet_grid(
    rows = vars(clcL2),
    scales = "free_y",
    space = "free_y",
    labeller = labeller(
      clcL2 = label_wrap_gen(width = 25)
    )
  ) +
  scale_fill_manual(
    values = cols,
    guide = "none"
  ) +
  scale_x_continuous(
    trans = "sqrt",
    breaks = seq(0, 50000, 10000),
    labels = scales::label_number(big.mark = ","),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    x = "Number of EVA plant surveys",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y = element_text(
      size = 9,
      colour = "grey30"
    ),
    axis.text.x = element_text(size = 9),
    axis.title.x = element_text(size = 12),

    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(
      colour = "grey88",
      linewidth = 0.5
    ),

    strip.background = element_rect(
      fill = "grey85",
      colour = "grey40",
      linewidth = 0.5
    ),
    strip.text.y.right = element_text(
      size = 10,
      face = "bold",
      angle = 0,
      margin = margin(5, 8, 5, 8)
    ),

    panel.border = element_rect(
      colour = "grey40",
      fill = NA,
      linewidth = 0.5
    ),

    panel.spacing.y = grid::unit(0.12, "cm"),

    plot.margin = margin(5, 5, 5, 5)
  )

ggsave(file.path(noteapp_dir,"Figure_01_app_barplot.png"), g_bar, dpi = 600, height = 15, width = 10, bg = "white", scale=0.7)

########################
# GSA Austria and Netherlands
nlSD <- as.data.frame(vect( file.path(data_dir, './nlsd_sp_GSA.gpkg')))
atSD <- as.data.frame(vect( file.path(data_dir, './auSD.sp_GSA.gpkg')))

atSD <- atSD[,which(colnames(atSD)%in%colnames(nlSD))]
nlSD <- nlSD[,which(colnames(nlSD)%in%colnames(atSD))]

names(nlSD)
names(atSD)

allSD <- rbind.data.frame(atSD,nlSD)

allSD$clcl2100 <- as.factor(allSD$clcl2100)
allSD$clcl2100 <- droplevels(allSD$clcl2100,exclude="na")
allSD[which(is.na(allSD$translated_name)),]$translated_name <- "Not-Agri"

levels(allSD$clcl2100) <- c("Arable land", "Green\nUrban Areas", "Forests", "Heterogeneous\nAgricultural Areas", "Natural Grasslands", "Open Spaces", "Pastures", "Shrubland", "Urban fabric", "Waters/Wetlands")

allSD$translated_name <- as.factor(allSD$translated_name)
allSD$translated_name <- droplevels(allSD$translated_name)

#Reclass
allSD$clcL2 <- factor(allSD$clcl2100, levels=c("Arable land", "Green\nUrban Areas", "Forests", "Heterogeneous\nAgricultural Areas", "Natural Grasslands", "Open Spaces", "Pastures", "Shrubland", "Urban fabric", "Waters/Wetlands"), ordered=T)
allSD$n <- 1
allSD.a <- aggregate(n~translated_name+clcL2+country, "sum", data=allSD)
allSD.a <-allSD.a[which(allSD.a$n>5),]

allSD.a$GSA <- as.factor(allSD.a$translated_name)
allSD.a$GSA <- droplevels(allSD.a$GSA)
# allSD.a$clcL3100n <- factor(allSD.a$clcL3100n, levels=unique(allSD.a$clcL3100n[order(allSD.a$clcL2)]), order=TRUE)
allSD.a <- allSD.a[allSD.a$clcL2%in%"Pastures",]
allSD.a$GSA <- droplevels(allSD.a$GSA)

levels(allSD.a$GSA) <- c(
  "Alm forage area", 
  "Natural Grassland\nwith main function nature",
  "Permanent\ngrassland",
  "Temporary\ngrassland",
  "Natural Grassland\nwith agri activities",
  "Natural Grassland:\nNature management (CAP)",
  "Litter meadow",
  "Meager pasture",
  "Mowing meadow |\nPasture (3+ uses)",
  "Mowing meadow |\nPasture (2 uses)",
  "Natural areas\n(incl. heath)",
  "Not-agri\ngrassland",
  "One-mown meadow",
  "Permanent pasture")

o <- order(allSD.a$country, partial=allSD.a$n)        # full, not partial, sort
allSD.a <- allSD.a[o, ]                       # reorder rows
allSD.a$GSA <- factor(allSD.a$GSA,levels = unique(allSD.a$GSA[o]), ordered = TRUE)

# Plot using ggplot2 with points and error bars
gg <- ggplot(
  allSD.a,
  aes( x = reorder_within(GSA, n, country),   # ① reorder per country
       y = n,
       fill = country, alpha=0.9)) +
  geom_col(show.legend = FALSE, colour = "grey") +              
  facet_wrap(~ country, scales = "free") +
  scale_x_reordered() +
  scale_fill_manual(values = c("#752982","#1A7736")) +
  scale_y_continuous(labels = scales::label_scientific()) +
  labs(
    x    = "GeoSpatial Application (GSA) land use categories",
    y    = "Number of EVA plant surveys",
    fill = "Country") +
  theme_minimal() +
  theme(
    axis.title   = element_text(size = 16),
    axis.text.x  = element_text(size = 14, angle = 25, hjust = 1),
    axis.text.y  = element_text(size = 16),
    strip.text.x = element_text(size = 18)
  ); gg

ggsave(file.path(noteapp_dir,"Figure_01_app_barplotb.png"), gg, width = 14, height = 6, units = "in", dpi = 300, scale=1.3, bg="white")

# -------------------------------------------------------------------------
# Appendix Figure 3: Environmental effects
# -------------------------------------------------------------------------

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
    geom_smooth(method = "loess", span = 0.9, se = TRUE, colour = "grey30", alpha = 0.25, linewidth=0.2) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
    labs(x = xLab, y = yLab) +
    theme_minimal()
}

####################
## Biogeographical region (iid)
bioGSum <- fitMain$summary.random$bioG
bioGSum[, 2:7] <- exp(bioGSum[, 2:7]) - 1

## Map iid index back to original bioGregions codes
bioGFac <- factor(euDatDf2$bioGregions)
bioGLook <- data.frame(
  ID = seq_along(levels(bioGFac)),
  bioGregions = as.numeric(as.character(levels(bioGFac)))
)

bioGSum$bioGregions <- bioGLook$bioGregions[match(bioGSum$ID, bioGLook$ID)]

bioGLabels <- c(
  `1`  = "Alpine",
  `4`  = "Atlantic",
  `5`  = "Black Sea",
  `6`  = "Boreal",
  `7`  = "Continental",
  `9`  = "Mediterranean",
  `11` = "Pannonian",
  `12` = "Steppic"
)

bioGSum$region <- bioGLabels[as.character(bioGSum$bioGregions)]

## Sample sizes per region (uses the model index, not the original codes)
bioGCounts <- aggregate(n ~ bioG, data = data.frame(bioG = mDat2$bioG, n = 1), sum)
bioGCounts$bioG <- bioGSum$bioGregions
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
  geom_text(aes(label = paste0("n=", samples)), y = 2.6, size = 2) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
  labs(x = "Biogeographical region", y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

####################
## Year (rw2)
yearSum <- summMarg(fitMain$marginals.random$year)
plotYear <- data.frame(
  year = 1990:2021,
  yearSum
)

ggYear <- plotRw(plotYear, "year", "Year")

####################
## Area (rw1; grouped)
areaSum <- summMarg(fitMain$marginals.random$area_g)
areaVals <- tapply(euDatDf2$area, mDat2$area_g, mean)

plotArea <- data.frame(
  area = as.numeric(areaVals),
  areaSum
)

ggArea <- plotRw(plotArea, "area", expression(Vegetation~Plot~Area~(m^2)))

####################
## Temperature (rw1; grouped)
tempSum <- summMarg(fitMain$marginals.random$temp_g)
tempVals <- tapply(euDatDf2$temp, mDat2$temp_g, mean)

plotTemp <- data.frame(
  temp = as.numeric(tempVals),
  tempSum
)

ggTemp <- plotRw(plotTemp, "temp", expression(Mean~Annual~Temperature~("\u00B0"*C)))

####################
## Precipitation (rw1; grouped)
prepSum <- summMarg(fitMain$marginals.random$prep_g)
prepVals <- tapply(euDatDf2$prec, mDat2$prep_g, mean)

plotPrep <- data.frame(
  prec = as.numeric(prepVals),
  prepSum
)

ggPrep <- plotRw(plotPrep, "prec", expression(Cumulative~Precipitation~(mm%.%year^{-1})))

####################
## Elevation (rw1; grouped)
elevSum <- summMarg(fitMain$marginals.random$elev_g)
elevVals <- tapply(euDatDf2$elev, mDat2$elev_g, mean)

plotElev <- data.frame(
  elev = as.numeric(elevVals),
  elevSum
)

ggElev <- plotRw(plotElev, "elev", "Elevation (m)")

####################
## pH (rw1; grouped)
pHSum <- summMarg(fitMain$marginals.random$pH_g)
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
  ncol = 3, nrow = 2,
  labels = c("a","    b","c","d","     e","f")
)
gArr <- annotate_figure(
  gArr,
  left = text_grob(
    "Relative change in 1/Simpson",
    rot = 90,
    size = 14
  )
); plot(gArr)

ggsave(file.path(noteapp_dir,"Figure_03_app_env_fact.png"), gArr, width = 14, height = 6, units = "in", dpi = 300, scale=1.3, bg="white")

####################

# -------------------------------------------------------------------------
# Appendix Figure 6: increase in plant survey from Atlantic biogeographical regions
# -------------------------------------------------------------------------
# Relative abundance of vegetation plots per region per year
regionRel <- data.frame(
  bioGregions = euDatDf2$bioGregions,
  year = euDatDf2$year,
  n = 1
)
regionRel$region <- bioGLabels[as.character(regionRel$bioGregions)]
regionRel$region2 <- ifelse(regionRel$region == "Atlantic", "Atlantic", "Other regions")
# regionRel$region2 <- ifelse(regionRel$region == "Continental", "Atlantic/Continental", "Other regions")

regionSum2 <- aggregate(n ~ region2 + year, data = regionRel, FUN = sum)

yearTot2 <- aggregate(n ~ year, data = regionRel, FUN = sum)
names(yearTot2)[2] <- "nYear"

regionSum2 <- merge(regionSum2, yearTot2, by = "year", all.x = TRUE)
regionSum2$nRel <- 100 * regionSum2$n / regionSum2$nYear

gg_biog <- ggplot(regionSum2, aes(x = year, y = nRel, colour = region2, group = region2)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_colour_viridis_d(option = "viridis", direction = 1) +
  scale_y_continuous(
    name = "Relative contribution of EVA vegetation plots (%)",
    limits = c(0, 100),
    expand = c(0, 0)
  ) +
  scale_x_discrete(
    name = "Sampling year"
    # breaks = pretty(regionSum2$year, n = 10)
  ) +
  labs(colour = NULL) +
    theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.text = element_text(size=15),
    axis.title = element_text(size=15)
  )

ggsave(file.path(noteapp_dir,"Figures_06_app_year_effect.png"), gg_biog, width = 14, height = 6, units = "in", dpi = 300, scale=1.3, bg="white")

###################################################

# -------------------------------------------------------------------------
# Appendix Figure 4: SPDE mesh + spatial field
# -------------------------------------------------------------------------

# First plot and save the SPDE mesh
euSD.sf <- st_as_sf(readRDS(file.path(data_dir, "evaVS_formodel.RDS")))

ggmesh <- ggplot() +
gg(euSpde2$mesh,lwd=0.1) +
# gg(poly.barrier, col="white", alpha=0.5) +
geom_sf(data=euSD.sf,col='purple',size=0.2,alpha=0.1) +
scale_x_continuous(limits=c(2350000, 6900000)) +
scale_y_continuous(limits=c(1201000, 5480000)) +
theme_minimal() +
labs(title = "Delaunay triangulation mesh constraining the spatial field",
  x = "Longitude", y = "Latitude") +
theme(legend.position = "right",
  plot.title = element_text(size = 14, face = "bold"),
  plot.subtitle = element_text(size = 12),
  axis.text = element_text(size = 10)); ggmesh

ggsave(file.path(noteapp_dir,"Figure_02_app_spdeA.png"), ggmesh, dpi = 600, height = 10, width = 10, bg = "white")

###################################################
proj <- inla.mesh.projector(
  euSpde2$mesh,
  xlim = range(euSpde2$mesh$loc[, 1]),
  ylim = range(euSpde2$mesh$loc[, 2]),
  dims = c(2000, 2000)
)

## Spatial field on the mesh (must be length euSpde2$f$n)
u <- fitMain$summary.random$nodes$mean
stopifnot(length(u) == euSpde2$f$n)

## Use the recommended evaluator (new INLA/fmesher)
spatial_values <- fmesher::fm_evaluate(proj, field = u)
spatial_raster <- rast(
  nrows = 2000, ncols = 2000,
  xmin = min(euSpde2$mesh$loc[, 1]), xmax = max(euSpde2$mesh$loc[, 1]),
  ymin = min(euSpde2$mesh$loc[, 2]), ymax = max(euSpde2$mesh$loc[, 2])
)

values(spatial_raster) <- as.vector(spatial_values)
spatial_raster <- flip(spatial_raster, direction = "vertical")

## CRS: prefer europeSimp/europe_sf CRS; left as provided
crs(spatial_raster) <- crs(europeSf)

spatial_raster_clipped <- mask(spatial_raster, vect(europeSf))

spatial_df <- as.data.frame(spatial_raster_clipped, xy = TRUE)
colnames(spatial_df) <- c("x", "y", "effect")

effect_min <- min(spatial_df$effect, na.rm = TRUE)
effect_max <- max(spatial_df$effect, na.rm = TRUE)

ggfield <- ggplot() +
  geom_tile(data = spatial_df, aes(x = x, y = y, col = effect, fill = effect)) +
  geom_sf(data = europeSf, fill = NA, color = "black", linewidth = 0.1) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0, limits = c(effect_min, effect_max),
    name = "Effect Size"
  ) +
  scale_colour_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0, limits = c(effect_min, effect_max),
    name = "Effect Size"
  ) +
  scale_x_continuous(limits = c(2350000, 6900000)) +
  scale_y_continuous(limits = c(1201000, 5480000)) +
  theme_minimal() +
  labs(title = "Estimated Random Spatial Field", x = "Longitude", y = "Latitude") +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold"),
    axis.text  = element_text(size = 10)
  )
ggfield

ggsave(file.path(noteapp_dir,"Figure_02_app_spdeB.png"), ggfield, dpi = 600, height = 10, width = 10, bg = "white")

###################################################
# PC priors for barrier-SPDE hyperparameters (Figure appendix)

U_range <- 50e3
alpha_r <- 0.10
U_sigma <- 1
alpha_s <- 0.10

lambda_r <- -log(1 - alpha_r) / U_range
lambda_s <- -log(alpha_s) / U_sigma

d_pc_range <- function(r) lambda_r * exp(-lambda_r * r)
d_pc_sigma <- function(s) lambda_s * exp(-lambda_s * s)

range_vals <- seq(0, 20 * U_range, length.out = 500)
sigma_vals <- seq(0, 4 * U_sigma, length.out = 500)

df <- rbind(
  data.frame(x = range_vals / 1e3, density = d_pc_range(range_vals), which = "Range (km)"),
  data.frame(x = sigma_vals,       density = d_pc_sigma(sigma_vals), which = "Sigma (σ)")
)

ggpriors <- ggplot(df, aes(x, density)) +
  geom_line(linewidth = 1.1, colour = "steelblue4") +
  geom_vline(
    data = subset(df, which == "Range (km)"),
    aes(xintercept = 50),
    linetype = "dashed", colour = "firebrick", linewidth = 0.7
  ) +
  geom_vline(
    data = subset(df, which == "Sigma (σ)"),
    aes(xintercept = 1),
    linetype = "dashed", colour = "firebrick", linewidth = 0.7
  ) +
  facet_wrap(~which, scales = "free", nrow = 1) +
  labs(
    x = "Parameter value",
    y = "PC-prior density",
    title = "PC priors for barrier-SPDE hyperparameters",
    subtitle = "Range: P(r < 50 km) = 0.10   |   σ: P(σ > 1) = 0.10"
  ) +
  scale_y_continuous(labels = scales::label_scientific()) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 15, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "none"
  )
ggpriors
ggsave(file.path(noteapp_dir,"Figure_02_app_spdeC.png"), ggpriors, dpi = 600, height = 10, width = 10, bg = "white")
