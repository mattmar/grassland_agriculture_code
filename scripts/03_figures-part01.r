### Figures for paper (based on Part 2 model output)
## Assumes Part 1–2 have been run and these objects exist:
## fitMain, euDatDf2, euDatSp2, mDat2, euMesh2, euSpde2, europeSimp (sf) or europeSf (sf)

## Packages (figures only)
library(INLA)
library(ggplot2)
library(sf)
library(terra)
library(tidyterra)
library(scales)
library(gridExtra)
library(viridis)
library(patchwork)
library(reshape2)
library(RColorBrewer)
library(ggdist)
library(ggbreak)
library(ggpubr)
library(dplyr)

###################################################
# Figure 02: Land use (raw data + model estimates)

## Land-use labels (CLC L2)
clcL2Labs <- c(
  "Arable land",
  "Forests",
  "Heterogeneous\nAgricultural Areas",
  "Natural Grasslands",
  "Open Spaces",
  "Pastures",
  "Shrubland",
  "Urban fabric",
  "Waters/Wetlands"
)

## Raw distributions (by land use)
oriData <- stats::aggregate(L ~ clcL2100, data = as.data.frame(euDatSp2), median)
names(oriData) <- c("landUse", "L")
oriData$landUse <- factor(oriData$landUse)
levels(oriData$landUse) <- clcL2Labs
oriData <- oriData[order(oriData$L, decreasing = FALSE), ]
oriData$landUse <- factor(oriData$landUse, levels = oriData$landUse, ordered = TRUE)

rawStats <- data.frame(
  landUse = euDatDf2$clcL2100,
  L       = euDatDf2$L
)
rawStats$landUse <- factor(rawStats$landUse)
levels(rawStats$landUse) <- clcL2Labs
rawStats$landUse <- factor(rawStats$landUse, levels = oriData$landUse, ordered = TRUE)
rawStats <- rawStats[-which(is.na(rawStats$landUse)), ]

medians <- aggregate(L ~ landUse, rawStats, "median")
Omedian <- median(rawStats$L)

## Land-use shares (these are external stats)
clcStats <- structure(
  list(
    landUse = clcL2Labs,
    count   = c(133234356, 171486903, 68844613, 21557169, 35092635, 43061005, 58089995, 24773642, 30174580),
    rela    = c(0.23, 0.29, 0.12, 0.04, 0.06, 0.07, 0.10, 0.04, 0.05)
  ),
  row.names = c(NA, -9L),
  class = "data.frame"
)
clcStats$landUse <- factor(clcStats$landUse, levels = oriData$landUse, ordered = TRUE)

g0 <- ggplot(rawStats, aes(x = L, y = landUse, fill = landUse)) +
  stat_slabinterval(
    justification  = 0.01,
    point_interval = NULL,
    orientation    = "horizontal",
    show.legend    = FALSE,
    slab_color     = "grey50",
    interval_alpha = 0.9,
    slab_linewidth = 0.5,
    slab_alpha     = 0.8
  ) +
  geom_segment(
    data = medians,
    aes(x = L, xend = L, y = as.numeric(landUse), yend = as.numeric(landUse) + 0.50),
    inherit.aes = FALSE,
    colour      = "black",
    linewidth   = 0.8,
    linetype    = 1
  ) +
  geom_vline(
    aes(xintercept = Omedian),
    colour    = "grey50",
    linewidth = 0.5,
    linetype  = 2
  ) +
  geom_text(
    data = clcStats,
    aes(x = 0.2, y = as.numeric(landUse), label = paste0("(", rela * 100, " %)"), size = 8),
    show.legend = FALSE
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.05)), limits = c(0, 45)) +
  scale_fill_viridis_d(option = "magma") +
  scale_x_break(
    breaks     = c(10, 20),
    space      = 0.2,
    scale      = 0.3,
    expand     = TRUE,
    ticklabels = c(25, 35, 45)
  ) +
  labs(
    x = "Grassland diversity (1 / Simpson index)",
    y = "Dominant Land Use Category"
  ) +
  theme_minimal() +
  theme(
    axis.text.y  = element_text(size = 20, lineheight = 1.0),
    axis.text.x  = element_text(size = 16, lineheight = 1.0),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.line.y  = element_blank(),
    axis.line.x  = element_blank()
  )
g0

## Model effects: land-use IID term
clcL2Eff <- fitMain$summary.random$clcL2

clcData <- data.frame(
  clcL2      = factor(clcL2Eff$ID),
  mean       = clcL2Eff$`0.5quant`,
  lower      = clcL2Eff$`0.025quant`,
  upper      = clcL2Eff$`0.975quant`,
  abundance  = as.numeric(summary(factor(mDat2$clcL2))[1:9])
)
levels(clcData$clcL2) <- clcL2Labs

clcData$group <- c(
  "High-Intensity", "Abandonment", "Low-Intensity", "Low-Intensity",
  "Other", "Low-Intensity", "Abandonment", "High-Intensity", "Other"
)

clcData <- clcData[order(clcData$clcL2), ]
clcData$clcL2 <- factor(clcData$clcL2, levels = levels(clcData$clcL2)[order(clcData$mean)], ordered = TRUE)

g1 <- ggplot(clcData, aes(x = clcL2, y = mean, shape = group, fill = ifelse(mean > 0, "#276419", "#8e0152"))) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.4, color = "grey") +
  geom_point(size = 10, stroke = 1.2, show.legend = TRUE) +
  geom_hline(yintercept = 0, lty = 2, col = "black") +
  scale_shape_manual(
    values = c(
      "Low-Intensity"  = 21,
      "High-Intensity" = 24,
      "Abandonment"    = 23,
      "Other"          = 22
    ),
    name = "Land Use Process"
  ) +
  scale_fill_viridis_d(option = "magma", guide = "none", direction = -1, begin = 0.15, end = 0.95) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Marginal effect on Grassland Plant Diversity (1 / Simpson)"
  ) +
  scale_y_continuous(labels = scales::label_percent()) +
  theme_minimal() +
  theme(
    axis.text.x       = element_text(size = 16),
    axis.text.y       = element_text(size = 18, margin = margin(r = -25)),
    axis.title.x      = element_text(size = 18),
    legend.position   = c(0.7, 0.25),
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
    legend.text       = element_text(size = 16),
    legend.title      = element_text(size = 18),
    legend.key.height = unit(2.5, "lines")
  ) +
  guides(shape = guide_legend(override.aes = list(size = 8)))
g1

axisBreaks <- seq(0, 50e3, 10e3)
axisLabels <- format(seq(0, 50e3, 10e3), scientific = TRUE)
axisLabels[1] <- ""

g2 <- ggplot(clcData) +
  geom_col(
    aes(x = clcL2, y = abundance, fill = abundance),
    position = position_dodge2(width = 0),
    width    = 7.6,
    show.legend = FALSE
  ) +
  scale_fill_gradient2(low = "white", mid = "grey", high = "black") +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.1)),
    breaks = axisBreaks,
    labels = NULL
  ) +
  theme_minimal() +
  theme(
    axis.title         = element_blank(),
    axis.text.y        = element_blank(),
    axis.text.x        = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin        = unit(c(-1.3, 0, 0, 0), "cm")
  ) +
  coord_flip() +
  annotate("text", x = 0.9, y = axisBreaks, label = axisLabels, size = 3, color = "black") +
  annotate("text", x = 0.6, y = 27.5e3, label = "Number of Surveys", size = 5, col = "black")
g2

## Export
g0Clean <- g0
g1Clean <- g1 + theme(axis.title.y = element_blank())
g2Clean <- g2 + theme(axis.title.y = element_blank())

layout <- matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE)

gArr <- grid.arrange(
  arrangeGrob(g1Clean, g2Clean, layout_matrix = layout, widths = c(4, 1), heights = c(1))
)

g2g <- arrangeGrob(
  gArr,
  left    = grid::textGrob(NULL, rot = 90, gp = grid::gpar(fontsize = 15)),
  padding = unit(1, "line")
)

g0Clean
plot(g2g)

###################################################
# Figure 03: Nitrogen effects (raw + posterior)

## Joint posterior samples of fixed effects
set.seed(2025)
postSamps <- inla.posterior.sample(n = 1000, result = fitMain, num.threads = 10)

feNames  <- rownames(fitMain$summary.fixed)
feLatent <- paste0(feNames, ":1")
betaMat  <- sapply(postSamps, function(ps) as.data.frame(ps$latent)[feLatent, ])
rownames(betaMat) <- feNames

## Scaling parameters
NDepMean <- attr(scale(euDatDf2$NDEP),   "scaled:center")
NDepSd   <- attr(scale(euDatDf2$NDEP),   "scaled:scale")
NBudMean <- attr(scale(euDatDf2$NFIELD), "scaled:center")
NBudSd   <- attr(scale(euDatDf2$NFIELD), "scaled:scale")

zeroBudScaled <- (0 - NBudMean) / NBudSd
zeroDepScaled <- (0 - NDepMean) / NDepSd

## Grids in scaled space (parameters unchanged)
n <- 1000
gD <- data.frame(
  NDep = seq(zeroDepScaled, max(mDat2$NDep, na.rm = TRUE), length = n),
  NBud = zeroBudScaled,
  cov  = "N deposition"
)
gB <- data.frame(
  NDep = zeroDepScaled,
  NBud = seq(zeroBudScaled, max(mDat2$NBud, na.rm = TRUE), length = n),
  cov  = "N fertiliser"
)
grid <- bind_rows(gD, gB) %>% mutate(intercept = 1)

X <- as.matrix(grid[, c("intercept","NBud", "NDep")])

pctCi <- function(X, X0Row, beta) {
  X0 <- matrix(rep(X0Row, each = nrow(X)), nrow = nrow(X))
  eta  <- X  %*% beta
  eta0 <- X0 %*% beta
  L  <- exp(eta)
  L0 <- exp(eta0)
  pct <- 100 * (L / L0 - 1)
  q <- t(apply(pct, 1, quantile, c(0.25, 0.5, 0.75)))
  colnames(q) <- c("lower", "median", "upper")
  q
}

X0Real <- c(
  intercept=1,
  NBud = (0 - NBudMean) / NBudSd,
  NDep = (0 - NDepMean) / NDepSd
)

ciMat <- pctCi(X, X0Real, betaMat)

plotDf <- bind_cols(grid, as.data.frame(ciMat)) %>%
  mutate(
    x = case_when(
      cov == "N deposition" ~ NDep * NDepSd + NDepMean,
      cov == "N fertiliser" ~ NBud * NBudSd + NBudMean
    )
  )

p1 <- ggplot(plotDf, aes(x = x, y = median, colour = cov, fill = cov)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, colour = NA, show.legend = FALSE) +
  geom_line(linewidth = 1, show.legend = FALSE) +
  facet_wrap(
    ~ cov, scales = "free_x",
    labeller = labeller(cov = c(
      `N deposition` = "Nitrogen deposition",
      `N fertiliser` = "Nitrogen fertiliser"
    ))
  ) +
  labs(
    x   = "Nitrogen input (kg ha⁻¹ yr⁻¹)",
    y   = "Percent change in Simpson Diversity (1/L)",
    tag = "b)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.tag.position = "topright",
    legend.title.position = "top",
    legend.position = "top",
    legend.key.spacing.x = unit(0.5, "cm"),
    legend.spacing = unit(0, "cm"),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    plot.title = element_text(size = 12)
  )

## Raw data (N vs L)
rawStatsN <- data.frame(NDep = euDatDf2$NDEP, NBud = euDatDf2$NFIELD, L = euDatDf2$L)
rawStatsNm <- reshape2::melt(rawStatsN, measure.vars = c("NDep", "NBud"))

labelSource <- c(NDep = "N deposition", NBud = "N fertiliser")

p2 <- ggplot(rawStatsNm, aes(y = L + 1e-6, x = value)) +
  stat_binhex(bins = 35, show.legend = TRUE) +
  facet_wrap(~ variable, scale = "free_x", labeller = labeller(variable = labelSource)) +
  viridis::scale_fill_viridis(
    option = "magma",
    name = "Count\n(log scale)",
    labels = function(x) round(x, 0),
    trans = "log10",
    limits = c(1, NA)
  ) +
  labs(
    x   = "Nitrogen input (kg ha⁻¹ yr⁻¹)",
    y   = "Simpson Diversity (1/L)",
    tag = "a)"
  ) +
  stat_cor(
    method = "pearson",
    p.accuracy = 0.001, r.accuracy = 0.01,
    label.x.npc = "middle",
    label.y.npc = 0.99,
    size = 4,
    color = "#0f993d"
  ) +
  geom_smooth(method = "lm", color = "#0f993d", linetype = "dashed", se = FALSE) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = c(0.45, 0.60),
    legend.background = element_rect(fill = alpha("white", 0.7), colour = NA),
    legend.key.width = unit(0.4, "cm"),
    legend.key.height = unit(0.4, "cm"),
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 8),
    axis.text.y  = element_text(size = 14),
    axis.text.x  = element_text(size = 14),
    axis.ticks.x = element_blank(),
    axis.title   = element_text(size = 14),
    plot.title   = element_text(size = 10),
    strip.text.x = element_text(size = 18),
    plot.tag.position = "topright"
  )

commonX <- grid::textGrob(
  "Nitrogen input (kg ha⁻¹ yr⁻¹)",
  gp = grid::gpar(fontsize = 16)
)

gArr <- gridExtra::grid.arrange(
  gridExtra::arrangeGrob(
    p2 + theme(axis.title.x = element_blank()),
    p1 + theme(axis.title.x = element_blank()),
    ncol = 1, nrow = 2, heights = c(1.2, 1)
  ),
  bottom = commonX
)

plot(gArr)

###################################################
# Figure 04: Barplot (country averages previously calculated)
predNuts <- read.csv("./data/nuts_averages.csv")
names(predNuts)[5] <- "Percent_decrease"
predNuts <- predNuts[complete.cases(predNuts), ]
idx <- order(predNuts$Percent_decrease, decreasing = TRUE)
predNuts$CNTR_CODE <- factor(predNuts$CNTR_CODE, levels = unique(predNuts$CNTR_CODE)[idx], order = TRUE)

g1 <- ggplot(predNuts, aes(x = CNTR_CODE, y = Percent_decrease, fill = Percent_decrease)) +
  geom_bar(aes(x = CNTR_CODE, y = Percent_decrease), stat = "identity", width = 0.95, show.legend = FALSE) +
  geom_hline(
    data = data.frame(y = quantile(predNuts$Percent_decrease, p = c(0.5))),
    aes(yintercept = y),
    lty = 2,
    col = c("red")
  ) +
  scale_y_reverse() +
  geom_label(
    aes(label = CNTR_CODE),
    fontface = "bold",
    fill = "white",
    colour = "black",
    linewidth = 0,
    label.padding = unit(0.05, "lines"),
    label.r = unit(0.05, "lines"),
    size = 5
  ) +
  labs(y = "Decrease in plant diversity (1/Simpson)", x = NULL) +
  scale_fill_viridis_c(option = "magma", direction = 1) +
  scale_colour_viridis_c(option = "magma", direction = 1) +
  theme_minimal() +
  theme(
    axis.text.x  = element_blank(),
    axis.text.y  = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid       = element_blank()
  )
g1
###################################################
# Figure A2 (appendix): Land-use share (CLC L3 within CLC L2)
euSD <- euDatSp

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

euSDm$n <- 1
euSDa <- aggregate(n ~ clcL3100n + clcL2, "sum", data = euSDm)
euSDa <- euSDa[which(euSDa$n > 20), ]

euSDa$clcL3100n <- factor(euSDa$clcL3100n)
euSDa$clcL3100n <- factor(euSDa$clcL3100n, levels = unique(euSDa$clcL3100n[order(euSDa$clcL2)]), order = TRUE)

gg <- ggplot(euSDa, aes(x = clcL3100n, y = n, fill = clcL2)) +
  geom_col(col = "grey", position = position_dodge2(width = 2.5), width = 1) +
  labs(
    fill = "CLC level 3",
    x = "Aggregated CLC categories",
    y = "Number of EVA plant surveys"
  ) +
  scale_y_continuous(labels = scales::label_scientific()) +
  facet_grid(~clcL2, switch = "x", scale = "free_x") +
  scale_fill_brewer(palette = "PRGn") +
  theme_minimal() +
  theme(
    axis.title   = element_text(size = 16),
    strip.text.x = element_text(size = 15),
    axis.text.x  = element_text(size = 12, angle = 55, hjust = 1, vjust = 1),
    axis.text.y  = element_text(size = 16)
  )
gg
###################################################
# Figure 01: Overview map
europe_sf <- st_read(".data/nuts0_simple_dissolved.gpkg")
gridS <- st_read(".data/euSD.spVOver25km_regular.gpkg")
hrl25 <- rast(".data/grassland_coverage25_updated3035.tif")
hrl25.dt <- as.data.frame(hrl25, xy = TRUE)

bbox <- st_bbox(europe_sf)
lat_limits <- c(bbox["ymin"], bbox["ymax"])
lon_limits <- c(bbox["xmin"], bbox["xmax"])
lat_limits[1] <- 1395000
lat_limits[2] <- 5150000
lon_limits[1] <- 2600000
lon_limits[2] <- 6400000

quantile_breaks <- quantile(gridS$L_50, probs = seq(0, 1, length.out = 10), na.rm = TRUE)

breaks <- sort(unique(as.numeric(quantile_breaks)))
breaks[length(breaks)] <- breaks[length(breaks)] + .Machine$double.eps
labels <- paste0(round(head(breaks, -1), 0), "–", round(tail(breaks, -1), 0))

gridS <- gridS |>
  dplyr::mutate(
    category = cut(
      L_50,
      breaks = breaks,
      include.lowest = TRUE,
      right = TRUE,
      labels = labels
    )
  )

gridS_diss <- gridS %>%
  group_by(category) %>%
  summarise(geometry = sf::st_union(geom), .groups = "drop")

main_map <- ggplot() +
  geom_sf(data = europe_sf, fill = "white", color = "grey90", size = 0.3) +
  geom_tile(data = hrl25.dt, aes(x = x, y = y), fill = "grey85", col = NA) +
  geom_sf(data = gridS_diss, aes(fill = category), colour = NA, linewidth = 0) +
  scale_fill_viridis_d(
    direction = 1,
    option = "viridis",
    name = "1/Simpson\nDiversity",
    limits = labels,
    drop = FALSE
  ) +
  guides(fill = "none") +
  coord_sf(xlim = lon_limits, ylim = lat_limits, expand = TRUE) +
  labs(caption = "Data: European Vegetation Archive (EVA) 1990-2023") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    legend.position = c(0.11, 0.25),
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = alpha("transparent", 0.7), colour = NA),
    legend.key = element_rect(fill = "white", colour = NA),
    plot.caption = element_text(size = 8, hjust = 1, color = "grey30", face = "italic")
  )

breaksH <- breaks
breaksH[which(breaksH > 20)] <- 20
y_mid <- (head(breaksH, -1) + tail(breaksH, -1)) / 2

histogram_inset <- ggplot(gridS[which(gridS$L_50 < 20), ], aes(x = L_50)) +
  stat_bin(
    closed = "right",
    bins = 50,
    show.legend = FALSE,
    aes(fill = after_stat(cut(..x.., breaks = breaksH, include.lowest = TRUE, pad = TRUE, labels = labels))),
    color = NA
  ) +
  labs(title = "Legend", subtitle = "1/Simpson diversity", y = NULL, x = NULL) +
  scale_fill_viridis_d(option = "viridis", direction = 1, limits = labels, drop = FALSE) +
  scale_x_continuous(breaks = y_mid, expand = c(0, 0), labels = labels) +
  theme_minimal(base_size = 12) +
  coord_flip() +
  theme(
    panel.background = element_rect(fill = NA, color = NA),
    plot.background  = element_rect(fill = NA, color = NA),
    plot.title       = element_text(size = 10, face = "bold", margin = margin(t = 0, b = 5)),
    plot.subtitle    = element_text(size = 8, margin = margin(t = 0, b = 0)),
    legend.background = element_rect(fill = NA, color = NA),
    panel.grid       = element_blank(),
    legend.position  = "none",
    axis.text.y      = element_text(angle = 0, margin = margin(t = 1), size = 7),
    axis.text.x      = element_text(size = 7)
  )

int_plot <- main_map +
  inset_element(histogram_inset, left = 0.75, bottom = 0.60, right = 0.99, top = 1)

int_plot

###################################################
# Figure B1 (appendix): SPDE mesh + spatial field
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