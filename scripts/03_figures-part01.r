# Figures of the main text
# Run from the repository root. Model/data objects are loaded from ./data
# and figures are written to ./Figures.

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

baseDir <- normalizePath("../git/grassland_agriculture_code/")
data_dir <- file.path(baseDir, "data")
fig_dir  <- file.path(baseDir, "Figures")

dir.create(fig_dir, showWarnings = FALSE)

fitMain <- readRDS(file.path(data_dir, "fitMain.RDS"))
euDatDf2 <- as.data.frame(readRDS(file.path(data_dir, "evaVS_formodel.RDS")))
mDat2 <- read.csv(file.path(data_dir, "mDat2.csv"))

# -------------------------------------------------------------------------
# Figure 1: Overview map
# -------------------------------------------------------------------------

europe_sf <- st_read(file.path(data_dir, "nuts0_simple_dissolved.gpkg"), quiet = TRUE)
gridS <- st_read(file.path(data_dir, "euSD.spVOver25km_regular.gpkg"), quiet = TRUE)
hrl25 <- rast(file.path(data_dir, "grassland_coverage25_updated3035.tif"))
hrl25_df <- as.data.frame(hrl25, xy = TRUE)

x_limits <- c(2600000, 6400000)
y_limits <- c(1395000, 5150000)

# Quantile classes of observed Simpson diversity
breaks <- quantile(
  gridS$L_50,
  probs = seq(0, 1, length.out = 10),
  na.rm = TRUE
)
breaks <- sort(unique(as.numeric(breaks)))
breaks[length(breaks)] <- breaks[length(breaks)] + .Machine$double.eps

class_labels <- paste0(
  round(head(breaks, -1), 0),
  "–",
  round(tail(breaks, -1), 0)
)

gridS <- gridS |>
  mutate(
    category = cut(
      L_50,
      breaks = breaks,
      include.lowest = TRUE,
      right = TRUE,
      labels = class_labels
    )
  )

# Dissolve grid cells within each diversity class
gridS_diss <- gridS |>
  group_by(category) |>
  summarise(.groups = "drop")

main_map <- ggplot() +
  geom_sf(data = europe_sf, fill = "white", color = "grey90", linewidth = 0.3) +
  geom_tile(data = hrl25_df, aes(x = x, y = y), fill = "grey80", color = NA) +
  geom_sf(data = gridS_diss, aes(fill = category), colour = NA, linewidth = 0) +
  scale_fill_viridis_d(
    direction = 1,
    option = "viridis",
    limits = class_labels,
    drop = FALSE,
    guide = "none"
  ) +
  coord_sf(xlim = x_limits, ylim = y_limits, expand = TRUE) +
  theme_minimal() +
  theme(axis.title = element_blank())

# Inset distribution, truncated at Simpson diversity = 20
breaks_hist <- breaks
breaks_hist[breaks_hist > 20] <- 20
breaks_hist <- unique(breaks_hist)

class_midpoints <- (head(breaks_hist, -1) + tail(breaks_hist, -1)) / 2
hist_labels <- paste0(
  format(round(head(breaks_hist, -1), 1), trim = TRUE),
  "–",
  format(round(tail(breaks_hist, -1), 1), trim = TRUE)
)

histogram_inset <- gridS |>
  filter(L_50 < 20) |>
  ggplot(aes(x = L_50)) +
  stat_bin(
    bins = 50,
    closed = "right",
    aes(
      fill = after_stat(
        cut(
          x,
          breaks = breaks_hist,
          include.lowest = TRUE,
          labels = hist_labels
        )
      )
    ),
    color = NA
  ) +
  scale_fill_viridis_d(
    option = "viridis",
    direction = 1,
    limits = hist_labels,
    drop = FALSE,
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = class_midpoints,
    labels = hist_labels,
    minor_breaks = NULL,
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_continuous(
    breaks = scales::breaks_pretty(n = 3),
    labels = scales::label_number(big.mark = ","),
    expand = expansion(mult = c(0, 0.03))
  ) +
  coord_flip() +
  labs(
    title = "Simpson diversity (1/D)\nquantile classes",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 9) +
  theme(
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 7),
    plot.title = element_text(
      size = 9,
      face = "bold",
      lineheight = 1,
      margin = margin(b = 4)
    ),
    plot.title.position = "plot",
    plot.margin = margin(1, 1, 1, 1),
    plot.background = element_rect(fill = NA, colour = NA),
    panel.background = element_rect(fill = NA, colour = NA)
  )

figure_01 <- main_map +
  inset_element(
    histogram_inset,
    left = 0.75,
    bottom = 0.60,
    right = 0.99,
    top = 1
  )

ggsave(
  file.path(fig_dir, "Figure_01.pdf"),
  plot = figure_01,
  width = 21,
  height = 21,
  units = "cm",
  device = cairo_pdf
)


# -------------------------------------------------------------------------
# Figure 2: Nitrogen effects
# -------------------------------------------------------------------------

# Draw joint posterior samples of fixed effects
set.seed(2025)
post_samps <- inla.posterior.sample(
  n = 1000,
  result = fitMain,
  num.threads = 30
)

fe_names <- rownames(fitMain$summary.fixed)
fe_latent <- paste0(fe_names, ":1")

beta_mat <- sapply(
  post_samps,
  function(ps) as.data.frame(ps$latent)[fe_latent, ]
)
rownames(beta_mat) <- fe_names

# Scaling parameters used in the fitted model
NDep_mean <- attr(scale(euDatDf2$NDEP), "scaled:center")
NDep_sd   <- attr(scale(euDatDf2$NDEP), "scaled:scale")
NBud_mean <- attr(scale(euDatDf2$NFIELD), "scaled:center")
NBud_sd   <- attr(scale(euDatDf2$NFIELD), "scaled:scale")

zero_bud_scaled <- (0 - NBud_mean) / NBud_sd
zero_dep_scaled <- (0 - NDep_mean) / NDep_sd

# Prediction grids in model-scaled space
n_grid <- 1000

grid_dep <- data.frame(
  NDep = seq(zero_dep_scaled, max(mDat2$NDep, na.rm = TRUE), length.out = n_grid),
  NBud = zero_bud_scaled,
  cov = "Atmospheric N deposition"
)

grid_bud <- data.frame(
  NDep = zero_dep_scaled,
  NBud = seq(zero_bud_scaled, max(mDat2$NBud, na.rm = TRUE), length.out = n_grid),
  cov = "Fertiliser N input"
)

pred_grid <- bind_rows(grid_dep, grid_bud) |>
  mutate(intercept = 1)

X <- as.matrix(pred_grid[, c("intercept", "NBud", "NDep")])

# Posterior median and interquartile range of percentage change from zero N
posterior_pct_interval <- function(X, X0_row, beta) {
  X0 <- matrix(rep(X0_row, each = nrow(X)), nrow = nrow(X))

  eta <- X %*% beta
  eta0 <- X0 %*% beta

  L <- exp(eta)
  L0 <- exp(eta0)

  pct_change <- 100 * (L / L0 - 1)

  q <- t(apply(
    pct_change,
    1,
    quantile,
    probs = c(0.25, 0.50, 0.75)
  ))
  colnames(q) <- c("lower", "median", "upper")
  q
}

X0_zero_n <- c(
  intercept = 1,
  NBud = zero_bud_scaled,
  NDep = zero_dep_scaled
)

posterior_interval <- posterior_pct_interval(X, X0_zero_n, beta_mat)

plot_df <- bind_cols(pred_grid, as.data.frame(posterior_interval)) |>
  mutate(
    x = case_when(
      cov == "Atmospheric N deposition" ~ NDep * NDep_sd + NDep_mean,
      cov == "Fertiliser N input" ~ NBud * NBud_sd + NBud_mean
    )
  )

p_posterior <- ggplot(
  plot_df,
  aes(x = x, y = median, colour = cov, fill = cov)
) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.2,
    colour = NA,
    show.legend = FALSE
  ) +
  geom_line(linewidth = 1, show.legend = FALSE) +
  facet_wrap(
    ~ cov,
    scales = "free_x",
    labeller = labeller(
      cov = c(
        `Atmospheric N deposition` = "Atmospheric N deposition",
        `Fertiliser N input` = "Fertiliser N input"
      )
    )
  ) +
  labs(
    x = "Nitrogen input (kg ha⁻¹ yr⁻¹)",
    y = "Change in Simpson diversity relative to zero N (%)",
    tag = "b"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.tag.position = "topleft",
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    plot.title = element_text(size = 12),
    plot.tag = element_text(face = "bold", size = 14)
  )

# Raw diversity observations and Pearson correlations
raw_stats_n <- data.frame(
  NDep = euDatDf2$NDEP,
  NBud = euDatDf2$NFIELD,
  L = euDatDf2$L
)

raw_stats_long <- reshape2::melt(
  raw_stats_n,
  measure.vars = c("NDep", "NBud")
)

nitrogen_labels <- c(
  NDep = "Atmospheric N deposition",
  NBud = "Fertiliser N input"
)

p_raw <- ggplot(raw_stats_long, aes(y = L + 1e-6, x = value)) +
  stat_binhex(bins = 35, show.legend = TRUE) +
  facet_wrap(
    ~ variable,
    scales = "free_x",
    labeller = labeller(variable = nitrogen_labels)
  ) +
  viridis::scale_fill_viridis(
    option = "magma",
    name = "Count\n(log scale)",
    labels = function(x) round(x, 0),
    trans = "log10",
    limits = c(1, NA)
  ) +
  labs(
    x = "Nitrogen input (kg ha⁻¹ yr⁻¹)",
    y = "Simpson Diversity",
    tag = "a"
  ) +
  stat_cor(
    method = "pearson",
    cor.coef.name = "r",
    p.accuracy = 0.001,
    r.accuracy = 0.01,
    label.x.npc = "middle",
    label.y.npc = 0.99,
    size = 4,
    color = "#0f993d"
  ) +
  geom_smooth(
    method = "lm",
    color = "#0f993d",
    linetype = "dashed",
    se = FALSE,
    linewidth = 1.2
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = c(0.47, 0.50),
    legend.background = element_rect(fill = alpha("white", 0), colour = NA),
    legend.key.width = unit(0.4, "cm"),
    legend.key.height = unit(0.4, "cm"),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 14),
    strip.text.x = element_text(size = 16),
    plot.tag.position = "topleft",
    plot.tag = element_text(face = "bold", size = 14)
  )

common_x <- grid::textGrob(
  "Nitrogen input (kg ha⁻¹ yr⁻¹)",
  gp = grid::gpar(fontsize = 12)
)

figure_02 <- gridExtra::grid.arrange(
  gridExtra::arrangeGrob(
    p_raw +
      theme(
        axis.title.x = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)
      ),
    p_posterior +
      theme(
        axis.title.x = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)
      ),
    ncol = 1,
    heights = c(1.2, 1)
  ),
  bottom = common_x
)

ggsave(
  file.path(fig_dir, "Figure_02.pdf"),
  plot = figure_02,
  width = 18,
  height = 14,
  units = "cm",
  device = cairo_pdf,
  scale = 1.5
)


# -------------------------------------------------------------------------
# Figure 3: Land use
# -------------------------------------------------------------------------

land_use_labels <- c(
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

# Order land-use categories by observed median diversity
land_use_medians <- stats::aggregate(L ~ clcL2100, data = euDatDf2, median)
names(land_use_medians) <- c("landUse", "L")

land_use_medians$landUse <- factor(land_use_medians$landUse)
levels(land_use_medians$landUse) <- land_use_labels

land_use_medians <- land_use_medians[order(land_use_medians$L), ]
land_use_medians$landUse <- factor(
  land_use_medians$landUse,
  levels = land_use_medians$landUse,
  ordered = TRUE
)

raw_land_use <- data.frame(
  landUse = euDatDf2$clcL2100,
  L = euDatDf2$L
)

raw_land_use$landUse <- factor(raw_land_use$landUse)
levels(raw_land_use$landUse) <- land_use_labels
raw_land_use$landUse <- factor(
  raw_land_use$landUse,
  levels = land_use_medians$landUse,
  ordered = TRUE
)
raw_land_use <- raw_land_use[!is.na(raw_land_use$landUse), ]

medians <- aggregate(L ~ landUse, raw_land_use, median)
overall_median <- median(raw_land_use$L)

# Relative land-use shares shown beside each distribution
clc_stats <- data.frame(
  landUse = land_use_labels,
  rela = c(0.23, 0.29, 0.12, 0.04, 0.06, 0.07, 0.10, 0.04, 0.05)
)
clc_stats$landUse <- factor(
  clc_stats$landUse,
  levels = land_use_medians$landUse,
  ordered = TRUE
)

p_land_raw <- ggplot(
  raw_land_use,
  aes(x = L, y = landUse, fill = landUse)
) +
  stat_slabinterval(
    justification = 0.01,
    point_interval = NULL,
    orientation = "horizontal",
    show.legend = FALSE,
    slab_color = "grey50",
    interval_alpha = 0.9,
    slab_linewidth = 0.5,
    slab_alpha = 0.8
  ) +
  geom_segment(
    data = medians,
    aes(
      x = L,
      xend = L,
      y = as.numeric(landUse),
      yend = as.numeric(landUse) + 0.50
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.8
  ) +
  geom_vline(
    xintercept = overall_median,
    colour = "grey50",
    linewidth = 0.5,
    linetype = 2
  ) +
  geom_text(
    data = clc_stats,
    aes(
      x = 0,
      y = as.numeric(landUse),
      label = paste0(rela * 100, "%")
    ),
    size = 4,
    show.legend = FALSE
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.1, 0.06)),
    limits = c(0, 45)
  ) +
  scale_fill_viridis_d(option = "magma") +
  scale_x_break(
    breaks = c(10, 20),
    space = 0.2,
    scale = 0.3,
    expand = TRUE,
    ticklabels = c(25, 35, 45)
  ) +
  labs(
    x = "Grassland diversity (1 / Simpson index)",
    y = "Dominant Land Use Category"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 14, lineheight = 3.0),
    axis.text.x = element_text(size = 14, lineheight = 3.0),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.line = element_blank()
  )

# Posterior medians and 95% credible intervals for land-use random effects
clc_effect <- fitMain$summary.random$clcL2

clc_data <- data.frame(
  clcL2 = factor(clc_effect$ID),
  median = clc_effect$`0.5quant`,
  lower = clc_effect$`0.025quant`,
  upper = clc_effect$`0.975quant`,
  abundance = as.numeric(summary(factor(mDat2$clcL2))[1:9])
)
levels(clc_data$clcL2) <- land_use_labels

clc_data$group <- c(
  "High-intensity",
  "Abandonment   ",
  "Low-intensity",
  "Low-intensity",
  "Other",
  "Low-intensity",
  "Abandonment   ",
  "High-intensity",
  "Other"
)

clc_data <- clc_data[order(clc_data$clcL2), ]
clc_data$clcL2 <- factor(
  clc_data$clcL2,
  levels = levels(clc_data$clcL2)[order(clc_data$median)],
  ordered = TRUE
)

p_land_effect <- ggplot(
  clc_data,
  aes(
    x = clcL2,
    y = median,
    shape = group,
    fill = ifelse(median > 0, "#276419", "#8e0152")
  )
) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.6,
    color = "grey70"
  ) +
  geom_point(size = 8, stroke = 1, show.legend = TRUE) +
  geom_hline(yintercept = 0, linetype = 2, colour = "black") +
  scale_shape_manual(
    values = c(
      "Low-intensity" = 21,
      "High-intensity" = 24,
      "Abandonment   " = 23,
      "Other" = 22
    ),
    name = "Land-use process:"
  ) +
  scale_fill_viridis_d(
    option = "magma",
    guide = "none",
    direction = -1,
    begin = 0.15,
    end = 0.95
  ) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Marginal effect on grassland plant diversity (%)"
  ) +
  scale_y_continuous(labels = scales::label_percent()) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14, margin = margin(r = -25)),
    axis.title.x = element_text(size = 16),
    legend.position = "top",
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    legend.key.height = unit(0.5, "lines"),
    plot.margin = unit(c(0, 0.2, 0, -0.5), "cm")
  ) +
  guides(
    shape = guide_legend(
      nrow = 2,
      byrow = TRUE,
      override.aes = list(size = 6)
    )
  )

axis_breaks <- seq(0, 50e3, 10e3)
axis_labels <- paste0(seq(0, 50, 10), "k")
axis_labels[1] <- ""

p_land_n <- ggplot(clc_data) +
  geom_col(
    aes(x = clcL2, y = abundance, fill = abundance),
    position = position_dodge2(width = 0),
    width = 8,
    show.legend = FALSE
  ) +
  scale_fill_gradient2(low = "white", mid = "grey", high = "black") +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.1)),
    breaks = axis_breaks,
    labels = NULL
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = unit(c(0.8, 0, 0, 0.1), "cm")
  ) +
  coord_flip() +
  annotate(
    "text",
    x = 0.8,
    y = axis_breaks,
    label = axis_labels,
    size = 5,
    color = "grey30"
  ) +
  annotate(
    "text",
    x = 0.54,
    y = 27.5e3,
    label = "Number of vegetation plots",
    size = 5.5,
    color = "black"
  )

p_land_raw <- p_land_raw +
  labs(tag = "a") +
  theme(
    plot.tag = element_text(face = "bold", size = 14),
    plot.tag.position = c(0.01, 0.99)
  )

p_land_effect <- p_land_effect +
  labs(tag = "b") +
  theme(
    axis.title.y = element_blank(),
    plot.tag = element_text(face = "bold", size = 14),
    plot.tag.position = c(0.01, 0.98)
  )

p_land_n <- p_land_n +
  labs(tag = "c") +
  theme(
    axis.title.y = element_blank(),
    plot.tag = element_text(face = "bold", size = 14),
    plot.tag.position = c(0.01, 1.021)
  )

figure_03 <- arrangeGrob(
  p_land_raw,
  p_land_effect,
  p_land_n,
  nrow = 1,
  widths = c(2.5, 2, 1)
)

ggsave(
  file.path(fig_dir, "Figure_03.pdf"),
  plot = figure_03,
  width = 21,
  height = 10,
  units = "cm",
  device = cairo_pdf,
  scale = 2.2
)

# -------------------------------------------------------------------------
# Figure 4: Country-level modelled diversity deficits
# -------------------------------------------------------------------------

pred_nuts <- read.csv(file.path(data_dir, "nuts_averages.csv"))
names(pred_nuts)[5] <- "Percent_decrease"

pred_nuts <- pred_nuts[complete.cases(pred_nuts), ]
pred_nuts$Diversity_deficit <- pred_nuts$Percent_decrease
pred_nuts <- pred_nuts[
  order(pred_nuts$Diversity_deficit),
]

pred_nuts$CNTR_CODE <- factor(
  pred_nuts$CNTR_CODE,
  levels = pred_nuts$CNTR_CODE
)

median_deficit <- median(pred_nuts$Diversity_deficit, na.rm = TRUE)

figure_04 <- ggplot(
  pred_nuts,
  aes(
    y = CNTR_CODE,
    x = Diversity_deficit,
    fill = Diversity_deficit
  )
) +
  geom_col(width = 0.8, show.legend = FALSE) +
  geom_vline(
    xintercept = median_deficit,
    linetype = 2,
    linewidth = 0.6,
    colour = "grey35"
  ) +
  geom_text(
    aes(label = sprintf("%.1f", Diversity_deficit)),
    hjust = -0.15,
    size = 3
  ) +
  scale_fill_viridis_c(
    option = "magma",
    direction = -1
  ) +
  scale_x_continuous(
    breaks = seq(0, 30, 5),
    labels = label_number(suffix = "%"),
    limits = c(0, 32),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "Modelled diversity deficit under current N inputs (%)",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 9, face = "bold"),
    axis.text.x = element_text(size = 9),
    axis.title.x = element_text(size = 11),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(5, 5, 5, 5)
  )

ggsave(
  file.path(fig_dir, "Figure_04barplot.png"),
  plot = figure_04,
  width = 15,
  height = 15,
  units = "cm",
  scale = 1.5
)

# -------------------------------------------------------------------------
# Figure 5: Species-level nitrogen responses
# -------------------------------------------------------------------------

fitted_data <- readRDS(file.path(data_dir, "nitrogen_species_fitted_output.RDS"))

# Retain fitted rows used in the species-level figure
fitted_data <- fitted_data[fitted_data$N.pval < 0.99, ]

eive_data <- read.csv(
  file.path(data_dir, "JD341_Dengler_et_al_2023_VCS_EIVE_1.0_SM8_EIVE.csv")
)

fitted_data_sub <- merge(
  fitted_data,
  eive_data[, c("TaxonConcept", "EIVEres.N", "EIVEres.N.nw3")],
  by.x = "Turboveg2.simple",
  by.y = "TaxonConcept"
)

fitted_data_sub$Ntype <- ifelse(
  fitted_data_sub$N.coef > 0,
  "Positive",
  "Negative"
)

pred <- fitted_data_sub
names(pred) <- c(
  "species",
  "Nscaled",
  "bioGregion",
  "year",
  "area",
  "clcL2100",
  "p_hat",
  "p_mu",
  "N",
  "N.coef",
  "N.pval",
  "rs",
  "paRatio",
  "bestM",
  "N_ind",
  "N_nic",
  "nitro_group"
)

# Species excluded from the final figure
excluded_species <- c(
  "Carex bigelowii",
  "Galium saxatile",
  "Serratula tinctoria",
  "Senecio nemorensis"
)
pred <- pred[!pred$species %in% excluded_species, ]

# One row per species for panel a
sp <- unique(pred$species)

species_summary <- data.frame(
  species = sp,
  nitro_group = as.character(
    tapply(
      pred$nitro_group,
      pred$species,
      function(z) z[match(TRUE, !is.na(z))]
    )
  ),
  N_ind = as.numeric(
    tapply(pred$N_ind, pred$species, function(z) mean(z, na.rm = TRUE))
  ),
  N_nic = as.numeric(
    tapply(pred$N_nic, pred$species, function(z) mean(z, na.rm = TRUE))
  ),
  slope = as.numeric(
    tapply(pred$N.coef, pred$species, function(z) mean(z, na.rm = TRUE))
  ),
  pval = as.numeric(
    tapply(pred$N.pval, pred$species, function(z) mean(z, na.rm = TRUE))
  ),
  rs = as.numeric(
    tapply(pred$rs, pred$species, function(z) mean(z, na.rm = TRUE))
  ),
  row.names = NULL,
  stringsAsFactors = FALSE
)

species_summary$nitro_group <- factor(species_summary$nitro_group)
species_summary$signif <- ifelse(
  species_summary$pval <= 0.05,
  "<0.05",
  "≥0.05"
)

# Species highlighted in red in panel a and shown in panel b
selected_species <- c(
  "Adonis vernalis",
  "Stipa pennata",
  "Veronica spicata",
  "Geum montanum",
  "Trifolium ochroleucon",
  "Alliaria petiolata",
  "Rorippa palustris",
  "Rumex obtusifolius",
  "Calystegia sepium",
  "Rubus caesius"
)

highlighted_species <- species_summary[
  species_summary$species %in% selected_species,
]

# Panel a: species nitrogen niche vs modelled nitrogen effect
p_species <- ggplot(
  species_summary,
  aes(x = N_ind, y = slope, colour = nitro_group)
) +
  geom_hline(
    yintercept = 0,
    linetype = 2,
    linewidth = 0.2,
    colour = "grey"
  ) +
  ggrepel::geom_text_repel(
    data = species_summary[
      species_summary$signif == "<0.05" &
        species_summary$nitro_group == "Positive",
    ],
    aes(label = species),
    max.overlaps = 10,
    xlim = c(-Inf, 10),
    ylim = c(-Inf, Inf),
    max.time = 20,
    min.segment.length = 0,
    segment.curvature = -0.1,
    segment.ncp = 3,
    segment.angle = 20,
    nudge_x = 0.15,
    box.padding = 0.45,
    nudge_y = 1.2,
    seed = 1987,
    show.legend = FALSE
  ) +
  ggrepel::geom_text_repel(
    data = species_summary[
      species_summary$signif == "<0.05" &
        species_summary$nitro_group == "Negative",
    ],
    aes(label = species),
    max.overlaps = 5,
    min.segment.length = 0,
    segment.curvature = -0.1,
    segment.ncp = 3,
    segment.angle = 20,
    nudge_x = 0.15,
    box.padding = 0.7,
    nudge_y = -1,
    seed = 1987,
    show.legend = FALSE
  ) +
  stat_smooth(
    data = species_summary[species_summary$signif == "<0.05", ],
    aes(x = N_ind, y = slope),
    method = "lm",
    colour = "black",
    se = TRUE,
    linewidth = 1,
    linetype = 1,
    alpha = 0.5,
    show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  geom_point(
    aes(
      shape = interaction(signif, nitro_group),
      size = signif
    ),
    alpha = 0.75
  ) +
  geom_point(
    data = highlighted_species,
    aes(
      shape = interaction(signif, nitro_group),
      size = signif
    ),
    colour = "red",
    alpha = 0.75,
    show.legend = FALSE
  ) +
  scale_shape_manual(
    values = c(
      "<0.05.Negative" = 25,
      "≥0.05.Negative" = 6,
      "<0.05.Positive" = 17,
      "≥0.05.Positive" = 2
    ),
    name = NULL
  ) +
  scale_size_manual(
    values = c(`≥0.05` = 1, `<0.05` = 3),
    name = "p-value"
  ) +
  scale_y_continuous(limits = c(-15, 5)) +
  scale_x_continuous(limits = c(0, 10)) +
  scale_colour_manual(values = c("#5a9eb4", "#d8b365")) +
  labs(
    x = "Nitrogen niche indicator value (EIVE)",
    y = "Nitrogen effect on occurrence (log odds ratio)",
    colour = "Direction of nitrogen\neffect on occurrence",
    tag = "a"
  ) +
  guides(
    colour = guide_legend(
      ncol = 1,
      nrow = 2,
      title.hjust = 1.0,
      override.aes = list(
        shape = c(24, 25),
        colour = c("#5a9eb4", "#d8b365"),
        fill = c("#5a9eb4", "transparent"),
        size = 3
      )
    ),
    size = guide_legend(
      ncol = 1,
      nrow = 2,
      override.aes = list(
        shape = c(24, 25),
        colour = "grey",
        fill = "grey"
      )
    ),
    shape = "none"
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.key.spacing.x = unit(0.5, "cm"),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    legend.title = element_text(
      size = 16,
      face = "bold",
      margin = margin(r = 10)
    ),
    legend.text = element_text(size = 16),
    plot.tag = element_text(size = 21),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Linear regression reported in the Figure 5 legend
lm1 <- lm(
  slope ~ N_ind,
  data = species_summary[species_summary$signif == "<0.05", ]
)
summary(lm1)
confint(lm1)

# Panel b: predicted occurrence curves for selected species
selected_curves <- pred[pred$species %in% selected_species, ]
selected_curves <- selected_curves[
  order(selected_curves$N, selected_curves$species),
]

# Rescale nitrophobic curves to display both groups with separate y-axis labels
max_p <- tapply(
  selected_curves$p_hat,
  selected_curves$nitro_group,
  max,
  na.rm = TRUE
)
scale_factor <- unname(max_p["Positive"] / max_p["Negative"])

selected_curves$nitro_group <- factor(
  selected_curves$nitro_group,
  levels = c("Positive", "Negative"),
  ordered = TRUE
)

y_max <- unname(max_p["Positive"]) + 0.01
y_breaks <- seq(0, y_max, 0.05)

p_curves <- ggplot(
  selected_curves,
  aes(x = N, y = p_hat, colour = nitro_group, group = species)
) +
  geom_line(
    aes(
      y = ifelse(
        nitro_group == "Positive",
        p_hat,
        p_hat * scale_factor
      )
    )
  ) +
  scale_y_continuous(
    limits = c(0, y_max),
    breaks = y_breaks,
    labels = label_scientific(
      round(y_breaks / scale_factor, 4)
    ),
    sec.axis = sec_axis(
      ~ .,
      name = "Predicted probability of presence (nitrophilous)"
    )
  ) +
  geomtextpath::geom_labelpath(
    show.legend = FALSE,
    aes(
      y = ifelse(
        nitro_group == "Positive",
        p_hat,
        p_hat * scale_factor
      ),
      label = paste0(species, " (", round(N_ind, 1), ")"),
      hjust = ifelse(nitro_group == "Positive", 0.99, 0.01)
    ),
    vjust = -0.1,
    text_smoothing = 30,
    size = 4,
    rich = TRUE,
    gap = TRUE,
    alpha = 0.9,
    fontface = "bold.italic",
    parse = FALSE,
    label.padding = grid::unit(1.1, "pt")
  ) +
  scale_colour_manual(values = c("#5a9eb4", "#d8b365")) +
  guides(
    colour = guide_legend(
      ncol = 1,
      nrow = 2,
      title.hjust = 1.0,
      override.aes = list(linewidth = 1.5)
    )
  ) +
  theme_minimal() +
  labs(
    x = "Total nitrogen input (kg ha⁻¹ yr⁻¹)",
    y = "Predicted probability of presence (nitrophobic)",
    colour = "Direction of nitrogen\neffect on occurrence",
    tag = "b"
  ) +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.title = element_text(
      size = 16,
      face = "bold",
      margin = margin(r = 10)
    ),
    legend.text = element_text(size = 16),
    plot.tag = element_text(size = 21),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "top",
    legend.key.spacing.x = unit(0.5, "cm")
  )

figure_05 <- gridExtra::arrangeGrob(
  p_species,
  p_curves,
  ncol = 2,
  widths = c(1.2, 0.8)
)

ggsave(
  file.path(fig_dir, "Figure_05.pdf"),
  plot = figure_05,
  width = 21,
  height = 10,
  units = "cm",
  device = cairo_pdf,
  scale = 2.2
)
