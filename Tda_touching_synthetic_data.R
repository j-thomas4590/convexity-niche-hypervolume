library(dplyr)
library(ggplot2)
library(TDAstats)
library(patchwork)
library(cowplot)
library(viridis)
library(purrr)

set.seed(123)

# ---- 0. Hypervolume colors ----
method_colors <- c(
  Raw = "#1b9e77",
  Gaussian = "#d95f02",
  Gaussian_half = "#7570b3",
  SVM = "#e7298a",
  TPD = "#66a61e",
  MASS = "#e6ab02"
)

# ---- 1. Thin points to 200 ----
thin_points <- function(df, n=200){
  if(nrow(df) > n) df <- df[sample(1:nrow(df), n), ]
  df
}

# ---- 2. PH calculation ----
calc_PH <- function(df){
  if(nrow(df) < 2) return(NULL)
  calculate_homology(as.matrix(df[, c("x","y")]), dim=1)  # only need up to dim=1 for 2D
}

# ---- 3. PH plotting function (H0 + H1 only) ----
make_PH_plot <- function(ph_data, global_min, global_max){
  if(is.null(ph_data) || nrow(as.data.frame(ph_data))==0) return(NULL)
  
  ph_df <- as.data.frame(ph_data) %>%
    filter(dimension %in% c(0,1)) %>%  # drop H2
    mutate(colour = case_when(
      dimension==0 ~ "firebrick2",
      dimension==1 ~ "steelblue2"
    )) %>%
    mutate(shape_val = ifelse(dimension==0, 15, 17))  # H0=15 square, H1=17 triangle
  
  ggplot(ph_df, aes(x=birth, y=death)) +
    geom_point(aes(shape=shape_val, color=colour), size=3) +
    geom_abline(slope=1, intercept=0, linetype="dashed", color="black") +
    scale_color_identity() +
    scale_shape_identity() +
    coord_cartesian(xlim=c(global_min, global_max), ylim=c(global_min, global_max)) +
    labs(x="Birth", y="Death") +
    theme_linedraw(base_size=14) +
    theme(panel.grid=element_blank(),
          legend.position="none",
          axis.title=element_text(face="bold", size=12),
          axis.text=element_text(size=10))
}

# ---- 4. Scatter plot function ----
make_scatter_plot <- function(df, title, color="#0072B2"){
  ggplot(df, aes(x=x, y=y)) +
    geom_point(color=color, alpha=0.6, size=2.5) +
    labs(title=title, x="X", y="Y") +
    theme_minimal(base_size=14) +
    theme(axis.title=element_text(face="bold", size=12),
          axis.text=element_text(size=10),
          plot.title=element_text(face="bold", size=14, hjust=0.5))
}

# ---- 5. Select the point clouds to plot (thinned) ----
df_raw <- thin_points(all_points$Touching$Raw)
df_hv1 <- thin_points(all_points$Touching$Gaussian)
df_hv2 <- thin_points(all_points$Touching$Gaussian_half)
df_hv3 <- thin_points(all_points$Touching$SVM)
df_hv4 <- thin_points(all_points$Touching$TPD)
df_hv5 <- thin_points(all_points$Touching$MASS)

all_dfs <- list(df_raw, df_hv1, df_hv2, df_hv3, df_hv4, df_hv5)

# ---- 6. Compute PH for all (safe for NULL) ----
all_PH <- lapply(all_dfs, calc_PH)

# ---- 7. Determine global PH limits (fixed for matrices) ----
ph_nonnull <- all_PH[!sapply(all_PH, is.null)]  # filter out NULLs

# convert each PH matrix to data.frame before extracting birth/death
all_births <- unlist(lapply(ph_nonnull, function(x) as.data.frame(x)$birth))
all_deaths <- unlist(lapply(ph_nonnull, function(x) as.data.frame(x)$death))

buffer <- 0.05 * diff(range(c(all_births, all_deaths)))
global_min <- min(c(all_births, all_deaths)) - buffer
global_max <- max(c(all_births, all_deaths)) + buffer

# ---- 8. Generate scatter plots with method colors ----
scatter_plots <- list(
  make_scatter_plot(df_raw, "Raw Data", color=method_colors["Raw"]),
  make_scatter_plot(df_hv1, "Gaussian HV", color=method_colors["Gaussian"]),
  make_scatter_plot(df_hv2, "Gaussian Half HV", color=method_colors["Gaussian_half"]),
  make_scatter_plot(df_hv3, "SVM HV", color=method_colors["SVM"]),
  make_scatter_plot(df_hv4, "TPD HV", color=method_colors["TPD"]),
  make_scatter_plot(df_hv5, "MASS HV", color=method_colors["MASS"])
)

# ---- 9. Generate PH plots ----
ph_plots <- map(all_PH, ~ make_PH_plot(.x, global_min, global_max))

# ---- 10. Combine into 3x4 grid (row-wise) ----
final_grid <- plot_grid(
  scatter_plots[[1]], ph_plots[[1]], scatter_plots[[2]], ph_plots[[2]],
  scatter_plots[[3]], ph_plots[[3]], scatter_plots[[4]], ph_plots[[4]],
  scatter_plots[[5]], ph_plots[[5]], scatter_plots[[6]], ph_plots[[6]],
  ncol=4,
  align='hv'
)

# ---- 11. Display final figure ----
final_grid
