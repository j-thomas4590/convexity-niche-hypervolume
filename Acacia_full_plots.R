# --------------------------------------------------------------
# Batch plotting of niche + PH diagrams for all species
# --------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)
library(viridis)
library(cowplot)
library(patchwork)

# --------------------------------------------------------------
# Helper functions
# --------------------------------------------------------------

make_PH_plot <- function(ph_data, global_min, global_max) {
  if(is.null(ph_data) || nrow(as.data.frame(ph_data))==0) return(NULL)
  
  max_dim <- 2
  # Threshold calculation
  get_sign_dims <- function(hom_data, max_dim = 2, reps = 500, cutoff = 0.975) {
    id_significant_adapted <- function(features, dim = 1, reps = 100, cutoff = 0.975) {
      colnames(features) <- c("dimension", "birth", "death")
      features <- features[features[, 1] == dim, ]
      if (nrow(features) < 1) return(NA)
      features$persist <- features$death - features$birth
      ans <- replicate(reps, mean(sample(features$persist, size = nrow(features), replace = TRUE)))
      stats::quantile(ans, cutoff, names = FALSE)
    }
    thresholds <- sapply(0:max_dim, function(d)
      id_significant_adapted(as.data.frame(hom_data), dim = d, reps = reps, cutoff = cutoff))
    data.frame(threshold = thresholds, dimension = 0:max_dim)
  }
  
  assign_cols_PH <- function(data){
    data <- data %>%
      mutate(temp = paste(dimension, colour)) %>%
      mutate(colvar = case_when(
        temp == "0 no" ~ "grey60",
        temp == "1 no" ~ "grey60",
        temp == "2 no" ~ "grey60",
        temp == "0 sig" ~ "firebrick2",
        temp == "1 sig" ~ "steelblue2",
        temp == "2 sig" ~ "orchid3"
      ))
    return(data)
  }
  
  ph_sig <- get_sign_dims(ph_data, max_dim=max_dim)
  ph_df <- as.data.frame(ph_data) %>%
    mutate(colour = case_when(
      dimension==2 ~ "sig",
      dimension %in% c(0,1) & (death-birth) > ph_sig$threshold[dimension+1] ~ "sig",
      TRUE ~ "no"
    )) %>%
    assign_cols_PH()
  
  ggplot(ph_df, aes(x=birth, y=death)) +
    geom_point(aes(shape=as.factor(dimension), color=colvar), size=3) +
    geom_abline(slope=1, intercept=0, linetype="dashed", color="black") +
    scale_color_identity() +
    scale_shape_manual(values=c(15,17,19)) +
    coord_cartesian(xlim=c(global_min, global_max), ylim=c(global_min, global_max)) +
    theme_linedraw(base_size=14) +
    theme(panel.grid=element_blank(),
          legend.position="none",
          axis.title=element_text(face="bold", size=12),
          axis.text=element_text(size=10))
}

make_scatter_plot <- function(pca_df, global_PC3=NULL){
  if(is.null(pca_df) || nrow(pca_df)==0) return(NULL)
  
  p <- ggplot(pca_df, aes(x=PC1, y=PC2, color=PC3)) +
    geom_point(alpha=0.8, size=3) +
    theme_minimal(base_size=14) +
    theme(axis.title=element_text(face="bold", size=12),
          axis.text=element_text(size=10),
          legend.position="none")
  
  if(!is.null(global_PC3)){
    p <- p + scale_color_viridis_c(option="C", limits=global_PC3)
  } else {
    p <- p + scale_color_viridis_c(option="C")
  }
  p
}

# --------------------------------------------------------------
# Single page function
# --------------------------------------------------------------
make_page <- function(species_to_plot, pca_scores_list, hv_list, PH_list_occ, PH_list_hv, global_PC3=NULL) {
  
  # Global birth/death range for PH plots
  all_births <- c()
  all_deaths <- c()
  for(sp in species_to_plot){
    if(!is.null(PH_list_occ[[sp]])) {
      ph_df <- as.data.frame(PH_list_occ[[sp]])
      all_births <- c(all_births, ph_df$birth)
      all_deaths <- c(all_deaths, ph_df$death)
    }
    if(!is.null(PH_list_hv[[sp]])) {
      ph_df <- as.data.frame(PH_list_hv[[sp]])
      all_births <- c(all_births, ph_df$birth)
      all_deaths <- c(all_deaths, ph_df$death)
    }
  }
  buffer <- 0.05 * diff(range(c(all_births, all_deaths)))
  global_min <- min(c(all_births, all_deaths), na.rm=TRUE) - buffer
  global_max <- max(c(all_births, all_deaths), na.rm=TRUE) + buffer
  
  # Species rows
  species_rows <- lapply(species_to_plot, function(sp){
    occ_scatter <- make_scatter_plot(pca_scores_list[[sp]][, c("PC1","PC2","PC3")], global_PC3)
    occ_PH <- make_PH_plot(PH_list_occ[[sp]], global_min, global_max)
    
    hv_scatter <- NULL
    if(!is.null(hv_list[[sp]])){
      hv_scatter <- make_scatter_plot(as.data.frame(hv_list[[sp]]@RandomPoints)[,1:3] %>%
                                        setNames(c("PC1","PC2","PC3")), global_PC3)
    }
    hv_PH <- make_PH_plot(PH_list_hv[[sp]], global_min, global_max)
    
    row_plots <- list(occ_scatter, occ_PH, hv_scatter, hv_PH) |> purrr::compact()
    
    widths <- sapply(row_plots, function(p){
      if(!is.null(p) && "birth" %in% names(ggplot_build(p)$data[[1]])) return(1.5) else return(1)
    })
    
    wrap_plots(row_plots, ncol=length(row_plots), widths=widths)
  })
  
  # Species labels
  species_rows_labeled <- mapply(function(row_plot, sp){
    plot_grid(
      ggplot() + annotate("text", x=0.5, y=0.5, label=sp, angle=90, size=6) + theme_void(),
      row_plot,
      ncol=2, rel_widths=c(0.08,1)
    )
  }, species_rows, species_to_plot, SIMPLIFY=FALSE)
  
  # Column headers
  col_header <- plot_grid(
    ggplot() + annotate("text", x=0.5, y=0.5, label="Occurrence Data", size=6) + theme_void(),
    ggplot() + theme_void(),
    ggplot() + annotate("text", x=0.5, y=0.5, label="Hypervolume", size=6) + theme_void(),
    ggplot() + theme_void(),
    ncol=4
  )
  
  main_grid <- wrap_plots(species_rows_labeled, ncol=1)
  
  # Legend (from first scatter)
  legend_plot <- get_legend(
    ggplot(pca_scores_list[[species_to_plot[1]]][, c("PC1","PC2","PC3")],
           aes(x=PC1, y=PC2, color=PC3)) +
      geom_point() +
      scale_color_viridis_c(option="C", name="PC3", limits=global_PC3) +
      theme_minimal(base_size=14) +
      theme(legend.position="right")
  )
  
  # Combine grid + legend
  page_grid <- plot_grid(
    plot_grid(col_header, main_grid, ncol=1, rel_heights=c(0.05,1)),
    legend_plot,
    ncol=2, rel_widths=c(4,0.4)
  )
  
  page_grid
}

# --------------------------------------------------------------
# Prepare batches of 6 species
# --------------------------------------------------------------
species_all <- sort(names(pca_scores_list))
species_batches <- split(species_all, ceiling(seq_along(species_all)/6))

# Compute global PC3 range across all species for consistent coloring
all_PC3 <- unlist(lapply(pca_scores_list, function(df) df$PC3))
global_PC3 <- range(all_PC3, na.rm=TRUE)

# Create pages folder
dir.create("pages", showWarnings = FALSE)
# Optionally delete old files
old_files <- list.files("pages", pattern="\\.png$", full.names=TRUE)
file.remove(old_files)

# --------------------------------------------------------------
# Generate and save all 12 pages
# --------------------------------------------------------------
for(i in seq_along(species_batches)){
  sp_list <- species_batches[[i]]
  cat("Rendering page", i, "with species:", paste(sp_list, collapse=", "), "\n")
  page_plot <- make_page(sp_list, pca_scores_list, hv_list, PH_list_occ, PH_list_hv, global_PC3)
  
  ggsave(
    filename = sprintf("pages/niche_page_%02d.png", i),
    plot = page_plot,
    width = 12, height = 16, dpi = 300
  )
}
