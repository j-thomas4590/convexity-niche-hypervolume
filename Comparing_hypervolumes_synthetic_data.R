############################################################
# HYPERVOLUME vs TPD vs MASS KDE vs SVM COMPARISON (POINT CLOUDS)
############################################################

library(hypervolume)
library(TPD)
library(MASS)
library(e1071)
library(ggplot2)
library(dplyr)
library(reshape2)

set.seed(123)

# ---- 1. Helper: generate two circular blobs ----
make_data <- function(n = 200, r = 1, centers = c(-1, 1)) {
  theta1 <- runif(n, 0, 2*pi)
  x1 <- r * cos(theta1) + centers[1] + rnorm(n, 0, 0.1)
  y1 <- r * sin(theta1) + rnorm(n, 0, 0.1)

  theta2 <- runif(n, 0, 2*pi)
  x2 <- r * cos(theta2) + centers[2] + rnorm(n, 0, 0.1)
  y2 <- r * sin(theta2) + rnorm(n, 0, 0.1)

  data.frame(x = c(x1, x2), y = c(y1, y2))
}

# ---- 2. Three datasets ----
data_touching  <- make_data(centers = c(-1, 1))
data_separated <- make_data(centers = c(-3, 3))
data_far       <- make_data(centers = c(-5, 5))

datasets <- list(Touching = data_touching,
                 Separated = data_separated,
                 Far = data_far)

# ---- 3. Hypervolumes ----
compute_hv_points <- function(data) {
  bw <- estimate_bandwidth(data[,1:2])
  
  hv_gauss <- tryCatch(hypervolume_gaussian(data[,1:2], kde.bandwidth = bw), 
                       error = function(e) NULL)
  hv_gauss_half <- tryCatch(hypervolume_gaussian(data[,1:2], kde.bandwidth = bw / 2),
                            error = function(e) NULL)
  
  hv_points <- function(hv, method) {
    if(is.null(hv) || is.null(hv@RandomPoints)) return(data.frame(x=numeric(0), y=numeric(0), method=method))
    as.data.frame(hv@RandomPoints) %>% rename(x=1, y=2) %>% mutate(method=method)
  }
  
  list(
    Gaussian = hv_points(hv_gauss, "Gaussian"),
    Gaussian_half = hv_points(hv_gauss_half, "Gaussian_half")
  )
}

# ---- 4. SVM points ----
compute_svm_points <- function(data, npoints=5000) {
  # Estimate bandwidth using default in hypervolume
  bw <- estimate_bandwidth(data[,1:2])
  
  # Compute hypervolume using SVM
  hv_svm <- tryCatch(
    hypervolume_svm(data[,1:2], 
                    name = "SVM", 
                    samples.per.point = 1000,  # number of random points for internal HV
                    verbose = FALSE),
    error = function(e) NULL
  )
  
  # Extract random points from the hypervolume
  if(!is.null(hv_svm) && !is.null(hv_svm@RandomPoints)) {
    points <- as.data.frame(hv_svm@RandomPoints) %>% 
      rename(x=1, y=2) %>% 
      mutate(method = "SVM")
    
    # Optional: subsample to npoints if too many
    if(nrow(points) > npoints) points <- points[sample(1:nrow(points), npoints),]
    
  } else {
    # fallback if HV fails
    points <- data.frame(x=numeric(0), y=numeric(0), method="SVM")
  }
  
  return(points)
}


# ---- 5. TPD and MASS points with higher threshold ----
compute_mass_points <- function(data, threshold=0.3, nbins=100) {  # higher threshold
  bw <- c(bandwidth.nrd(data$x), bandwidth.nrd(data$y))
  kde_res <- kde2d(data$x, data$y, n=nbins, h=bw)
  z <- kde_res$z
  z_thresh <- max(z) * threshold
  idx <- which(z >= z_thresh, arr.ind=TRUE)
  points <- data.frame(
    x = kde_res$x[idx[,1]],
    y = kde_res$y[idx[,2]],
    method = "MASS"
  )
  return(points)
}

compute_tpd_points <- function(data, threshold=0.3) {  # higher threshold
  traits <- as.matrix(data)
  species <- rep("sp1", nrow(traits))
  TPDsp_result <- TPDs(species = species, traits = traits)
  tpd_vec <- TPDsp_result$TPDs[[1]]
  grid_size <- round(sqrt(length(tpd_vec)))
  trait1_grid <- seq(min(traits[,1]), max(traits[,1]), length.out = grid_size)
  trait2_grid <- seq(min(traits[,2]), max(traits[,2]), length.out = grid_size)
  prob_matrix <- matrix(tpd_vec, nrow = grid_size, ncol = grid_size, byrow = TRUE)
  thresh <- max(prob_matrix) * threshold
  idx <- which(prob_matrix >= thresh, arr.ind=TRUE)
  
  points <- data.frame(
    x = trait1_grid[idx[,2]],  # FIX: column -> trait1
    y = trait2_grid[idx[,1]],  # row -> trait2
    method = "TPD"
  )
  return(points)
}

# ---- 6. Assemble all points ----
all_points <- list()
for(name in names(datasets)) {
  data <- datasets[[name]]
  hv_pts <- compute_hv_points(data)
  svm_pts <- compute_svm_points(data)
  tpd_pts <- compute_tpd_points(data)
  mass_pts <- compute_mass_points(data)
  
  all_points[[name]] <- list(
    Raw = data %>% mutate(method="Raw"),
    Gaussian = hv_pts$Gaussian,
    Gaussian_half = hv_pts$Gaussian_half,
    SVM = svm_pts,
    TPD = tpd_pts,
    MASS = mass_pts
  )
}


# ---- 7. Combine into one data.frame with ordered factors ----
all_points_df <- bind_rows(lapply(names(all_points), function(scenario) {
  bind_rows(lapply(names(all_points[[scenario]]), function(method) {
    df <- all_points[[scenario]][[method]]
    if(nrow(df) > 0) {  # only assign factors if df is non-empty
      df$scenario <- factor(scenario, levels = c("Touching", "Separated", "Far"))
      df$method <- factor(method, levels = c("Raw","Gaussian","Gaussian_half","SVM","TPD","MASS"))
    } else {  # create empty df with correct columns for bind_rows
      df <- data.frame(x = numeric(0), y = numeric(0),
                       scenario = factor(levels = c("Touching","Separated","Far")),
                       method = factor(levels = c("Raw","Gaussian","Gaussian_half","SVM","TPD","MASS")))
    }
    df
  }))
}))


# ---- 8. Dynamic limits ----
xlim <- range(all_points_df$x, na.rm=TRUE) + c(-1,1)
ylim <- range(all_points_df$y, na.rm=TRUE) + c(-1,1)

# ---- 9. Color palette ----
method_colors <- c(
  Raw = "#1b9e77",
  Gaussian = "#d95f02",
  Gaussian_half = "#7570b3",
  SVM = "#e7298a",
  TPD = "#66a61e",
  MASS = "#e6ab02"
)

# ---- 10. Faceted plot ----
ggplot(all_points_df, aes(x=x, y=y, color=method)) +
  geom_point(size=0.2, alpha=0.5) +
  facet_grid(scenario ~ method) +
  scale_color_manual(values=method_colors) +
  coord_equal(xlim = xlim, ylim = ylim) +
  theme_void() +  # removes all axes, labels, ticks
  theme(
    strip.text = element_text(face="bold", size=12),
    legend.position = "none"
  ) +
  labs(title = "Comparison of Hypervolume, TPD, SVM, and KDE Point Clouds")
