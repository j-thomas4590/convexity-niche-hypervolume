# --------------------------------------------------------------
# Global PCA + Hypervolume + TDA convexity analysis (corrected)
# --------------------------------------------------------------

library(raster)
library(terra)
library(dplyr)
library(hypervolume)
library(geometry)
library(TDAstats)
library(ggplot2)
library(tidyr)

set.seed(4000)

# --- 0. Load data ---
data(acacia_pinus)
acacia_occ <- acacia_pinus[grep("^Acacia", acacia_pinus$species), ]

# --- Load environmental layers ---
bio_path <- "data/worldclim/bio/"
bio_files <- list.files(bio_path, pattern = "\\.tif$", full.names = TRUE)
bio_stack <- rast(bio_files)

# elevation raster
elev_path <- "data/worldclim/elev/elevation.tif"
elev <- rast(elev_path)

# combine
env_stack <- c(bio_stack, elev)

# rename layer
names(env_stack)[20] <- "elevation"


# --- Extract environmental values at occurrence points ---
coords <- acacia_occ[, c("longitude","latitude")]
env_vals <- terra::extract(env_stack, coords)
acacia_env_data <- cbind(acacia_occ["species"], env_vals)
acacia_env_data <- na.omit(acacia_env_data)

# --- Standardize environmental variables (once) ---
acacia_scaled <- acacia_env_data %>% mutate(across(-species, scale))

# --- GLOBAL PCA across all species ---
env_points_all <- acacia_scaled[, -1]
pca_global <- prcomp(env_points_all, center = TRUE, scale. = FALSE)  # data already scaled

# --- Variance explained ---
var_explained <- summary(pca_global)$importance[2, 1:3]
cum_var_explained <- summary(pca_global)$importance[3, 3] * 100
cat("Variance explained by PC1–3:\n")
print(var_explained)
cat("Cumulative variance explained:", round(cum_var_explained, 2), "%\n\n")

# --- Project all points into PCA space ---
pca_scores_all <- as.data.frame(pca_global$x[, 1:3])
pca_scores_all$species <- acacia_scaled$species

# --- Split by species ---
species_list <- split(pca_scores_all, pca_scores_all$species)

# --- Output containers ---
hv_list <- list()
pca_scores_list <- list()
PH_list_hv <- list()     # NEW: raw PH results for HV
PH_list_occ <- list()    # NEW: raw PH results for PCA
n_thin <- 10000

# ------------------------------------------------------------------
# 1. Build hypervolumes and store PCA scores
# ------------------------------------------------------------------
for (sp in names(species_list)) {
  cat("Processing species:", sp, "\n")
  env_points <- species_list[[sp]][, 1:3]

  if (nrow(env_points) > 5) {
    hv <- hypervolume_gaussian(env_points, samples.per.point = 50)
    hv_thinned <- hypervolume_thin(hv, num.points = n_thin)
    hv_list[[sp]] <- hv_thinned

    pca_scores_list[[sp]] <- cbind(env_points, species = sp)
  } else {
    cat("Skipping", sp, "- not enough points\n")
  }
}

# ------------------------------------------------------------------
# 2. Helper function for TDA significance thresholds
# ------------------------------------------------------------------
id_significant_adapted <- function(features, dim = 1, reps = 100, cutoff = 0.975) {
  colnames(features) <- c("dimension","birth","death")
  features <- features[features[,1] == dim,]
  if (nrow(features) < 1) return(0)

  features$persist <- features$death - features$birth
  ans <- replicate(reps, mean(sample(features$persist, size = nrow(features), replace = TRUE)))
  stats::quantile(ans, cutoff, names = FALSE)
}

# ------------------------------------------------------------------
# 3. TDA on hypervolume points
# ------------------------------------------------------------------
tda_hv_results <- data.frame(species = character(), convex_tda_hv = logical())

for (sp in names(hv_list)) {
  hv <- hv_list[[sp]]
  pts <- hv@RandomPoints
  if (is.null(pts) || nrow(pts) < 5) next

  if (nrow(pts) > 300) pts <- pts[sample(nrow(pts), 300), ]

  ph <- calculate_homology(pts, dim=2, format="cloud")
  PH_list_hv[[sp]] <- ph          # SAVE FOR CODE 2

  ph_df <- as.data.frame(ph)

  H1_thresh <- id_significant_adapted(ph_df, dim=1, reps=500)
  H2_thresh <- id_significant_adapted(ph_df, dim=2, reps=500)

  has_H1 <- any((ph_df$dimension==1) & ((ph_df$death - ph_df$birth) > H1_thresh))
  has_H2 <- any((ph_df$dimension==2) & ((ph_df$death - ph_df$birth) > H2_thresh))

  tda_hv_results <- rbind(tda_hv_results,
                          data.frame(species=sp,
                                     convex_tda_hv = !(has_H1 | has_H2)))
}

# ------------------------------------------------------------------
# 4. TDA on PCA occurrence points
# ------------------------------------------------------------------
tda_pca_results <- data.frame(species = character(), convex_tda_pca = logical())

for (sp in names(pca_scores_list)) {
  pts <- as.matrix(pca_scores_list[[sp]][, 1:3])
  storage.mode(pts) <- "double"
  if (nrow(pts) < 5) next
  if (nrow(pts) > 300) pts <- pts[sample(nrow(pts), 300), ]

  ph <- calculate_homology(pts, dim=2, format="cloud")
  PH_list_occ[[sp]] <- ph         # SAVE FOR CODE 2

  ph_df <- as.data.frame(ph)

  H1_thresh <- id_significant_adapted(ph_df, dim=1, reps=500)
  H2_thresh <- id_significant_adapted(ph_df, dim=2, reps=500)

  has_H1 <- any((ph_df$dimension==1) & ((ph_df$death - ph_df$birth) > H1_thresh))
  has_H2 <- any((ph_df$dimension==2) & ((ph_df$death - ph_df$birth) > H2_thresh))

  tda_pca_results <- rbind(tda_pca_results,
                           data.frame(species=sp,
                                      convex_tda_pca = !(has_H1 | has_H2)))
}

# ------------------------------------------------------------------
# 5. Combine results
# ------------------------------------------------------------------
tda_combined <- merge(tda_hv_results, tda_pca_results, by="species", all=TRUE)

# ------------------------------------------------------------------
# 6. Summary table
# ------------------------------------------------------------------
summary_counts <- data.frame(
  Test = c("TDA_Hypervolume","TDA_PCA"),
  Convex = c(sum(tda_combined$convex_tda_hv, na.rm=TRUE),
             sum(tda_combined$convex_tda_pca, na.rm=TRUE)),
  NonConvex = c(sum(!tda_combined$convex_tda_hv, na.rm=TRUE),
                sum(!tda_combined$convex_tda_pca, na.rm=TRUE))
)

# ------------------------------------------------------------------
# 7. Bar plot summary
# ------------------------------------------------------------------
summary_long <- pivot_longer(summary_counts,
                             cols=c("Convex","NonConvex"),
                             names_to="Type",
                             values_to="Count")

ggplot(summary_long, aes(x=Test, y=Count, fill=Type)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_manual(values=c("Convex"="steelblue","NonConvex"="firebrick")) +
  theme_minimal(base_size=14) +
  labs(title="Convex vs Non-Convex Across Methods",
       x="Method", y="Number of Species", fill="Classification")

# ------------------------------------------------------------------
# 8. Print final results
# ------------------------------------------------------------------
cat("\nCumulative variance explained (PC1–3):", round(cum_var_explained, 2), "%\n")
tda_combined
