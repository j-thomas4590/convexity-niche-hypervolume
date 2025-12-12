# convexity-niche-hypervolume

# Niche hypervolumes are not necessarily convex

## 1. Overview

This repository contains all R scripts and code needed to reproduce the analyses for *“Niche hypervolumes are not necessarily convex”*.

It includes:

* Analyses of **Acacia species occurrences** using PCA, hypervolumes, and TDA.
* Synthetic 2D datasets comparing hypervolume estimation methods (Gaussian KDE, SVM, TPD, MASS KDE).
* Persistent homology (PH) and bottleneck distance computations.
* Publication-quality plots of species niches and TDA results.

---

## 2. Dependencies

All scripts require **R ≥ 4.2** and the following packages:

```r
# Data manipulation
install.packages(c("dplyr", "tidyr", "purrr"))

# Spatial / raster
install.packages(c("raster", "terra"))

# Hypervolume & KDE
install.packages("hypervolume")
install.packages("MASS")
install.packages("e1071")  # for SVM

# TDA
install.packages("TDA")
install.packages("TDAstats")

# Visualization
install.packages(c("ggplot2", "cowplot", "patchwork", "viridis", "geometry"))

# Trait probability densities
install.packages("TPD")  # if not available on CRAN, use devtools::install_github()
```

---

## 3. Acacia occurrence data

* Species occurrence data is provided by the `acacia_pinus` dataset in the `hypervolume` package.
* Environmental layers are downloaded from **WorldClim v2.1** (2.5 arc-minutes): [WorldClim v2.1](https://www.worldclim.org/data/worldclim21.html)

**Required layers:**

* BIO1–BIO19 (temperature & precipitation) in `.tif` format
* Elevation (`wc2.1_2.5m_elev.tif`)

## 4. Script descriptions

### 4.1 `Tda_convex_test_Acacia.R` 

* **Purpose:** Computes global PCA on Acacia occurrence points, builds **hypervolumes per species**, and performs **persistent homology** to test for convexity.
* **Outputs:**

  * `tda_combined`: dataframe with TDA convexity classification per species
  * Summary bar plot comparing convex vs non-convex species
  * PCA scores and thinned hypervolumes for downstream plotting

### 4.2 `Acacia_full_plots.R` 

* **Purpose:** Batch plotting of **occurrence points**, **hypervolumes**, and **persistent homology diagrams** for all Acacia species.
* **Key functions:**

  * `make_scatter_plot()`: plots PCA scores or hypervolume points
  * `make_PH_plot()`: plots persistent homology (birth–death) diagrams
  * `make_page()`: combines species plots into pages for visualization
* **Outputs:** PNG pages in `pages/` folder.

### 4.3 `Comparing_hypervolumes_synthetic_data.R` 

* **Purpose:** Creates **synthetic 2D datasets** of two circular blobs at varying distances.
* **Generates hypervolume estimates** using:

  * Gaussian KDE (`hypervolume_gaussian`)
  * Half-bandwidth Gaussian KDE
  * SVM-based hypervolume (`hypervolume_svm`)
  * MASS KDE (`kde2d`)
  * TPD
* **Outputs:** `all_points_df` containing all point clouds.

### 4.4 `Tda_touching_synthetic_data.R`

* **Purpose:** Runs **TDA (persistent homology)** on the synthetic datasets produced in Code 4.
* **Tasks:**

  * Thins points for plotting (≤200 points)
  * Computes PH (H0 + H1)
  * Generates scatter + PH plots in a grid for each method
* **Outputs:** Final combined grid showing raw points, hypervolume points, and PH diagrams.

### 4.5 `bottleneck_synthetic.R` 

* **Purpose:** Computes **bottleneck distances** between raw PH and PH from different hypervolume methods for synthetic datasets.
* **Outputs:** `bottleneck_results` dataframe with distances per method and dimension (H0 and H1).

---



---



