library(TDA)

# ---- 1. Names of methods ----
method_names <- c("Gaussian", "Gaussian_half", "SVM", "TPD", "MASS")
PH_methods <- all_PH[2:6]  # assume all_PH[[1]] is Raw

# ---- 2. Convert Raw PH to numeric matrix ----
PH_raw_df <- as.matrix(
  data.frame(
    dimension = as.integer(all_PH[[1]][, "dimension"]),
    Birth = as.numeric(all_PH[[1]][, "birth"]),
    Death = as.numeric(all_PH[[1]][, "death"])
  )
)

# ---- 3. Convert method PHs to numeric matrices ----
PH_methods_df <- lapply(PH_methods, function(ph) {
  if(is.null(ph)) return(NULL)
  as.matrix(
    data.frame(
      dimension = as.integer(ph[, "dimension"]),
      Birth = as.numeric(ph[, "birth"]),
      Death = as.numeric(ph[, "death"])
    )
  )
})
names(PH_methods_df) <- method_names

# ---- 4. Compute bottleneck distances for H0 and H1 ----
bottleneck_results <- data.frame(
  method = character(),
  dimension = integer(),
  distance = numeric()
)

for(method in method_names){
  ph_method <- PH_methods_df[[method]]
  if(is.null(ph_method)) next  # skip if PH failed
  
  for(dim in 0:1){
    dist <- bottleneck(PH_raw_df, ph_method, dimension = dim)
    bottleneck_results <- rbind(
      bottleneck_results,
      data.frame(method = method, dimension = dim, distance = dist)
    )
  }
}

# ---- 5. Display results ----
bottleneck_results
