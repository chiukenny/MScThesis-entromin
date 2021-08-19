# Code for NYT articles analysis
# https://archive.ics.uci.edu/ml/datasets/bag+of+words
# ---------------------------------------------------------


set.seed(1)

library(data.table)
library(Matrix)
rm(list=ls())

source("functions/dataset.R")
source("functions/orth_rotation.R")
source("functions/util.R")

# Control parameters
load_raw_data = T # Load raw data from saved file
load_data     = T # Load PCs from saved file
save_raw_data = F # Save raw data if downloaded
save_data     = F # Save PCs to saved file
save_plots    = T # Save plots and related data

# Rotation parameters
k_init = 50 # Initial number of PCs 
k      = 8  # Final number of PCs

tau           = -1 # Regularization parameter for bigPCA
normalization = T  # Normalize data before PCA
centering     = T  # Center data before PCA
recenter      = T  # Recenter rotated factors
fast          = F  # Compute rotation based on 10% of rows

max_iters = 1000 # Max number of rotation iterations
verbose   = F    # Print intermediate progress

# Plotting parameters
num_samples = 5000 # Number of samples to plot
p_height    = 6.5  # Plot height
p_width     = 6.5  # Plot width



# Data preprocessing
# ---------------------------------------------------------

# Read raw data
A = get_NYT("data/nytA.RData", from_saved=load_raw_data, save_data=save_raw_data)

# Get PCs
if (load_data)
{
  load(file="outputs/data/pcs.RData")
} else {
  pcs = bigPCA(A, k_init, tau, normalization, centering)
  if (save_data)
  {
    save(pcs, file="outputs/data/pcs.RData")
  }
}

# Check if PCs are localized on a few documents following Rohe (2020)
U = pcs$U[,1:k_init]
L4 = colSums(U^4)^(1/4)
if (save_plots)
{
  # Save L4 norms
  df_loc = data.frame(k=1:k_init, L4=L4)
  write.table(df_loc, "localized.txt", quote=F, row.names=F)
}

# Remove PCs with L4 norm greater than 0.15
localized = L4 > 0.15
if (save_plots)
{
  # Save data for scree plot
  unloc_ind = which(!localized)
  df_scree = data.frame(k=1:length(unloc_ind), scree=pcs$scree[unloc_ind])
  write.table(df_scree, "scree.txt", quote=F, row.names=F)
}

# Keep leading k=8 PCs (eigengap at 8)
nonloc_cols = which(!localized)[1:8]
kpcs = keep_factors(pcs, nonloc_cols)

# Select samples to plot
if (save_plots)
{
  p = sqrt(rowSums(U[,nonloc_cols]^2))
  m = min(num_samples, sum(p>0))
  samples = sample(nrow(U), m, replace=F, p)
}



# Factor rotation methods
# ---------------------------------------------------------

# Varimax
system.time(
  var_factors <- varimax_rotate(kpcs, recenter=recenter, fast=fast,
                                max_iters=max_iters, verbose=verbose)
)
get_stats(kpcs$U %*% var_factors$Rz)

# Entromin2
system.time(
 entr2_factors <- entromin2_rotate(kpcs, recenter=recenter, fast=fast,
                                   max_iters=max_iters, verbose=verbose)
)
get_stats(kpcs$U %*% entr2_factors$Rz)

# Entromin
system.time(
  entr_factors <- entromin_rotate(kpcs, recenter=recenter, fast=fast,
                                  max_iters=max_iters, verbose=verbose)
)
get_stats(kpcs$U %*% entr_factors$Rz)



# Plot results
# ---------------------------------------------------------

if (save_plots)
{
  # Unrotated vs Varimax rotated
  png(file="outputs/plots/nyt_raw_rotated.png",
      height=p_height, width=p_width, units="in", res=200)
  plot_rot_unrot(U[,nonloc_cols], var_factors$Z, k, samples)
  dev.off()
  
  
  # Entromin vs Varimax
  entr_col = c(6,2,5,3,8,4,7,1) # Match Varimax PCs manually
  P_entr = diag(k)[,entr_col]
  png(file="outputs/plots/nyt_entr_var.png",
      height=p_height, width=p_width, units="in", res=200)
  plot_diff_pairs(var_factors$Z %*% var_factors$P,
                  entr_factors$Z %*% P_entr,
                  k, samples)
  dev.off()
  
  
  # Entromin vs Varimax
  entr2_col = c(7,6,5,3,8,4,2,1) # Match Varimax PCs manually
  P_entr2 = diag(k)[,entr2_col]
  png(file="outputs/plots/nyt_entr2_var.png",
      height=p_height, width=p_width, units="in", res=200)
  plot_diff_pairs(var_factors$Z %*% var_factors$P,
                  entr2_factors$Z %*% P_entr2,
                  k, samples)
  dev.off()
}