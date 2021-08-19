# Code for visual sparsity analysis on the simulated
# sparse Gaussian dataset
# ---------------------------------------------------------


set.seed(1)

library(latex2exp)
rm(list=ls())

source("functions/dataset.R")
source("functions/orth_rotation.R")
source("functions/util.R")

# Control parameters
save_plots = T # Save plots
first_row  = T # Include method names as titles

# Dataset parameters
p0 = 0.3  # Probability
n  = 1000 # Number of samples
d  = 3    # Number of columns

# Rotation parameters
k = 2 # Final number of PCs

tau           = -1 # Regularization parameter for bigPCA
normalization = F  # Normalize data before PCA
centering     = F  # Center data before PCA
recenter      = F  # Recenter rotated factors
fast          = F  # Compute rotation based on 10% of rows

max_iters = 1000 # Max number of rotation iterations
verbose   = F    # Print intermediate progress

# Plotting parameters
p_height    = 1.45 # Plot height
p_width     = 5.5  # Plot width

# File parameter
plot_dest = "outputs/plots/"



# Data preprocessing
# ---------------------------------------------------------

A = get_sparse_normal(p0, n=n, k=d)
pcs = bigPCA(A, k, tau, normalization, centering)



# Factor rotation methods
# ---------------------------------------------------------

# Varimax
var_factors = varimax_rotate(pcs, recenter=recenter, fast=fast,
                             max_iters=max_iters, verbose=verbose)
get_stats(pcs$U %*% var_factors$Rz, 1e-4)

# Entromin2
entr2_factors = entromin2_rotate(pcs, recenter=recenter, fast=fast,
                                 max_iters=max_iters, verbose=verbose)
get_stats(pcs$U %*% entr2_factors$Rz, 1e-4)

# Entromin
entr_factors = entromin_rotate(pcs, recenter=recenter, fast=fast,
                               max_iters=max_iters, verbose=verbose)
get_stats(pcs$U %*% entr_factors$Rz, 1e-4)



# Plot results
# ---------------------------------------------------------

if (save_plots)
{
  p0s = format(p0, digits=2)
  if (first_row)
  {
    # Adjust height for first row
    png(file=paste(plot_dest,p0s,"norm_comparison.png",sep=""),
        height=p_height+0.35, width=p_width, units="in", res=200)
    par(mfrow=c(1,4), mai=c(0.05,0.05,0.4,0.05))
    
    Ulim = max(abs(range(pcs$U)))
    Zlim = max(abs(range(rbind(var_factors$Z, entr_factors$Z, entr2_factors$Z))))
    
    plot_pair(rotate_random(pcs$U), xlim=c(-Ulim,Ulim), ylim=c(-Ulim,Ulim),
              main=TeX("Random"), cex.main=1.9)
    plot_pair(var_factors$Z, xlim=c(-Zlim,Zlim), ylim=c(-Zlim,Zlim),
              main=TeX("Varimax"), cex.main=1.9)
    plot_pair(entr2_factors$Z, xlim=c(-Zlim,Zlim), ylim=c(-Zlim,Zlim),
              main=TeX("Entromin$_2$"), cex.main=1.9)
    plot_pair(entr_factors$Z, xlim=c(-Zlim,Zlim), ylim=c(-Zlim,Zlim),
              main=TeX("Entromin"), cex.main=1.9)
    dev.off()
  } else {
    png(file=paste(plot_dest,p0s,"norm_comparison.png",sep=""),
        height=p_height, width=p_width, units="in", res=200)
    par(mfrow=c(1,4), mai=c(0.05,0.05,0.05,0.05))
    
    Ulim = max(abs(range(pcs$U)))
    Zlim = max(abs(range(rbind(var_factors$Z, entr_factors$Z, entr2_factors$Z))))
    
    plot_pair(rotate_random(pcs$U), xlim=c(-Ulim,Ulim), ylim=c(-Ulim,Ulim))
    plot_pair(var_factors$Z, xlim=c(-Zlim,Zlim), ylim=c(-Zlim,Zlim))
    plot_pair(entr2_factors$Z, xlim=c(-Zlim,Zlim), ylim=c(-Zlim,Zlim))
    plot_pair(entr_factors$Z, xlim=c(-Zlim,Zlim), ylim=c(-Zlim,Zlim))
    dev.off()
  }
}