# Code for soft zero sparsity analysis on the simulated
# sparse Gaussian dataset

# Note: may want to suppress rotation convergence
#       messages before running to avoid spam
# ---------------------------------------------------------


set.seed(1)

rm(list=ls())

source("functions/dataset.R")
source("functions/orth_rotation.R")
source("functions/util.R")

# Control parameter
save_data = F # Save data to file

# Dataset parameters
N = 100  # Number of datasets
n = 1000 # Number of samples
d = 3    # Number of columns

# Rotation parameters
k = 2 # Final number of PCs

tau           = -1 # Regularization parameter for bigPCA
normalization = F  # Normalize data before PCA
centering     = F  # Center data before PCA
recenter      = F  # Recenter rotated factors
fast          = F  # Compute rotation based on 10% of rows

max_iters = 1000 # Max number of rotation iterations
verbose   = F    # Print intermediate progress



# Simulate datasets
# ---------------------------------------------------------

# Probabilities of zero
p0 = 0:9/10

raw_soft0 = c()
raw_sd = c()
var_soft0 = c()
var_sd = c()
entr2_soft0 = c()
entr2_sd = c()
entr_soft0 = c()
entr_sd = c()

for (p in p0)
{
  sprintf("Probability: %.1f",p) %>% print
  
  raw_p = c()
  var_p = c()
  entr2_p = c()
  entr_p = c()
  
  for (i in 1:N)
  {
    # Simulate new dataset
    A = get_sparse_normal(p, n=n, k=d)
    pcs = bigPCA(A, k, tau, normalization, centering)
    
    raw_p = c(raw_p, count_soft_zeros(pcs$U))
    
    # Rotate and count soft zeros
    R = varimax_rotate(pcs, recenter=recenter, fast=fast,
                       max_iters=max_iters, verbose=verbose)$Rz
    var_p = c(var_p, count_soft_zeros(pcs$U%*%R))
    
    R = entromin2_rotate(pcs, recenter=recenter, fast=fast,
                         max_iters=max_iters, verbose=verbose)$Rz
    entr2_p = c(entr2_p, count_soft_zeros(pcs$U%*%R))
    
    R = entromin_rotate(pcs, recenter=recenter, fast=fast,
                        max_iters=max_iters, verbose=verbose)$Rz
    entr_p = c(entr_p, count_soft_zeros(pcs$U%*%R))
  }
  
  raw_soft0 = c(raw_soft0, mean(raw_p))
  raw_sd = c(raw_sd, sd(raw_p))
  var_soft0 = c(var_soft0, mean(var_p))
  var_sd = c(var_sd, sd(var_p))
  entr2_soft0 = c(entr2_soft0, mean(entr2_p))
  entr2_sd = c(entr2_sd, sd(entr2_p))
  entr_soft0 = c(entr_soft0, mean(entr_p))
  entr_sd = c(entr_sd, sd(entr_p))
}

# Clean up and save data
df = data.frame(p0 = p0,
                # Mean soft zeros
                raw = raw_soft0,
                var = var_soft0,
                entr2 = entr2_soft0,
                entr = entr_soft0,
                # Standard deviation
                raw_sd = raw_sd,
                var_sd = var_sd,
                entr2_sd = entr2_sd,
                entr_sd = entr_sd)
if (save_data) {write.table(df,"sparsity.txt",quote=F,row.names=F)}