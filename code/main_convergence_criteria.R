# Code for convergence criterion analysis on the simulated
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
save_data = T # Save data to file

# Dataset parameters
p0 = 0.5   # Probability of zero
N  = 100   # Number of datasets
n  = 50000 # Number of samples
d  = 20    # Number of columns

# Rotation parameters
k          = 20 # Final number of PCs
drop_const = T  # Drop constants from objectives

tau           = -1 # Regularization parameter for bigPCA
normalization = F  # Normalize data before PCA
centering     = F  # Center data before PCA
recenter      = F  # Recenter rotated factors
fast          = F  # Compute rotation based on 10% of rows

max_iters = 1000 # Max number of rotation iterations
verbose   = F    # Print intermediate progress



# Convergence criteria analysis
# ---------------------------------------------------------

# C1
var_df = NULL
entr2_df = NULL
entr_df = NULL

# C2
varo_df = NULL
entr2o_df = NULL
entro_df = NULL

for (i in 1:N)
{
  sprintf("Starting dataset %i",i) %>% print
  
  # Simulate new dataset
  A = get_sparse_normal(p0, n=n, k=d, all_sparse=T)
  pcs = bigPCA(A, k, tau, normalization, centering)
  
  # Varimax
  var_factors = varimax_rotate(pcs, recenter=recenter, drop_const=drop_const,
                               fast=fast, max_iters=max_iters, verbose=verbose)
  var_stats = save_stats(pcs$U%*%var_factors$Rz, drop_const=drop_const)
  var_df = rbind(var_df, unlist(var_stats))
  
  var_factors = varimax_rotate(pcs, recenter=recenter,
                               drop_const=drop_const, convergence="objective",
                               fast=fast, max_iters=max_iters, verbose=verbose)
  var_stats = save_stats(pcs$U%*%var_factors$Rz, drop_const=drop_const)
  varo_df = rbind(varo_df, unlist(var_stats))
  
  # Entromin2
  entr2_factors = entromin2_rotate(pcs, recenter=recenter, drop_const=drop_const,
                                   fast=fast, max_iters=max_iters, verbose=verbose)
  entr2_stats = save_stats(pcs$U%*%entr2_factors$Rz, drop_const=drop_const)
  entr2_df = rbind(entr2_df, unlist(entr2_stats))
  
  entr2_factors = entromin2_rotate(pcs, recenter=recenter,
                                   drop_const=drop_const, convergence="objective",
                                   fast=fast, max_iters=max_iters, verbose=verbose)
  entr2_stats = save_stats(pcs$U%*%entr2_factors$Rz, drop_const=drop_const)
  entr2o_df = rbind(entr2o_df, unlist(entr2_stats))
  
  # Entromin
  entr_factors = entromin_rotate(pcs, recenter=recenter, drop_const=drop_const,
                                 fast=fast, max_iters=max_iters, verbose=verbose)
  entr_stats = save_stats(pcs$U%*%entr_factors$Rz, drop_const=drop_const)
  entr_df = rbind(entr_df, unlist(entr_stats))
  
  entr_factors = entromin_rotate(pcs, recenter=recenter,
                                 drop_const=drop_const, convergence="objective",
                                 fast=fast, max_iters=max_iters, verbose=verbose)
  entr_stats = save_stats(pcs$U%*%entr_factors$Rz, drop_const=drop_const)
  entro_df = rbind(entro_df, unlist(entr_stats))
}

# Save data
if (save_data)
{
  write.table(var_df,"var_c1.txt",row.names=F,quote=F)
  write.table(entr2_df,"entr2_c1.txt",row.names=F,quote=F)
  write.table(entr_df,"entr_c1.txt",row.names=F,quote=F)
  
  write.table(varo_df,"var_c2.txt",row.names=F,quote=F)
  write.table(entr2o_df,"entr2_c2.txt",row.names=F,quote=F)
  write.table(entro_df,"entr_c2.txt",row.names=F,quote=F)
}