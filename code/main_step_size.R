# Code for step size analysis on the simulated sparse
# Gaussian dataset
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
verbose   = T    # Print intermediate progress



# Redefine entromin() for convenience
# ---------------------------------------------------------

# Output intermediate objectives, step sizes and iterations
entromin_ss = function(U, convergence="gradient", alpha=1, eps=1e-5,
                       max_iters=100, verbose=F, ...)
{
  n = nrow(U)
  k = ncol(U)
  if (k<2) {return(U)}
  if (convergence != "gradient" & convergence != "objective")
  {
    stop("Convergence criterion unrecognized. Must be one of {\"objective\",\"gradient\"}")
  }
  
  # Perform gradient projection
  R = diag(k)
  Z = U %*% R
  if (convergence=="gradient") {d=0}
  H = entromin_objective(Z)
  a = alpha
  
  objs = c(H)
  iters = c(0)
  alphas = c(a)
  for (i in 1:max_iters)
  {
    GH = Z*log(Z^2)
    G = t(U) %*% (-replace(GH,is.nan(GH),0)-Z)
    
    # Find appropriate step size
    # a = a*2 # For this analysis, do not try larger step sizes
    for (j in 1:10)
    {
      s = La.svd(R - a*G)
      R_new = s$u %*% s$vt
      Z_new = U %*% R_new
      H_new = entromin_objective(Z_new)
      if (H_new<H) {break}
      # Step size too large
      a = a/2
    }
    
    # Check convergence
    converged = F
    if (convergence == "objective")
    {
      if (max(H/H_new,H_new/H)<1+eps) {converged=T}
    } else {
      # Gradient
      d_past = d
      d = sum(s$d)
      if (d<d_past*(1+eps)) {converged=T}
    }
    
    # Update progress
    R = R_new
    Z = Z_new
    H = H_new
    
    objs = c(objs, H)
    iters = c(iters, i)
    alphas = c(alphas, a)
    
    # Print progress
    if (verbose)
    {
      sprintf("Iter %i: alpha=%.5f ; H=%.5f",i,a,round(H,5)) %>% print
    }
    
    if (converged)
    {
      sprintf("Converged after %i iterations",i) %>% print
      break
    }
    if (i == max_iters)
    {
      sprintf("Failed to converge after %i iterations",i) %>% print
    }
  }
  return(list(R=R, objs=objs, iters=iters, alphas=alphas))
}



# Data preprocessing
# ---------------------------------------------------------

A = get_sparse_normal(p0, n=n, k=d, all_sparse=T)
pcs = bigPCA(A, k, tau, normalization, centering)




# Step size analysis
# ---------------------------------------------------------

ents1 = entromin_ss(pcs$U, alpha=1, max_iters=max_iters, verbose=verbose)
df1 = data.frame(iter=ents1$iters, obj=ents1$objs, alpha=ents1$alphas)

ents2 = entromin_ss(pcs$U, alpha=0.1, max_iters=max_iters, verbose=verbose)
df2 = data.frame(iter=ents2$iters, obj=ents2$objs, alpha=ents2$alphas)

if (save_data)
{
  write.table(df1, "a1.txt", quote=F, row.names=F)
  write.table(df2, "a01.txt", quote=F, row.names=F)
}