# This file implements the following orthogonal rotation
# methods:
#   1. Varimax
#   2. Entromin
#   3. Entromin2 (second-order Entromin approximation)

# Both pairwise optimization and gradient projection
# implementations are available

# Varimax and general rotation procedure code is based on
# code from https://github.com/RoheLab/vsp-paper
# ---------------------------------------------------------


library(rARPACK)
library(dplyr)

source("functions/rohe2020.R")
source("functions/util.R")



# Shared functions
# ---------------------------------------------------------
# These functions are intended to be private and not called outside of the file


# Computes the optimal rotation of a matrix given a rotation method
# output: rotation matrix
get_rotation = function(rotation_f, U, fast=T, ...) {
  # rotation_f: rotation method function
  # U: matrix to rotate
  # fast: T to compute rotation based on subset of U; F otherwise
  # ...: max_iters, verbose
  
  if (fast)
  {
    # Compute rotation based on top 10% of rows with greatest row norms
    lev = rowSums(U^2)
    R = rotation_f(U[lev>quantile(lev,0.9),], ...)
  } else {
    R = rotation_f(U, ...)
  }
  
  # Make sure columns of rotated matrix skew positive
  R = skew_positive(U,R)
  
  return(R)
}


# Wrapper function to rename vsp_recenter
# output: matrix of means for recentering
recenter_rotated = function(means, V, scree, Rz) {
  # means: means used in centering
  # V: right singular values from SVD
  # scree: singular values from SVD
  # Rz: rotation matrix
  
  return(vsp_recenter(means, V, scree, Rz))
}


# Computes the rotated principal components (PC)
# output:
#   Z: rotated matrix
#   Rz: rotation matrix
#   P: permutation matrix based on given sorting function
rotate = function(rotation_f, sort_f, pcs, recenter=F, ...) {
  # rotation_f: rotation method function
  # sort_f: function to sort columns of rotated PCs
  # pcs: PC object
  # recenter: T to recenter rotated matrix; F otherwise
  # ...: fast, max_iters, verbose
  
  U = pcs$U
  n = nrow(U)
  k = ncol(U)
  
  # Compute rotation matrix and factor estimates
  Rz = get_rotation(rotation_f, U, ...)
  Z = U %*% Rz * sqrt(n)
  
  # Create permutation matrix for sorting factor estimates
  P = diag(k)
  if (!is.null(sort_f))
  {
    sort_obj = apply(Z, 2, sort_f)
    P = P[sort(sort_obj,decreasing=T,index.return=T)$ix,]
  }
  
  # Recenter the factor estimates if requested
  if (is.null(recenter))
  {
    recenter = pcs$centered
  }
  if (recenter)
  {
    muZ = recenter_rotated(means=pcs$csA/n,pcs$V,pcs$scree,Rz) / sqrt(n)
    Z = Z + matrix(muZ,nrow=n,ncol=k,byrow=T)
  }
  
  rot_factors = list(Z=Z, Rz=Rz, P=P)
  return(rot_factors)
}


# Computes rotation based on gradient ascent pairwise optimization
# Note: only gradient convergence criterion is implemented
# output: rotation matrix
pairwise_gd = function(grad_f, U, eps=1e-5, max_iters=100, verbose=T, ...) {
  # grad_f: function for computing the gradient
  # U: matrix to rotate
  # eps: tolerance parameter
  # max_iters: maximum number of iterations
  # verbose: T to print intermediate progress; F otherwise
  # ...: drop_const
  
  n = nrow(U)
  k = ncol(U)
  R = diag(k)
  Z = U
  
  # Precompute gradients
  grads = matrix(0,k,k)
  for (j in 1:(k-1))
  {
    for (l in (j+1):k)
    {
      grads[j,l] = grad_f(Z[,c(j,l)])
    }
  }
  
  # Rotate pair of columns that have largest gradient
  for (i in 1:max_iters)
  {
    inds = which(abs(grads)==max(abs(grads)), arr.ind=T)
    j = inds[1]
    l = inds[2]
    
    # Perform gradient ascent
    Z2 = Z[,c(j,l)]
    iter = 1
    th = 0
    dt = grads[j,l]
    Rt = diag(2)
    while (abs(dt) > eps)
    {
      th = th + (1/sqrt(iter))*dt
      Rt = matrix(c(cos(th),sin(th),-sin(th),cos(th)), 2, 2)
      dt = grad_f(Z2, Rt)
      iter = iter + 1
    }
    grads[j,l] = dt
    
    # Update overall rotation matrix and rotated matrix
    R2 = diag(k)
    R2[c(j,l),c(j,l)] = Rt
    R = R %*% R2
    Z = Z %*% R2
    
    # Update pre-computed gradients
    for (q in 1:k)
    {
      if (q>j) {grads[j,q]=grad_f(Z[,c(j,q)])}
      if (q<l) {grads[q,l]=grad_f(Z[,c(q,l)])}
    }
    
    # Print progress
    if (verbose)
    {
      sprintf("Iter %i: cols=(%i,%i) ; V=%e ; H2=%.5f ; H=%.5f",
              i,j,l,
              varimax_objective(Z,...),
              round(entromin2_objective(Z,...),5),
              round(entromin_objective(Z),5)) %>% print
    }
    
    # Check convergence
    if (max(abs(grads)) < eps)
    {
      sprintf("Converged after %i iterations",i) %>% print
      break
    }
    if (i == max_iters)
    {
      sprintf("Failed to converge after %i iterations",i) %>% print
    }
  }
  return(R)
}



# Varimax-specific functions
# ---------------------------------------------------------

# Computes the value of the Varimax objective
# output: Varimax objective
varimax_objective = function(U, R=NULL, drop_const=F) {
  # U: input matrix
  # R: optional rotation matrix to apply to U
  # drop_const: T to drop constant terms from the objective; F otherwise
  
  n = nrow(U)
  k = ncol(U)
  if (is.null(R))
  {
    Z = U
  } else {
    Z = U %*% R
  }
  if (drop_const)
  {
    return(sum(Z^4))
  }
  Z2 = Z^2
  return(sum(Z2^2)/n - sum(colMeans(Z2)^2))
}


# Varimax rotation based on modified version of base R's varimax()
# output: rotation matrix
varimax = function(U, convergence="gradient", eps=1e-05,
                   max_iters=100, verbose=T, ...) {
  # U: matrix to rotate
  # convergence: {"objective","gradient"}
  # eps: tolerance parameter
  # max_iters: maximum number of iterations
  # verbose: T to print intermediate progress; F otherwise
  # ...: drop_const
  
  n = nrow(U)
  k = ncol(U)
  if (k<2) {return(U)}
  if (convergence!="gradient" & convergence!="objective")
  {
    stop("Convergence criterion unrecognized. Must be one of {\"objective\",\"gradient\"}")
  }
  
  # Helper function for computing objective at most once per iteration
  obj = NULL
  get_obj = function(Z, ...)
  {
    if (!is.null(obj)) {return(obj)}
    return(varimax_objective(Z, ...))
  }
  
  # Perform gradient projection
  R = diag(k)
  d = 0
  Z = U %*% R
  for (i in 1:max_iters)
  {
    B = t(U) %*% (Z^3 - Z%*%diag(drop(rep(1,n)%*%Z^2))/n)
    s = La.svd(B)
    R = s$u %*% s$vt
    Z = U %*% R
    obj = NULL
    
    # Print progress
    if (verbose)
    {
      obj = get_obj(Z, ...)
      sprintf("Iter %i: V=%e",i,obj) %>% print
    }
    
    # Check convergence
    converged = F
    dpast = d
    if (convergence == "objective")
    {
      d = get_obj(Z, ...)
      if (max(abs(d/dpast),abs(dpast/d))<1+eps) {converged=T}
    } else {
      # Gradient
      d = sum(s$d)
      if (d<dpast*(1+eps)) {converged=T}
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
  return(R)
}


# Computes the gradient of the Varimax objective for pairwise optimization
# output: gradient of Varimax objective
dvarimax = function(XY, R=NULL) {
  # XY: n x 2 matrix to compute gradient at
  # R: optional rotation matrix
  
  if (is.null(R))
  {
    XYR = XY
  } else {
    XYR = XY %*% R
  }
  return(4*sum(t(XYR[,1]*XYR[,2]) %*% XYR^2 %*% matrix(c(1,-1),2)))
}


# Wrapper function for finding Varimax rotation using pairwise optimization
# output: rotation matrix
varimax_pairwise = function(U, ...) {
  # U: matrix to rotate
  # ...: eps, max_iters, verbose
  
  return(pairwise_gd(dvarimax, U, ...))
}


# Wrapper function for computing the rotated PCs
# output: Z, Rz, P (see rotate())
varimax_rotate = function(pcs, pairwise=F, ...) {
  # pcs: PC object
  # pairwise: T to compute rotation using pairwise optimization; F otherwise
  # ...: recenter, fast, max_iters, verbose
  if (pairwise)
  {
    rotation_f = varimax_pairwise
  } else {
    rotation_f = varimax
  }
  var_factors = rotate(rotation_f, var, pcs, ...)
  return(var_factors)
}



# Entromin2-specific functions
# ---------------------------------------------------------

# Computes the value of the Entromin2 objective
# output: Entromin2 objective
entromin2_objective = function(U, R=NULL, drop_const=F) {
  # U: input matrix
  # R: optional rotation matrix to apply to U
  # drop_const: T to drop constant terms from the objective; otherwise F
  
  if (is.null(R))
  {
    Z = U
  } else {
    Z = U %*% R
  }
  const = 0
  if (!drop_const)
  {
    const = 3/2 * ncol(U)
  }
  return(sum(Z^6/2 - 2*Z^4) + const)
}


# Entromin2 rotation based on gradient projection algorithm
# output: rotation matrix
entromin2 = function(U, convergence="gradient", eps=1e-5,
                     max_iters=100, verbose=T, ...) {
  # U: matrix to rotate
  # convergence: {"objective","gradient"}
  # eps: tolerance parameter
  # max_iters: maximum number of iterations
  # verbose: T to print intermediate progress; otherwise F
  # ...: drop_const
  
  n = nrow(U)
  k = ncol(U)
  if (k<2) {return(U)}
  if (convergence!="gradient" & convergence!="objective")
  {
    stop("Convergence criterion unrecognized. Must be one of {\"objective\",\"gradient\"}")
  }
  
  # Helper function for computing objective at most once per iteration
  obj = NULL
  get_obj = function(Z, ...)
  {
    if (!is.null(obj)) {return(obj)}
    return(entromin2_objective(Z, ...))
  }
  
  # Perform gradient projection
  R = diag(k)
  Z = U %*% R
  d = 0
  for (i in 1:max_iters)
  {
    G = t(U) %*% (3*Z^5 - 8*Z^3)
    s = La.svd(G)
    R = s$u %*% s$vt
    Z = U %*% R
    obj = NULL
    
    # Print progress
    if (verbose)
    {
      obj = get_obj(Z, ...)
      sprintf("Iter %i: H2=%.5f",i,obj) %>% print
    }
    
    # Check convergence
    converged = F
    d_past = d
    if (convergence == "objective")
    {
      d = get_obj(Z, ...)
      if (max(abs(d/d_past),abs(d_past/d))<1+eps) {converged=T}
    } else {
      # Gradient
      d = sum(s$d)
      if (d<d_past*(1+eps)) {converged=T}
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
  return(R)
}


# Computes the gradient of the Entromin2 objective for pairwise optimization
# output: gradient of Entromin2 objective
dentromin2 = function(XY, R=NULL) {
  # XY: n x 2 matrix to compute gradient at
  # R: optional rotation matrix
  
  if (is.null(R))
  {
    XYR = XY
  } else {
    XYR = XY %*% R
  }
  tt = 3*XYR^5-8*XYR^3
  return(t(XYR[,2])%*%tt[,1] - t(XYR[,1]%*%tt[,2]))
}


# Wrapper function for finding Entromin2 rotation using pairwise optimization
# output: rotation matrix
entromin2_pairwise = function(U, ...) {
  # U: matrix to rotate
  # ...: eps, max_iters, verbose
  
  dentr2 = function(XY,R=NULL) {-dentromin2(XY,R)}
  return(pairwise_gd(dentr2, U, ...))
}


# Wrapper function for computing the rotated PCs
# output: Z, Rz, P (see rotate())
entromin2_rotate = function(pcs, pairwise=F, ...) {
  # pcs: PC object
  # pairwise: T to compute rotation using pairwise optimization; F otherwise
  # ...: recenter, fast, max_iters, verbose
  
  if (pairwise)
  {
    rotation_f = entromin2_pairwise
  } else {
    rotation_f = entromin2
  }
  entr2_factors = rotate(rotation_f, NULL, pcs, ...)
  return(entr2_factors)
}



# Entromin-specific functions
# ---------------------------------------------------------

# Computes the value of the Entromin objective
# output: Entromin objective
entromin_objective = function(U, R=NULL) {
  # U: input matrix
  # R: optional rotation matrix to apply to U
  
  if (is.null(R))
  {
    Z = U
  } else {
    Z = U %*% R
  }
  Z2 = Z^2
  Ze = Z2*log(Z2)
  return(-sum(replace(Ze,is.nan(Ze),0)))
}


# Entromin rotation based on gradient projection algorithm
# output: rotation matrix
entromin = function(U, convergence="gradient", alpha=1, eps=1e-5,
                    max_iters=100, verbose=F, ...) {
  # U: matrix to rotate
  # convergence: {"objective","gradient"}
  # alpha: step size parameter
  # eps: tolerance parameter
  # max_iters: maximum number of iterations
  # verbose: T to print intermediate progress; otherwise F
  # ...: drop_const
  
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
  for (i in 1:max_iters)
  {
    GH = Z*log(Z^2)
    G = t(U) %*% (-replace(GH,is.nan(GH),0)-Z)
    
    # Find appropriate step size
    a = a*2
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
  return(R)
}


# Computes the gradient of the Entromin objective for pairwise optimization
# output: gradient of Entromin objective
dentromin = function(XY, R=NULL) {
  # XY: n x 2 matrix to compute gradient at
  # R: optional rotation matrix
  
  if (is.null(R))
  {
    XYR = XY
  } else {
    XYR = XY %*% R
  }
  logXYR2 = log(XYR^2)
  return(t(XYR[,1]*XYR[,2]) %*% replace(logXYR2,-Inf,0) %*% matrix(c(-2,2),2))
}


# Wrapper function for finding Entromin rotation using pairwise optimization
# output: rotation matrix
entromin_pairwise = function(U, ...) {
  # U: matrix to rotate
  # ...: eps, max_iters, verbose
  
  dentr = function(XY,R=NULL) {-dentromin(XY,R)}
  return(pairwise_gd(dentr, U, ...))
}


# Wrapper function for computing the rotated PCs
# output: Z, Rz, P (see rotate())
entromin_rotate = function(pcs, pairwise=F, ...) {
  # pcs: PC objective
  # pairwise: T to compute rotation using pairwise optimization; F otherwise
  # ...: recenter, fast, max_iters, verbose
  
  if (pairwise)
  {
    rotation_f = entromin_pairwise
  } else {
    rotation_f = entromin
  }
  entr_factors = rotate(rotation_f, NULL, pcs, ...)
  return(entr_factors)
}