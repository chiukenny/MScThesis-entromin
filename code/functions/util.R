# This file implements all other functions (helper,
# plotting, etc.)

# Some plotting code is based on code from
# https://github.com/RoheLab/vsp-paper
# ---------------------------------------------------------



# Utility
# ---------------------------------------------------------

# Rotates a matrix randomly
# output: rotated matrix
rotate_random = function(X) {
  # X: matrix to rotate
  
  n = dim(X)[1]
  k = dim(X)[2]
  for (i in 1:(k-1))
  {
    for (j in (i+1):k)
    {
      # Randomly rotate every pair of columns
      t = runif(1) * 2*pi
      R = matrix(c(cos(t), sin(t), -sin(t), cos(t)), 2, 2)
      X[,c(i,j)] = X[,c(i,j)] %*% R
    }
  }
  return(X)
}


# Sets sign of each column to have positive third sample moment
# output: one of
#   X: sign-flipped matrix
#   R: sign-flipped rotation matrix (if initial rotation matrix given)
skew_positive = function(X, R=NULL) {
  # X: matrix to skew positive
  # R: optional rotation matrix to skew positive
  
  k = dim(X)[2]
  if (is.null(R))
  {
    X = X %*% diag(sign(colSums(X^3)))
    return(X)
  } else {
    R = R %*% diag(sign(colSums((X %*% R)^3)))
    return(R)
  }
}


# Discards unneeded PCs
# output: updated PC object
keep_factors = function(pcs, keep) {
  # pcs: PC object
  # keep: indices of PCs OR number of leading PCs to keep
  
  if (length(keep) == 1)
  {
    keep = 1:keep
  }
  k = length(keep)
  
  keep_pcs = pcs
  keep_pcs$U = pcs$U[,keep]
  keep_pcs$V = pcs$V[,keep]
  keep_pcs$scree = pcs$scree[keep]
  return(keep_pcs)
}


# Counts soft zeros in given matrix
# output: number of soft zeros
count_soft_zeros = function(X, eps=1e-5) {
  # X: matrix to count
  # eps: soft zero threshold
  
  return(sum(abs(X) < eps))
}


# Print objective values for given matrix
get_stats = function(X, eps=1e-5, drop_const=F) {
  # X: matrix to compute values for
  # eps: threshold for soft zeros
  # drop_const: T to drop constants from objectives; F otherwise
  
  print(paste("Varimax:", varimax_objective(X,drop_const=drop_const)))
  print(paste("Entromin2:", entromin2_objective(X,drop_const=drop_const)))
  print(paste("Entromin:", entromin_objective(X)))
  print(paste("Soft zeros:", count_soft_zeros(X,eps=eps)))
}


# Save objective values for given matrix
# output:
#   var: Varimax
#   entr2: Entromin2
#   entr: Entromin
#   zeros: soft zeros
save_stats = function(X, eps=1e-5, drop_const=F) {
  # X: matrix to compute values for
  # eps: threshold for soft zeros
  # drop_const: T to drop constants from objectives; F otherwise
  
  return(list(var=varimax_objective(X,drop_const=drop_const),
              entr2=entromin2_objective(X,drop_const=drop_const),
              entr=entromin_objective(X),
              zeros=count_soft_zeros(X,eps=eps)))
}



# Plotting
# ---------------------------------------------------------

# Plot red axes and subset of points
plot_pair = function(X, samples=NULL, ...)
{
  if (is.null(samples))
  {
    samples = 1:nrow(X)
  }
  x = X[samples,1]
  y = X[samples,2]
  
  plot(NULL, xlab="", ylab="", xaxt="n", yaxt="n", ...)
  abline(h=0, col="red", ...)
  abline(v=0, col="red", ...)
  points(x, y, pch=".", cex=1.5, ...)
}


# Helper function to plot grey axes and points
plot_scatter = function(x, y, ...)
{
  abline(h=0, col="grey", ...)
  abline(v=0, col="grey", ...)
  points(x, y, pch=".", ...)
}


# Plot different pairwise data for upper and lower with diagonal
plot_diff_pairs = function(X, Y, k, samples)
{
  nX = nrow(X[samples,])
  nY = nrow(Y[samples,])
  XY = rbind(X[samples,1:k], Y[samples,1:k])
  
  pairs_highlight(XY, labels=NULL,
        lower.panel=function(x,y,...) {
          Xx = x[seq_len(nX)]
          Xy = y[seq_len(nX)]
          plot_scatter(Xx, Xy, cex=1.5, ...)
        },
        upper.panel=function(x,y,...) {
          Yx = x[(nX+1):(nX+nY)]
          Yy = y[(nX+1):(nX+nY)]
          plot_scatter(Yx, Yy, cex=1.5, ...)
        },
        diag.panel=function(x,...) {
          Xx = x[seq_len(nX)]
          Yx = x[(nX+1):(nX+nY)]
          abline(a=0, b=1, col="grey")
          plot_scatter(Xx, Yx, cex=1.5, ...)
        },
        xaxt="n", yaxt="n", oma=c(0,0,0,0))
}


# Plot different pairwise data for upper and lower without diagonal
plot_rot_unrot = function(X, rotX, k, samples)
{
  nX = nrow(X[samples,])
  nrX = nrow(rotX[samples,])
  XrX = rbind(X[samples,1:k], rotX[samples,1:k])
  
  pairs_no_diag(XrX, labels=NULL, label.pos=NULL,
                lower.panel=function(x,y,...) {
                  Xx = x[seq_len(nX)]
                  Xy = y[seq_len(nX)]
                  usr = par("usr")
                  on.exit(par(usr))
                  i = par("mfg")[1]
                  j = par("mfg")[2]
                  par(usr=c(range(X[,j]),range(X[,2:i])))
                  plot_scatter(Xx, Xy, cex=1.5, ...)
                },
                upper.panel=function(x,y,...) {
                  rXx = x[(nX+1):(nX+nrX)]
                  rXy = y[(nX+1):(nX+nrX)]
                  usr = par("usr")
                  on.exit(par(usr))
                  i = par("mfg")[1]
                  j = par("mfg")[2]
                  par(usr=c(range(rotX[,j]),range(rotX[,(i+1):k])))
                  plot_scatter(rXx, rXy, cex=1.5, ...)
                },
                xaxt="n", yaxt="n", oma=c(0,0,0,0))
}


# Custom definition of pairs() to remove the boxes on
# the diagonal + coloured borders
pairs_no_diag <- 
  function (x, labels, panel = points, ..., lower.panel = panel, 
            upper.panel = panel, diag.panel = NULL, text.panel = textPanel, 
            label.pos = 0.5 + has.diag/3, cex.labels = NULL, font.labels = 1, 
            row1attop = TRUE, gap = 1) 
  {
    textPanel <- function(x = 0.5, y = 0.5, txt, cex, font) text(x, 
                                                                 y, txt, cex = cex, font = font)
    localAxis <- function(side, x, y, xpd, bg, col = NULL, main, 
                          oma, ...) {
      if (side%%2 == 1) 
        Axis(x, side = side, xpd = NA, ...)
      else Axis(y, side = side, xpd = NA, ...)
    }
    localPlot <- function(..., main, oma, font.main, cex.main) plot(...)
    localLowerPanel <- function(..., main, oma, font.main, cex.main) lower.panel(...)
    localUpperPanel <- function(..., main, oma, font.main, cex.main) upper.panel(...)
    localDiagPanel <- function(..., main, oma, font.main, cex.main) diag.panel(...)
    dots <- list(...)
    nmdots <- names(dots)
    if (!is.matrix(x)) {
      x <- as.data.frame(x)
      for (i in seq_along(names(x))) {
        if (is.factor(x[[i]]) || is.logical(x[[i]])) 
          x[[i]] <- as.numeric(x[[i]])
        if (!is.numeric(unclass(x[[i]]))) 
          stop("non-numeric argument to 'pairs'")
      }
    }
    else if (!is.numeric(x)) 
      stop("non-numeric argument to 'pairs'")
    panel <- match.fun(panel)
    if ((has.lower <- !is.null(lower.panel)) && !missing(lower.panel)) 
      lower.panel <- match.fun(lower.panel)
    if ((has.upper <- !is.null(upper.panel)) && !missing(upper.panel)) 
      upper.panel <- match.fun(upper.panel)
    if ((has.diag <- !is.null(diag.panel)) && !missing(diag.panel)) 
      diag.panel <- match.fun(diag.panel)
    if (row1attop) {
      tmp <- lower.panel
      lower.panel <- upper.panel
      upper.panel <- tmp
      tmp <- has.lower
      has.lower <- has.upper
      has.upper <- tmp
    }
    nc <- ncol(x)
    if (nc < 2) 
      stop("only one column in the argument to 'pairs'")
    has.labs <- TRUE
    if (missing(labels)) {
      labels <- colnames(x)
      if (is.null(labels)) 
        labels <- paste("var", 1L:nc)
    }
    else if (is.null(labels)) 
      has.labs <- FALSE
    oma <- if ("oma" %in% nmdots) 
      dots$oma
    else NULL
    main <- if ("main" %in% nmdots) 
      dots$main
    else NULL
    if (is.null(oma)) {
      oma <- c(4, 4, 4, 4)
      if (!is.null(main)) 
        oma[3L] <- 6
    }
    opar <- par(mfrow = c(nc, nc), mar = rep.int(gap/2, 4), oma = oma)
    on.exit(par(opar))
    dev.hold()
    on.exit(dev.flush(), add = TRUE)
    for (i in if (row1attop) 
      1L:nc
      else nc:1L) for (j in 1L:nc) {
        localPlot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, 
                  type = "n", ...)
        if (i == j || (i < j && has.lower) || (i > j && has.upper)) {
          if (i == 1 && (!(j%%2) || !has.upper || !has.lower))
          localAxis(1 + 2 * row1attop, x[, j], x[, i], ...)
          if (j == nc && (i%%2 || !has.upper || !has.lower))
            localAxis(4, x[, j], x[, i], ...)
          mfg <- par("mfg")
          if (i == j) {
            if (has.diag) 
              localDiagPanel(as.vector(x[, i]), ...)
            if (has.labs) {
              par(usr = c(0, 1, 0, 1))
              if (is.null(cex.labels)) {
                l.wid <- strwidth(labels, "user")
                cex.labels <- max(0.8, min(2, 0.9/max(l.wid)))
              }
              text.panel(0.5, label.pos, labels[i], cex = cex.labels, 
                         font = font.labels)
            }
          }
          else if (i < j) {
            box(col="orange", lwd=2)
            localLowerPanel(as.vector(x[, j]), as.vector(x[, 
                                                           i]), ...)
          }
          else {
            box(col="blue", lwd=2)
            localUpperPanel(as.vector(x[, j]), as.vector(x[, 
                                                              i]), ...)
          }
          if (any(par("mfg") != mfg)) 
            stop("the 'panel' function made a new plot")
        }
        else par(new = FALSE)
      }
    if (!is.null(main)) {
      font.main <- if ("font.main" %in% nmdots) 
        dots$font.main
      else par("font.main")
      cex.main <- if ("cex.main" %in% nmdots) 
        dots$cex.main
      else par("cex.main")
      mtext(main, 3, 3, TRUE, 0.5, cex = cex.main, font = font.main)
    }
    invisible(NULL)
  }


# Custom definition of pairs() to colour upper and lower borders
pairs_highlight <- 
  function (x, labels, panel = points, ..., lower.panel = panel, 
            upper.panel = panel, diag.panel = NULL, text.panel = textPanel, 
            label.pos = 0.5 + has.diag/3, cex.labels = NULL, font.labels = 1, 
            row1attop = TRUE, gap = 1) 
  {
    textPanel <- function(x = 0.5, y = 0.5, txt, cex, font) text(x, 
                                                                 y, txt, cex = cex, font = font)
    localAxis <- function(side, x, y, xpd, bg, col = NULL, main, 
                          oma, ...) {
      if (side%%2 == 1) 
        Axis(x, side = side, xpd = NA, ...)
      else Axis(y, side = side, xpd = NA, ...)
    }
    localPlot <- function(..., main, oma, font.main, cex.main) plot(...)
    localLowerPanel <- function(..., main, oma, font.main, cex.main) lower.panel(...)
    localUpperPanel <- function(..., main, oma, font.main, cex.main) upper.panel(...)
    localDiagPanel <- function(..., main, oma, font.main, cex.main) diag.panel(...)
    dots <- list(...)
    nmdots <- names(dots)
    if (!is.matrix(x)) {
      x <- as.data.frame(x)
      for (i in seq_along(names(x))) {
        if (is.factor(x[[i]]) || is.logical(x[[i]])) 
          x[[i]] <- as.numeric(x[[i]])
        if (!is.numeric(unclass(x[[i]]))) 
          stop("non-numeric argument to 'pairs'")
      }
    }
    else if (!is.numeric(x)) 
      stop("non-numeric argument to 'pairs'")
    panel <- match.fun(panel)
    if ((has.lower <- !is.null(lower.panel)) && !missing(lower.panel)) 
      lower.panel <- match.fun(lower.panel)
    if ((has.upper <- !is.null(upper.panel)) && !missing(upper.panel)) 
      upper.panel <- match.fun(upper.panel)
    if ((has.diag <- !is.null(diag.panel)) && !missing(diag.panel)) 
      diag.panel <- match.fun(diag.panel)
    if (row1attop) {
      tmp <- lower.panel
      lower.panel <- upper.panel
      upper.panel <- tmp
      tmp <- has.lower
      has.lower <- has.upper
      has.upper <- tmp
    }
    nc <- ncol(x)
    if (nc < 2) 
      stop("only one column in the argument to 'pairs'")
    has.labs <- TRUE
    if (missing(labels)) {
      labels <- colnames(x)
      if (is.null(labels)) 
        labels <- paste("var", 1L:nc)
    }
    else if (is.null(labels)) 
      has.labs <- FALSE
    oma <- if ("oma" %in% nmdots) 
      dots$oma
    else NULL
    main <- if ("main" %in% nmdots) 
      dots$main
    else NULL
    if (is.null(oma)) {
      oma <- c(4, 4, 4, 4)
      if (!is.null(main)) 
        oma[3L] <- 6
    }
    opar <- par(mfrow = c(nc, nc), mar = rep.int(gap/2, 4), oma = oma)
    on.exit(par(opar))
    dev.hold()
    on.exit(dev.flush(), add = TRUE)
    for (i in if (row1attop) 
      1L:nc
      else nc:1L) for (j in 1L:nc) {
        localPlot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, 
                  type = "n", ...)
        if (i == j || (i < j && has.lower) || (i > j && has.upper)) {
          if (i == 1 && (!(j%%2) || !has.upper || !has.lower))
            localAxis(1 + 2 * row1attop, x[, j], x[, i], ...)
          if (j == nc && (i%%2 || !has.upper || !has.lower))
            localAxis(4, x[, j], x[, i], ...)
          mfg <- par("mfg")
          if (i == j) {
            box(col="darkgreen", lwd=2)
            if (has.diag) 
              localDiagPanel(as.vector(x[, i]), ...)
            if (has.labs) {
              par(usr = c(0, 1, 0, 1))
              if (is.null(cex.labels)) {
                l.wid <- strwidth(labels, "user")
                cex.labels <- max(0.8, min(2, 0.9/max(l.wid)))
              }
              text.panel(0.5, label.pos, labels[i], cex = cex.labels, 
                         font = font.labels)
            }
          }
          else if (i < j) {
            box(col="orange", lwd=2)
            localLowerPanel(as.vector(x[, j]), as.vector(x[, 
                                                           i]), ...)
          }
          else {
            box(col="blue", lwd=2)
            localUpperPanel(as.vector(x[, j]), as.vector(x[, 
                                                           i]), ...)
          }
          if (any(par("mfg") != mfg)) 
            stop("the 'panel' function made a new plot")
        }
        else par(new = FALSE)
      }
    if (!is.null(main)) {
      font.main <- if ("font.main" %in% nmdots) 
        dots$font.main
      else par("font.main")
      cex.main <- if ("cex.main" %in% nmdots) 
        dots$cex.main
      else par("cex.main")
      mtext(main, 3, 3, TRUE, 0.5, cex = cex.main, font = font.main)
    }
    invisible(NULL)
  }