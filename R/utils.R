
# make a P-spline basis with K basis functions

psBasis <- function(x, K = length(x), spline.degree = 3, diff.ord = 2,
                    knots = NULL) {
  if (is.null(knots)) {
    knots.no <- K - spline.degree + 1
    xl <- min(x)
    xr <- max(x)
    xmin <- xl - (xr - xl) / 100
    xmax <- xr + (xr - xl) / 100
    dx <- (xmax - xmin) / (knots.no - 1)
    knots <- seq(xmin - spline.degree * dx, xmax + spline.degree * dx, by = dx)
  }
  X <- splines::spline.des(knots, x, spline.degree + 1, outer.ok = TRUE)$design
  P <- diag(K) # precision
  if (diff.ord > 0) {
    for (d in 1:diff.ord) P <- diff(P)
    P <- crossprod(P)
  }
  return(list(
    X = X, P = P, knots = knots, K = K, spline.degree = spline.degree,
    diff.ord = diff.ord
  ))
}


# repeat each row i of matrix X by ntimes[i] times

rowrep <- function(X, ntimes){

  X[rep(seq_along(ntimes), ntimes), ]

}

# compute posterior mean and (1-alpha)\% pointwise credible intervals from mcmc draws of function

postf <- function(f_post, alpha = .05){

  f_ci_l <- apply(f_post[,], 1, function(x) stats::quantile(x, c(alpha/2)))
  f_ci_u <- apply(f_post[,], 1, function(x) stats::quantile(x, c(1 - alpha/2)))
  f_mean <- rowMeans(f_post[,])

  output <- list(f_mean = f_mean,
                 f_ci_l = f_ci_l,
                 f_ci_u = f_ci_u)
  output
}


#' Compute Simultaneous Credible Bands
#'
#' Compute (1-alpha)\% credible BANDS for a function based on MCMC samples using Crainiceanu et al. (2007)
#'
#' @param sampFuns \code{Nsims x m} matrix of \code{Nsims} MCMC samples and \code{m} points along the curve
#' @param alpha confidence level
#'
#' @return \code{m x 2} matrix of credible bands; the first column is the lower band, the second is the upper band
#'
#' @note The input needs not be curves: the simultaneous credible "bands" may be computed
#' for vectors. The resulting credible intervals will provide joint coverage at the (1-alpha)%
#' level across all components of the vector.
#'
#' @export

credBands = function(sampFuns, alpha = .05){

  N = nrow(sampFuns); m = ncol(sampFuns)

  # Compute pointwise mean and SD of f(x):
  Efx = colMeans(sampFuns); SDfx = apply(sampFuns, 2, stats::sd)

  # Compute standardized absolute deviation:
  Standfx = abs(sampFuns - tcrossprod(rep(1, N), Efx))/tcrossprod(rep(1, N), SDfx)

  # And the maximum:
  Maxfx = apply(Standfx, 1, max)

  # Compute the (1-alpha) sample quantile:
  Malpha = stats::quantile(Maxfx, 1-alpha)

  # Finally, store the bands in a (m x 2) matrix of (lower, upper)
  t(cbind(Efx - Malpha*SDfx, Efx + Malpha*SDfx))
}
