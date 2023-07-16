
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



rowrep <- function(X, ntimes){

  X[rep(seq_along(ntimes), ntimes), ]

}

