#----------------------------------------------------------------------------
# Generate simulated longitudinal functional data
#----------------------------------------------------------------------------

#' Simulated Longitudinal Functional Data
#'
#' Generate simulated longitudinal functional data and scalar covariates.
#'
#' @param N number of subjects.
#' @param Mi integer or vector of length N specifying the number of replicates of each subject.
#' @param L number of non-zero scalar covariates.
#' @param Tn total number of time points.
#' @param K number of basis functions in final resulting basis (with K-2 reparameterized P-spline functions).
#' @param sig_noise observation level variance.
#' @param sig_with within-curve smooth variance.
#' @param sig_bet between-subjects smooth variance.
#' @param sig_alpha fixed effects variance.
#' @param seed (optional) sets the seed.
#'
#' @return list with the following elements:
#' * Y: simulated longitudinal functional data
#' * X: simulated design matrix
#' * z: vector specifying group memberships of each curve
#' * alphaf_true: true functional regression coefficients from which the data was generated from
#' @export
#'
#' @examples
#' simdata <- sim_lfosr_data(N = 20, Mi = 5, L = 3, Tn = 100)
#'
#' @importFrom stats poly rnorm
sim_lfosr_data <- function(N, Mi, L, Tn, K = 5, sig_noise = .1, sig_with = 1, sig_bet = 1, sig_alpha = 1, seed = NULL){

  if(length(Mi) == 1){
    Mi <- rep(Mi, N)
    Mi <<- Mi
  }

  tau <- seq(0, 1, by = 1/(Tn-1))

  b = fda::create.bspline.basis(rangeval = c(1, Tn), breaks = seq(1, Tn, 10), norder = 4)
  Bmat = fda::eval.basis(1:Tn, b)
  Kbmat <- ncol(Bmat)
  P <- psBasis(1:Kbmat)$P
  B <- cbind(1/sqrt(Tn), poly(tau, 1), eigen(Bmat %*% ((MASS::ginv(P)) %*% t(Bmat)), symmetric = T)$vectors[,1:(K-2)])

  sig_w <- sig_with

  sig_ga <- sig_bet

  sig_alpha <- sig_alpha

  set.seed(seed)
  alpha_true <- matrix(NA, nrow = K, ncol = L+1)


  alpha_true[,1] <- 1/K

  for(i in 2:(L+1)){

    a <- mgcv::rmvn(1, c(rep(0, K)), (sig_alpha)*diag(K))
    alpha_true[,i] <- a

  }
  alphaf_true <- B%*%alpha_true

  w_true <- matrix(rnorm(sum(Mi)*K, mean = 0, sqrt(sig_w)), nrow = sum(Mi), ncol = K)
  ga_true <- matrix(rnorm(N*K, mean = 0, sqrt(sig_ga)), nrow = N, ncol = K)
  gaf_true <-  B%*%t(ga_true)
  noise <- matrix(rnorm(sum(Mi)*Tn, mean = 0, sd = sqrt(sig_noise)), nrow = Tn, ncol = sum(Mi))

  X <- cbind(rep(1, N), scale(matrix(rnorm(N*L), nrow = N, ncol = L)))

  z <- rep(1:N, Mi)

  Y <- B%*%t(w_true + rowrep(X%*%t(alpha_true), Mi)) + t(rowrep(t(gaf_true), Mi)) + noise

  simdata <- list(Y = Y,
                  X = X,
                  z = z,
                  alphaf = alphaf_true)

  return(simdata)

}
