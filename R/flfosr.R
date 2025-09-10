#----------------------------------------------------------------------------
# Fast longitudinal function-on-scalar regression (FLFOSR)
#----------------------------------------------------------------------------

#' Fast Longitudinal Function-On-Scalar Regression
#'
#' Performs MCMC estimation for Bayesian longitudinal function-on-scalar regression.
#'
#' @param Y matrix of functional responses. Each column corresponds to one curve from a subject and each row is a measurement at a single time point on a common grid.
#' @param X design matrix of fixed effects using `model.matrix`. It should include column of 1s for the intercept.
#' @param z vector specifying group memberships of each curve, e.g. if there are 20 subjects with 5 repeated measurements each, `z=rep(1:20, each = 5)`.
#' @param k number of basis functions
#' @param S number of total MCMC iterations
#' @param S_burn number of initial MCMC iterations to discard for warm-up
#' @param varhyp a list containing hyperparameter values for the gamma variance priors. It must contain the following named elements:
#' * a_a, b_a: shape and rate hyperparameters for fixed effects variance.
#' * a_g, b_g: shape and rate hyperparameters for between-groups variance.
#' * a_w, b_w: shape and rate hyperparameters for within-group variance.
#'
#' @return A list with the following elements:
#' * X: the inputed design matrix
#' * B: orthogonalized basis matrix
#' * w_post: posterior draws of curve specific random effects coefficients
#' * ga_post: posterior draws of the subject specific random effects coefficients
#' * alpha_post: posterior draws of the fixed effect coefficients
#' * alphaf_post: posterior draws of the fixed effect functions
#' @export
#'
#' @examples
#' simdata <- sim_lfosr_data(N = 20, Mi = 5, L = 3, Tn = 100)
#'
#' m1 <- flfosr(simdata$Y, simdata$X, simdata$z)
#'
#' @importFrom stats poly rgamma rnorm
flfosr <- function(Y, X, z, k = 10, S = 2000, S_burn = S/2,
                       varhyp = list(a_a = .1, b_a = .1, a_g = .1, b_g = .1, a_w = .1, b_w = .1)){


  Mi <- table(z)
  MM <- sum(Mi)

  N <- length(Mi)

  Tn <- nrow(Y)


  # tau is the obs points (length Tn)
  tau <- seq(0, 1, by = 1/(Tn-1))
  P <- psBasis(1:Tn, K = k)$P

  B = cbind(1/sqrt(Tn), poly(tau, 1),
            spikeSlabGAM::sm(tau, K = k, rankZ = .99999999,  spline.degree = 3, diff.ord = 2, centerBase = T))
  B = B/sqrt(sum(diag(crossprod(B))))
  K <- ncol(B)
  L <- ifelse(is.matrix(X) == T, ncol(X) - 1, 0)
  Dk <- diag(crossprod(B))
  Bk <- B%*%diag(1/sqrt(Dk))
  Yk <- crossprod(Y, Bk)
  group1 <- z

  Nx <- nrow(X)

  if(Nx == N){
    ZX <- rowrep(X, Mi)
  }else{
    ZX <- X
  }

  ## MCMC

  a_alph <- varhyp$a_a
  b_alph <- varhyp$b_a
  a_omega <- varhyp$a_w
  b_omega <- varhyp$b_w
  a_ga <- varhyp$a_g
  b_ga <- varhyp$b_g

  e <- 1
  alpha <- matrix(0, nrow = L+1 , K)
  ga <- (rowsum(t(Y - Bk%*%t(ZX%*%alpha)), group1)/c(Mi))%*%Bk
  w <- t(Y - Bk%*%t(ZX%*%alpha + rowrep(ga,Mi)))%*%Bk

  sig_e <- .1
  sig_alpha <- c(matrix(.1, nrow = L+1))
  sig_ga <- sum(ga^2)/(N*K)
  sig_w <- rowsum(rowSums(w^2), z)/(K*c(Mi))

  Sk <- S-S_burn

  w_post <- list()
  ga_post <- list()
  alpha_post <- list()
  sig_e_post <- rep(NA, S)
  sig_alpha_post <- matrix(NA, nrow = S, ncol = L+1)
  sig_ga_post <- matrix(NA, nrow = S, ncol = N)
  sig_w_post <- matrix(NA, nrow = S, ncol = MM)

  nog <- which(Mi == 1)

  progress <- floor(seq(1, S, length.out= 11))

  for(s in 1:S){

    eG <- matrix(rep(1/(sig_e + sig_w), K), ncol = K)

    edG <- matrix(1/(sig_e + sig_w), nrow=N, ncol=K, byrow=FALSE)

    W <- edG^2/(c(Mi)*edG + 1/sig_ga)

    alpha <- sapply(1:K, function(x) {

      XVWk <- rowrep(X*((edG - c(Mi)*W)[,x]), Mi)

      Q_alpha <- diag(1/sig_alpha) + crossprod(XVWk, ZX)
      l_alpha <- crossprod(XVWk, Yk[,x])

      ch_Q <- chol(matrix(Q_alpha, nrow=L+1, ncol=L+1))
      alpha <- backsolve(ch_Q,
                         # forwardsolve(t(ch_Q), l_alpha[,x]) +
                         forwardsolve(t(ch_Q), l_alpha) +
                           rnorm((L+1)))

      alpha

      # if(L <= Nx){
      #
      #   ch_Q <- chol(matrix(Q_alpha[,x], nrow=L+1, ncol=L+1))
      #   alpha <- backsolve(ch_Q,
      #                      forwardsolve(t(ch_Q), l_alpha[,x]) +
      #                        rnorm((L+1)))
      #
      #   alpha
      #
      # }else{
      #   H <- Map('*', 1/((1/sig_ga)+sumeG[,x]), tapply(eG[,x],  group1, function(v)
      #     if(length(v) > 1) {outer(v, v)}
      #     else{v}))
      #
      #   u <- rnorm(L+1, 0, sqrt(sig_alpha))
      #   delta <- rnorm(N)
      #
      #   Phi <- X*(sqrt(sumeG[,x] -  sapply(H, sum)))
      #   v = Phi%*%u + delta
      #
      #   pw <- solve(tcrossprod(Phi*rep(sqrt(sig_alpha), each = N)) + diag(N),
      #               (rowsum(Yk[,x],group1)/c(Mi))*(sqrt(sumeG[,x] -  sapply(H, sum))) - v)
      #   alpha <- u + sig_alpha*t(Phi)%*%pw
      #
      #
      #   alpha
      # }
    })

    Q_gaij <- 1/sig_ga + edG*c(Mi)

    l_gaij <- rowsum(rowrep(edG,Mi)*(Yk - ZX%*%alpha), group1)
    ga <- matrix(rnorm(N*K, (1/Q_gaij)*l_gaij, sqrt(1/Q_gaij)), nrow = N , K)
    ga[nog,] <- 0


    Q_wij <- rep(1/sig_w + 1/sig_e, Mi)

    if(Nx != N){
      w <- matrix(rnorm(MM*K, (1/(sig_e))*(Yk -   rowrep(X%*%alpha + ga, Mi))*c(1/Q_wij), sd = sqrt(1/Q_wij)  ), nrow = MM, ncol = K)
    }else{
      w <- matrix(rnorm(MM*K, (1/(sig_e))*(Yk -   ZX%*%alpha - rowrep(ga, Mi))*c(1/Q_wij), sd = sqrt(1/Q_wij)  ), nrow = MM, ncol = K)
    }

    sig_w <- 1/rgamma(N,  a_omega + b_omega + MM*K/2, rate = w^2/2)

    sig_ga <- 1/rgamma(1, a_ga + N*K/2, rate = b_ga + sum(((ga))^2)/2)

    sig_alpha <- 1/rgamma((L+1), a_alph + K/2, rate = b_alph + rowSums((alpha)^2)/2)

    if(MM > 50000){
      sig_e <- 1/rgamma(1, Tn*MM/2, rate = sum((Y - tcrossprod(B, (w + rowrep((X%*%alpha + ga), Mi)) ))^2)/2)
    }else{
      if(Nx == N){
        sig_e <- 1/rgamma(1, Tn*MM/2, rate = sum((Y - tcrossprod(Bk, (w +  rowrep(X%*%alpha + ga, Mi)) ))^2)/2)
      }else{
        sig_e <- 1/rgamma(1, Tn*MM/2, rate = sum((Y - tcrossprod(Bk, (w +  ZX%*%alpha + rowrep(ga, Mi)) ))^2)/2)
      }
    }

    if(s > S_burn){
      sk <- s - S_burn
      w_post[[sk]] <- w
      ga_post[[sk]] <- ga
      alpha_post[[sk]] <- alpha
      sig_e_post[sk] <- sig_e
      sig_alpha_post[sk,] <- sig_alpha
      sig_ga_post[sk] <- sig_ga
      sig_w_post[sk,] <- sig_w
    }
    if(s %in% progress){
      cat(paste0("MCMC draws: [", s, "/", S, "]\n"))
    }

  }


  #store MCMC draws of fixed effects functions
  alphaf_post <- list()
  for(i in 1:(L+1)){
    alphaf_post[[i]] <- (Bk)%*%t(do.call(rbind, lapply(alpha_post, function(x) x[i,]))[,])
  }

  m1 <- list(X = X,
             B = Bk,
             w_post = w_post,
             ga_post = ga_post,
             alpha_post = alpha_post,
             alphaf_post = alphaf_post)

  return(m1)

}
