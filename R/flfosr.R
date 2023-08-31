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
#' * X: the inputted design matrix
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
  ga <- (rowsum(t(Y - Bk%*%t(ZX%*%alpha)), z)/c(Mi))%*%Bk
  w <- t(Y - Bk%*%t(ZX%*%alpha + rowrep(ga,Mi)))%*%Bk

  sig_e <- .1
  sig_alpha <- c(matrix(.1, nrow = L+1))
  sig_ga <- sum(ga^2)/(N*K)
  sig_w <- rowsum(rowSums(w^2), z)/(K*c(Mi))

  Sk <- S-S_burn

  w_post <- list()
  ga_post <- list()
  alpha_post <- list()
  sig_e_post <- rep(NA, Sk)
  sig_alpha_post <- matrix(NA, nrow = Sk, ncol = L+1)
  sig_ga_post <- matrix(NA, nrow = Sk, ncol = 1)
  sig_w_post <- matrix(NA, nrow = Sk, ncol = N)

  progress <- floor(seq(1, S, length.out= 11))

  cat(paste0("Beginning MCMC with ", Sk, " draws after burn-in\n"))
  for(s in 1:S){
    edG <- sapply(Dk, function(x) 1/(sig_e + (x)*sig_w))

    edGdF <- sapply(Dk, function(x) 1/(x*sig_ga*Mi + x*c(sig_w) + sig_e))

    eGh <- rep(Dk, each = MM)*rowrep(edGdF, Mi)

    l_alpha <- crossprod((ZX), Yk*eGh/rep(Dk, each = MM))

    alpha <- sapply(1:K, function(x) {

      if(L <= Nx){

        Q_alpha <- diag(1/sig_alpha) + crossprod(ZX, ZX*eGh[,x])

        ch_Q <- chol(matrix(Q_alpha, nrow=L+1, ncol=L+1))
        alpha <- backsolve(ch_Q,
                           forwardsolve(t(ch_Q), l_alpha[,x]) +
                             rnorm((L+1)))

        alpha

      }else{
        u <- rnorm(L+1, 0, sqrt(sig_alpha))
        delta <- rnorm(MM)

        Phi <- ZX*sqrt(eGh[,x])
        v = Phi%*%u + delta

        pw <- solve(tcrossprod(Phi*rep(sqrt(sig_alpha), each = MM)) + diag(MM),
                    Yk[,x]*sqrt(eGh[,x])/(Dk[x]) - v)
        alpha <- u + sig_alpha*t(Phi)%*%pw


        alpha
      }
    })
    alpha

    Q_gaij <- 1/sig_ga + rep(Dk, each = N)*edG*c(Mi)

    l_gaij <- rowsum(rowrep(edG,Mi)*(Yk - ZX%*%alpha*rep(Dk, each = MM)), z)

    ga <- matrix(rnorm(N*K, (1/Q_gaij)*l_gaij, sqrt(1/Q_gaij)), nrow = N , K)


    Q_wij <- rowrep(sapply(Dk, function(x) x/sig_e + 1/sig_w), Mi)

    if(Nx != N){
      w <- matrix(rnorm(MM*K, (1/sig_e)*(Yk -   rowrep(X%*%alpha + ga, Mi)*rep(Dk, each = MM))*(1/Q_wij), sd = sqrt(1/Q_wij)  ), nrow = MM, ncol = K)
    }else{
      w <- matrix(rnorm(MM*K, (1/sig_e)*(Yk -   (ZX%*%alpha + rowrep(ga, Mi))*rep(Dk, each = MM))*(1/Q_wij), sd = sqrt(1/Q_wij)  ), nrow = MM, ncol = K)
    }

    sig_w <- 1/rep(rgamma(1, a_omega + MM*K/2, rate = b_omega + w^2/2), N)

    sig_ga <- 1/rgamma(1, a_ga + N*K/2, rate = b_ga + sum(((ga))^2)/2)

    sig_alpha <- 1/rgamma((L+1), a_alph + K/2, rate = b_alph + rowSums((alpha)^2)/2)

    if(MM > 50000){
      sig_e <- 1/rgamma(1, Tn*MM/2, rate = sum((Y - tcrossprod(B, (w + rowrep((X%*%alpha + ga), Mi)) ))^2)/2)
    }else{
      if(Nx == N){
        sig_e <- 1/rgamma(1, Tn*MM/2, rate = sum((Y - tcrossprod(B, (w +  rowrep(X%*%alpha + ga, Mi)) ))^2)/2)
      }else{
        sig_e <- 1/rgamma(1, Tn*MM/2, rate = sum((Y - tcrossprod(B, (w +  ZX%*%alpha + rowrep(ga, Mi)) ))^2)/2)
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
    alphaf_post[[i]] <- (B)%*%t(do.call(rbind, lapply(alpha_post, function(x) x[i,]))[,])
  }

  m1 <- list(X = X,
             B = B,
             w_post = w_post,
             ga_post = ga_post,
             alpha_post = alpha_post,
             sig_e_post = sig_e_post,
             sig_alpha_post = sig_alpha_post,
             sig_ga_post = sig_ga_post,
             sig_w_post = sig_w_post,
             alphaf_post = alphaf_post)

  return(m1)

}
