# FLFOSR

`FLFOSR` (Fast Longitudinal Function-on-scalar Regression) implements a highly efficient Bayesian function-on-scalar regression model for longitudinal or repeated measurements functional data.

## Background

Consider a sample of $i = 1, \dots, N$ subjects where for each subject $i$ we have $j = 1, \dots, M_i$ functional observations $Y_{i,j}(\tau)$ along with some scalar variables ${x}_{i,j} = \{x_{i,j,1}, \dots, x_{i,j,L}\}$. Then we can write the following regression model:

$$
Y_{i,j}(\tau) = {\alpha}_0(\tau) + \sum_{\ell=1}^{L} x_{i,j,\ell}{\alpha}_{\ell}(\tau) + {\gamma}_{i}(\tau) + {\omega}_{i,j}(\tau) + \epsilon_{i,j}\left(\tau\right), \quad \epsilon_{i,j}\left(\tau\right) \overset{iid}{\sim} N(0, \sigma^2_{\epsilon}).
$$

In this model we are regressing the functional response variable $Y_{i,j}(\tau)$ against the scalar predictors ${x}_{i,j}$. Here the regression coefficient ${\alpha}_{\ell}(\tau)$ represents the expected increase in $Y_{i,j}$ at time $\tau$ given a one unit increase in $x_{i,j,\ell}$.

Additionally, we have the subject-specific functional random intercepts ${\gamma}_{i}(\tau)$ which represent the subject-level mean deviation from the overall functional mean (global intercept) ${\alpha}_0(\tau)$. Similarly, the curve-specific functional random intercepts ${\omega}_{i,j}(\tau)$ represent curve-level deviations. Finally, the error term $\epsilon_{i,j}\left(\tau\right)$ accounts for noisy errors.

**Purpose**: `FLFOSR` efficiently estimates all regression terms in a fully Bayesian model using a highly optimized MCMC algorithm. This method has been demonstrated to achieve high point accuracy and accurate interval coverage while having better scalability than even Variational Bayes or Frequentist methods as $N$, $M_i$, or $L$ increases. See <https://arxiv.org/abs/2306.07168> for details.


## Installation

You can install and load the package using the following lines:
```{r install}
# install.packages("devtools")
# devtools::install_github("thomasysun/FLFOSR")
library(FLFOSR) 
```

## Using `FLFOSR`


The primary function in `FLFOSR` is `flfosr()`, which runs the model. It requires three inputs:

- `Y` : the $T \times M$ matrix of functional observations, where $T$ is the total number of observed points along the domain and $M$ is the total number of curves, $M = \sum^{N}_{i=1} M_i$.
- `X` : the matrix of scalar observations, which can be either $N \times L+1$ or $M \times L+1$. It should include a column of ones for the intercept. 
- `z` : vector specifying group memberships of each curve. For example, if there are 20 subjects with 5 repeated measurements each, then `z=rep(1:20, each = 5)`.


