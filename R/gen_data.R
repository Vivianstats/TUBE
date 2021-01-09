#' Generate test and example data
#' @param model A character specifying the model used to generate data.
#' Should be one of "gaussian", "t", or "mixture".
#' @param seed An integer specifying the seed (optional).
#' @param n An integer specifying the sample size.
#' @param d An intger specifying the feature number/dimension.
#' @param pi A numeric specifying the proportion of class 0 sample. Defaults to 0.5.
#' @param nleave An integer specifying the size of left-out sample. For debug only.
#' @return A named list of generated data. It has two elements: \code{x} and \code{y}.
#' @export
#' @import MASS
#' @import mvtnorm
gen_data = function(model, seed = NULL, n, d, pi = 0.5, nleave = 0) {
  if(!is.null(seed)) {set.seed(seed)}
  modellist = c("gaussian", "t", "mixture")
  if(! (model %in% modellist)) stop("Unrecognized data model!")

  if(model == "gaussian"){
    mu0 = c(c(0,0), rep(0, d-2)); mu1 = c(c(1.5,1.5), rep(0, d-2))
    cov = diag(d)
    cov[abs(row(cov) - col(cov)) == 1] = 0.5
    sigma1 = cov
    sigma0 = diag(d)
    y = sample(0:1, size=n, replace=TRUE, prob=c(pi, 1-pi))
    n1 = sum(y)
    n0 = n - n1 + nleave
    x = rep(0, n)
    if(d != length(mu0) | d != nrow(sigma0)) stop("wrong dimension!")
    S0 = list(x=mvrnorm(n=n0, mu0, sigma0), y=rep(0, n0))
    S1 = list(x=mvrnorm(n=n1, mu1, sigma1), y=rep(1, n1))
  }
  if(model %in% c("t")){
    if(d >= 3){
      mu0 = c(0,0); mu1 = c(2.5,2.5)

      y = sample(0:1, size=n, replace=TRUE, prob=c(pi, 1-pi))
      n1 = sum(y)
      n0 = n - n1 + nleave

      x01 = rmvt(n0, df = 3, delta = mu0, type = "shifted")
      x02 = matrix(rnorm(n0*(d-2)), nrow = n0)
      x0 = cbind(x01, x02)
      x11 = rmvt(n1, df = 3, delta = mu1, type = "shifted")
      x12 = matrix(rnorm(n1*(d-2)), nrow = n1)
      x1 = cbind(x11, x12)
      S0 = list(x=x0, y=rep(0, n0))
      S1 = list(x=x1, y=rep(1, n1))
    }else{stop("d<3 in t distribution!")}
  }

  if(model == "mixture"){
    a = 2/sqrt(d)
    mu01 = rep(a, d); mu02 = rep(-a, d)
    mu1 = rep_len(c(a, -a), length.out = d)

    nn = rmultinom(1,n,c(pi,1-pi))
    n0 = rmultinom(1, nn[1,1] + nleave, c(1/2,1/2))
    n1 = nn[2,1]

    X1 = rbind(t(matrix(mu01,d,n0[1])), t(matrix(mu02,d,n0[2]))) + matrix(rnorm(sum(n0)*d),sum(n0),d)
    X2 = t(matrix(mu1,d,n1)) + matrix(rnorm(n1*d),n1,d)
    S0 = list(x=X1, y=rep(0,sum(n0)))
    S1 = list(x=X2, y=rep(1,n1))
  }
  if(nleave > 0){
    xtrain0 = S0$x
    ind_leave = sample(1:nrow(xtrain0), nleave)
    xleave = xtrain0[ind_leave, , drop = FALSE]
    x = rbind(xtrain0[-ind_leave, , drop = FALSE], S1$x)
    y = c(S0$y[-ind_leave], S1$y)
  }else{
    xleave = NULL
    x = rbind(S0$x, S1$x)
    y = c(S0$y, S1$y)
  }
  return(list(x=x, y=y, xleave = xleave))
}


