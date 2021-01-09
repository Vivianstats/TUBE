
#' TUBEc-assisted cost-sensitive classification
#' @param data A named list specifying the training data. It should have two elements:
#' \code{data$x} should be a matrix with columns representing features and rows representing observations.
#' \code{data$y} should be a binary vector specifying the corresponding labels of the observations.
#' @param classify_fun A function specifying the cost-sensitive classification algorithm.
#' @param alpha A numeric specifying the target type I error (between 0 and 1). Defaults to 0.05.
#' @param delta A numeric specifying the target violation rate (between 0 and 1). Defaults to 0.1.
#' @param t_cs A numeric specifying the threshold applied to classification scores. Defaults to 0.5.
#' @param nleave A integer specifying the sample size of left-out class 0 data. Should be smaller than the
#' number of class 0 cases in \code{data}.
#' @param cost A numeric vector specifying the candidate type I error costs. Should be ordered increasingly.
#' Defaults to \code{seq(0.1, 0.99, 0.01)}.
#' @param B An integer specifying the number of bootstrap samples. Defaults to 200.
#' @return A list with two elements. \code{c0} gives the selected type I error cost;
#' \code{yhat} gives the predicted class labels for \code{data$x} given \code{c0}.
#' @export
#' @importFrom stats quantile
#' @author Wei Vivian Li, \email{vivian.li@rutgers.edu}
#' @author Xin Tong, \email{xtong001@gmail.com>}
#' @author Jingyi Jessica Li, \email{jli@stat.ucla.edu}
#' @references Li WV, Tong X, Li JJ. Bridging Cost-sensitive and Neyman-Pearson Paradigms for Asymmetric Binary Classification.
#' arXiv preprint arXiv:2012.14951. 2020 Dec 29. \url{https://arxiv.org/abs/2012.14951}
TUBEc = function(data, classify_fun, alpha = 0.05, delta = 0.1,
                 t_cs = 0.5, nleave = NULL,
                 cost = seq(0.1, 0.99, 0.01), B = 200){
  XX = data$x
  YY = data$y
  x_c0 = XX[YY == 0, ]
  x_c1 = XX[YY == 1, ]
  if(is.null(nleave)){
    nleave = floor(sum(YY == 0)/2)
  }
  ind_leave = sample(1:nrow(x_c0), nleave, replace = FALSE)
  xleave = x_c0[ind_leave, ]
  xtrain = rbind(x_c0[-ind_leave, ], x_c1)
  ytrain = c(rep(0, nrow(x_c0)-nleave), rep(1, nrow(x_c1)))

  c = 1
  wcost = cost[c]
  alp_emp = 1
  while(alp_emp > alpha){
    if(wcost < max(cost)){
      c = c+1; wcost = cost[c]
    }else{break}
    prob1 = classify_fun(xtrain, ytrain, cost = wcost, xnew = xleave)
    alp_bs = sapply(1:B, function(b2){
      set.seed(b2)
      probb = prob1[sample(1:nleave, nleave, replace = TRUE)]
      # find k_cs
      if(sum(probb <= t_cs) == 0){
        alp = 1
      }else{
        k_cs = sum(probb <= t_cs, na.rm = TRUE)
        u = 1-(k_cs)/nleave
        alp = 1+u-(delta+u^nleave)^(1/nleave)
        alp = max(alp, u)
      }
      return(alp)
    })
    alp_emp = quantile(alp_bs, 1-delta)
  }
  w0 = wcost

  prob1 = classify_fun(xtrain, ytrain, cost = w0, xnew = XX)

  yhat = as.numeric(prob1 > t_cs)
  return(list(c0 = w0, yhat = yhat))
}


#' TUBE-assisted cost-sensitive classification
#' @param data A named list specifying the training data. It should have two elements:
#' \code{data$x} should be a matrix with columns representing features and rows representing observations.
#' \code{data$y} should be a binary vector specifying the corresponding labels of the observations.
#' @param classify_fun A function specifying the cost-sensitive classification algorithm.
#' @param alpha A numeric specifying the target type I error (between 0 and 1). Defaults to 0.05.
#' @param delta A numeric specifying the target violation rate (between 0 and 1). Defaults to 0.1.
#' @param t_cs A numeric specifying the threshold applied to classification scores. Defaults to 0.5.
#' @param nleave A integer specifying the sample size of left-out class 0 data. Should be smaller than the
#' number of class 0 cases in \code{data}.
#' @param cost A numeric vector specifying the candidate type I error costs. Should be ordered increasingly.
#' Defaults to \code{seq(0.1, 0.99, 0.01)}.
#' @param B1 An integer specifying the number of random data splitting. Defaults to 50.
#' @param B2 An integer specifying the number of bootstrap samples. Defaults to 200.
#' @return A list with two elements. \code{c0} gives the selected type I error cost;
#' \code{yhat} gives the predicted class labels for \code{data$x} given \code{c0}.
#' @export
#' @importFrom stats quantile
#' @author Wei Vivian Li, \email{vivian.li@rutgers.edu}
#' @author Xin Tong, \email{xtong001@gmail.com>}
#' @author Jingyi Jessica Li, \email{jli@stat.ucla.edu}
#' @references Li WV, Tong X, Li JJ. Bridging Cost-sensitive and Neyman-Pearson Paradigms for Asymmetric Binary Classification.
#' arXiv preprint arXiv:2012.14951. 2020 Dec 29. \url{https://arxiv.org/abs/2012.14951}
TUBE = function(data, classify_fun, alpha = 0.05, delta = 0.1,
                 t_cs = 0.5, nleave = NULL,
                 cost = seq(0.1, 0.99, 0.01), B1 = 50, B2 = 200){
  xobs = data$x
  yobs = data$y

  if(is.null(nleave)){
    nleave = floor(sum(yobs == 0)/2)
  }

  c = 1
  wcost = cost[c]
  alp_emp = 1
  while(alp_emp > alpha){
    if(wcost < max(cost)){
      c = c+1; wcost = cost[c]
    }else{break}
    prob_train = classify_fun(xobs, yobs, cost = wcost, xnew = xobs)
    est_all_emp = sum(prob_train > t_cs & yobs == 0, na.rm = TRUE)/sum(!is.na(prob_train))

    exps = sapply(1:B1, function(b1){
      set.seed(b1)

      ind0 = which(yobs == 0)
      ind_leave = sample(ind0, size = nleave)
      xtrain = xobs[-ind_leave, , drop = FALSE]
      ytrain = yobs[-ind_leave]
      xleave = xobs[ind_leave, , drop = FALSE]

      # leave-out sample scores
      xmerge = rbind(xleave, xtrain)
      prob_merge = classify_fun(xtrain, ytrain, cost = wcost, xnew = xmerge)
      prob1 = prob_merge[1:nleave]

      # Bootstrap
      alp_bs = sapply(1:B2, function(b2){
        probb = prob1[sample(1:nleave, nleave, replace = TRUE)]
        # find k_cs
        if(sum(probb <= t_cs) == 0){
          alp = 1
        }else{
          k_cs = sum(probb <= t_cs, na.rm = TRUE)
          u = 1-(k_cs)/nleave
          alp = 1+u-(delta+u^nleave)^(1/nleave)
          alp = max(alp, u)
        }
        return(alp)
      })

      prob_train = prob_merge[-(1:nleave)]
      est_split_emp = sum(prob_train > t_cs & ytrain == 0, na.rm = TRUE)/sum(!is.na(prob_train))

      est_split_tube = quantile(alp_bs, 1-delta)
      est_all_tube =  est_split_tube - (est_split_emp - est_all_emp)
      return(est_all_tube)
    })
    alp_emp = quantile(exps, delta)
  }
  w0 = wcost

  prob1 = classify_fun(xobs, yobs, cost = w0, xnew = xobs)

  yhat = as.numeric(prob1 > t_cs)
  return(list(c0 = w0, yhat = yhat))
}



