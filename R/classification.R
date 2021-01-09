#' Base classification algorithm
#' @param xtrain A matrix with columns representing features and rows representing observations.
#' @param ytrain A binary vector specifying the corresponding labels of the observations.
#' @param xnew A matrix specifying new observations whose class labels need to be predicted.
#' @param method A character specifying the base classification method to be used. Should be one of
#' "LR" (logistic regression), "penLR" (penalized logistic regression)
#' "RF" (random forest), "GB" (gradient boosting), and "NB" (naive Bayes).
#' @param ... Additional arguments parsed to the algorithms.
#' @return A vector of classification scores for \code{xnew}.
#' @export
#' @importFrom glmnet glmnet
#' @importFrom randomForest randomForest
#' @importFrom xgboost xgboost xgb.DMatrix
#' @importFrom naivebayes naive_bayes
#' @import stats
classify_base = function(xtrain, ytrain, xnew, method, ...){
  if(method == "LR"){
    dtrain = data.frame(ytrain = ytrain, xtrain)
    obj = glm(ytrain~., family=binomial(link='logit'), data = dtrain, ...)
    if(is.null(xnew)){
      prob = predict(obj, type = "response")
    }else{
      prob = predict(obj, data.frame(xnew), type = "response")
    }
  }
  if(method == "penLR"){
    obj = glmnet(x = xtrain, y = ytrain, family = "binomial",
                 alpha = 1, ...)
    if(is.null(xnew)){
      prob = predict(obj, type = "response", s = obj$lambda)
    }else{
      prob = predict(obj, newx = xnew, type = "response", s = obj$lambda)
    }
    prob = as.numeric(prob[, 1])
  }
  if(method == "RF"){
    obj = randomForest(x = xtrain, y = factor(ytrain), ...)
    if(is.null(xnew)){
      prob = predict(obj, type = "prob")[, "1"]
    }else{
      prob = predict(obj, newdata = xnew, type = "prob")[, "1"]
    }
  }
  if(method == "GB"){
    dtrain = xgb.DMatrix(data = xtrain, label = ytrain)
    obj = xgboost(data = dtrain, nthread = 1, nrounds = 2,
                  objective = "binary:logistic", verbose = 0, ...)
    if(is.null(xnew)){
      prob = predict(obj)
    }else{
      prob = predict(obj, xnew)
    }
  }
  if(method == "NB"){
    colnames(xtrain) = colnames(xnew) = paste0("V", 1:ncol(xtrain))
    obj = naive_bayes(x = xtrain, y = as.character(ytrain))
    if(is.null(xnew)){
      prob = predict(obj, type = "prob")[,"1"]
    }else{
      prob = predict(obj, xnew, type = "prob")[,"1"]
    }
  }
  return(prob)
}

#' Stratify training data
#' @param xtrain A matrix with columns representing features and rows representing observations.
#' @param ytrain A binary vector specifying the corresponding labels of the observations.
#' @param wcost A numeric specifying the type I error cost used for stratification.
#' @return A named list specifying the stratified data. It has two elements: \code{x} and \code{y}.
#' @export
stratification = function(xtrain, ytrain, wcost){
  ind0 = which(ytrain == 0)
  ind1 = which(ytrain == 1)
  if(wcost >= 0.5){
    n1 = sum(ytrain == 1)
    new_n0 = round(n1 * wcost / (1-wcost))
    bsind = sample(ind0, size = new_n0, replace = TRUE)
    new_xtrain = xtrain[c(bsind, ind1), ]
    new_ytrain = ytrain[c(bsind, ind1)]
  }else{
    n0 = sum(ytrain == 0)
    new_n1 = round(n0 * (1-wcost) / wcost)
    bsind = sample(ind1, size = new_n1, replace = TRUE)
    new_xtrain = xtrain[c(bsind, ind0), ]
    new_ytrain = ytrain[c(bsind, ind0)]
  }
  return(list(x = new_xtrain, y = new_ytrain))
}


