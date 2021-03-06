Bridging Cost-sensitive and Neyman-Pearson Paradigms for Asymmetric Binary Classification
================

<!-- README.md is generated from README.Rmd. Please edit that file -->
## Latest News

> 2021/01/08: Version 0.1.0 released!

## Introduction

To bridge cost-sensitive and Neyman-Pearson paradigms for asymmetric binary classification, we have developed two algorithms, TUBEc and TUBE (estimating the Type I error Upper Bound of a cost-sensitive classifiEr). This package contains functions of TUBE-assisted and TUBEc-assisted cost-sensitive classification methods for selecting the type I error cost so that the resulting cost-sensitive classifier has its population type I error under a target upper bound with high probability.

Any suggestions on the package are welcome! For technical problems, please report to [Issues](https://github.com/Vivianstats/TUBE/issues). For suggestions and comments on the package, please contact Vivian (<vivian.li@rutgers.edu>).

A detailed introduction to the methods and algorithms implemented in this package is available in our [manuscript](https://arxiv.org/abs/2012.14951).

## Installation

You can install from Github with:

``` r
# install.packages("devtools")
devtools::install_github("Vivianstats/TUBE")
```

## Quick Start

The packages contains two main functions, `TUBE` and `TUBEc`, corresponding to the TUBE-assisted or TUBEc-assisted cost-sensitive classification paradigms decribed in the manuscript.

To use either function, the first step is to set up a cost-sensitive classification function. We provide convenient ways to apply the stratification approach to commonly used classification algorithms. For example, if the user wants to use cost-sensitive logistic regression, then the function can be set up as follows

``` r
library(TUBE)

classify_fun = function(xtrain, ytrain, cost, xnew, method = "LR", ...){
  data_str = stratification(xtrain, ytrain, cost)
  score = classify_base(xtrain = data_str$x, ytrain = data_str$y, xnew, method, ...)
  return(score)
}
```

The `method` argument can be one "LR" (logistic regression), "penLR" (penalized logistic regression) "RF" (random forest), "GB" (gradient boosting), and "NB" (naive Bayes). If the user wants to use a cost-sensitive approach other than stratification or a base classification algorithm not belonging to the above list, they can define their own `classify_fun` function that takes `xtrain` (feature matrix of training data), `ytrain` (class labels of training data), `cost` (type I error cost), `xnew` (feature matrix of new data) as arguments and returns the classification scores for .

As an example, we use the `gen_data` function to generate data based on Gaussian distributions (as described in the manuscript).

``` r
set.seed(1)
data = gen_data(model = "gaussian",  # could also be "t" or "mixture"
               n = 1000, # sample size
               d = 30, # feature dimension
               pi = 0.5 # proportion of class 0 data
               )
dim(data$x)
#> [1] 1000   30
length(data$y)
#> [1] 1000
```

Then, we can use the `TUBE` or `TUBEc` function to select the type I error cost based on a target type I error.

``` r
res_tube = TUBE(data, # training data
                classify_fun, # cost-sensitive
                alpha = 0.1, # target type I error 
                delta = 0.1, # target violation rate
                t_cs = 0.5, # threshold on classification scores
                cost = seq(0.1, 0.99, 0.01) # candidate type I error costs
                )
                
res_tubec = TUBEc(data, classify_fun, alpha = 0.1, delta = 0.1, t_cs = 0.5, 
                  cost = seq(0.1, 0.99, 0.01))

res_tube$c0
#> [1] 0.72

res_tubec$c0
#> [1] 0.76
```

For detailed usage, please refer to the package [manual](https://github.com/Vivianstats/TUBE/blob/master/inst/docs/).
