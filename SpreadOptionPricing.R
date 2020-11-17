

# Initialization ----------------------------------------------------------
if (1) {
    # 
    setwd('~/Documents/0_ongoing/fe5222_project2/')
    # 
    library(plyr)
    library(tidyverse)
    library(doParallel)
    library(ggplot2)
    # 
    set.seed(2020)
}


Spread_Call <- function (
    S_0_1 = 100, S_0_2 = 90, sigma_1 = .2, sigma_2 = .2, rho = .4, 
    Strike = 1, rf = .08, mT = 1, 
    flag_AV = F, n = 1000000, seed = 2020
) {
    #' Spread Call Monte Carlo Pricing
    #'
    #' @param S_0_1     Spot price of Stock 1
    #' @param S_0_2     Spot price of Stock 2
    #' @param sigma_1   Vol of Stock 1
    #' @param sigma_2   Vol of Stock 2
    #' @param rho       instantaneous Correlation
    #' @param Strike 
    #' @param rf        risk free rate
    #' @param mT        Time to Maturity
    #' @param flag_AV   indicator to use Antithetic Variates
    #' @param n         Number of paths
    #' @param seed      random seed to reproduce
    #'
    set.seed(seed)
    VarCov <- matrix(c(sigma_1^2, rep(rho * sigma_1 * sigma_2, 2), sigma_2^2),
                     nrow = 2, byrow = T)
    if (flag_AV) {
        temp <- MASS::mvrnorm(n = round(n/2), 
                              mu = rep(0, 2), Sigma = VarCov)
        temp <- rbind(temp, - temp)
    } else {
        temp <- MASS::mvrnorm(n = n, 
                              mu = rep(0, 2), Sigma = VarCov)
    }
    # 
    seq_S_T_1 <- S_0_1 * exp(
        (rf - .5 * sigma_1 ^ 2) * mT + sigma_1 * sqrt(mT) * temp[, 1]
    )
    seq_S_T_2 <- S_0_2 * exp(
        (rf - .5 * sigma_2 ^ 2) * mT + sigma_2 * sqrt(mT) * temp[, 2]
    )
    # 
    seq_Payoff <- pmax(seq_S_T_1 - seq_S_T_2 - Strike, 0)
    Value <- mean(seq_Payoff) * exp(- rf * mT)
    return(Value)
}


Spread_Call_Kirk <- function (
    S_0_1 = 100, S_0_2 = 90, sigma_1 = .2, sigma_2 = .2, rho = .4, 
    Strike = 1, rf = .08, mT = 1
) {
    temp <- S_0_2/(S_0_2 + Strike * exp(- rf * mT))
    sigma <- sqrt(
        sigma_1^2 + sigma_2^2 * temp^2 - 2 * rho * sigma_1 * sigma_2 * temp
    )
    temp <- temp / S_0_2 * S_0_1
    d_1 <- log(temp) / (sigma * sqrt(mT)) + .5 * sigma * sqrt(mT)
    d_2 <- d_1 - sigma * sqrt(mT)
    temp <- S_0_1 / temp
    Value <- S_0_1 * pnorm(d_1) - temp * pnorm(d_2)
    return(Value)
}


Spread_Call_Kirk()
fExoticOptions::SpreadApproxOption(TypeFlag = 'c', S1 = 100, S2 = 90, X = 1, 
                                   Time = 1, r = .08, sigma1 = .2, sigma2 = .2, rho = .4)
Spread_Call(n = 10000)

