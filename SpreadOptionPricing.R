

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


# Spread Call & Kirk's Approximation --------------------------------------

Spread_Call <- function (
    S_0_1 = 100, S_0_2 = 90, sigma_1 = .2, sigma_2 = .2, rho = .4, 
    Strike = 1, rf = .08, mT = 1, 
    flag_AV = F, n = 10000, seed = 2020
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
    # Generate for the terminal only, since the payoff is not path-dependent
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
    # stddev modification is embedded in the above VarCov
    seq_S_T_1 <- S_0_1 * exp(
        (rf - .5 * sigma_1 ^ 2) * mT + sqrt(mT) * temp[, 1] # * sigma_1
    )
    seq_S_T_2 <- S_0_2 * exp(
        (rf - .5 * sigma_2 ^ 2) * mT + sqrt(mT) * temp[, 2] # * sigma_1
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

if (1) {
    # Histogram example
    seq_temp <- seq_len(100000) %>% sapply(
        FUN = function(iter_temp){Spread_Call(seed = iter_temp, 
                                              flag_AV = T, n = 100)})
    hist(seq_temp, breaks = 17, probability = T, main = 'MC Hist', xlab = 'Spread Call MC ($)')
    lines(density(seq_temp, bw = 'SJ'), col = 'red', lwd = 1.5)
    abline(v = Spread_Call_Kirk(), col = 'blue', lwd = 3, lty = 2)
    legend('topright', legend = c('Kirk\'s Approximation' , 'KDE'), 
           col = c('blue', 'red'), lty = c(2, 1), lwd = 2)
}



# Grid Comparison for pricing parameters ----------------------------------

if (1) {
    temp_seq <- paste0('seq_', c('Spot_Diff', 'Strike', 'sigma', 'r', 'mT'))
    temp_seq %>% print()
    combn(length(temp_seq), 2)
    # pricing parameters
    seq_Spot <- 1:20 * 10
    seq_Strike <- 1:20 * 10
    seq_sigma <- 1:20 / 40
    seq_r <- 1:20 / 200
    seq_mT <- 1:20 / 10
}


