

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
    S_0_1 = 110, S_0_2 = 100, sigma_1 = .2, sigma_2 = .2, rho = .4, 
    Strike = 5, rf = .08, mT = 1, 
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
    S_0_1 = 110, S_0_2 = 100, sigma_1 = .2, sigma_2 = .2, rho = .4, 
    Strike = 5, rf = .08, mT = 1
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

if (0) {
    # Histogram example
    seq_temp <- seq_len(10000) %>% sapply(
        FUN = function(iter_temp){Spread_Call(seed = iter_temp, 
                                              flag_AV = T, n = 10000)})
    hist(seq_temp, breaks = 17, probability = T, 
         main = 'MC Hist', xlab = 'Spread Call MC ($)')
    lines(density(seq_temp, bw = 'SJ'), col = 'red', lwd = 1.5)
    abline(v = Spread_Call_Kirk(), col = 'blue', lwd = 3, lty = 2)
    legend('topright', legend = c('Kirk\'s Approximation' , 'KDE'), 
           col = c('blue', 'red'), lty = c(2, 1), lwd = 2)
}



# Grid Comparison for pricing parameters ----------------------------------

if (1) {
    head(Spread_Call)
    head(Spread_Call_Kirk)
    # Spot_Diff, sigma_1, sigma_2, rho, StrikePercent of SpotDiff
    temp_seq <- c('SpotDiff', 'sigma_1', 'sigma_2', 'rho', 'StrikePercent')
    combn(temp_seq, 2)
    # pricing parameters
    seq_SpotDiff <- seq(from = 0, to = 200, by = 10)
    seq_sigma_1 <- seq(from = 0, to = 1, by = .05)
    seq_sigma_2 <- seq(from = 0, to = 1, by = .05)
    seq_rho <- seq(from = -1, to = 1, by = .1)
    seq_StrikePercent <- seq(from = 0, to = 2, by = .1)
    # Strike as in percentages of Spot_Diff
}

ContourPlots <- function(
    temp_df,
    # iter_X, iter_Y, iter_Z, 
    fill_label = 'Value Diff\n(Kirk\'s - MC, in $)',
    xlab = 'X', ylab = 'Y', fig_title = '', fig_subtitle = '',
    scale_range = NULL,
    flag_X = F,
    flag_plot = F, flag_save_plot = T, save_folder = ''
) {
    #' Title
    #' @param iter_X                        x seq
    #' @param iter_Y                        y seq
    #' @param iter_Z                        z matrix
    #' @param fill_label                    legend title
    # Need to plot
    require(ggplot2)
    lut_16 <- matrix(
        data = c(0,0,0,    1,1,171,    1,1,224,    0,110,255,    1,171,254,
                 1,224,254,    1,254,1,    190,255,0,    255,255,0,    255,224,0,
                 255,141,0,    250,94,0,    245,0,0,    245,0,172,    222,180,222),
        ncol = 3, byrow = T)
    lut_16 <- rgb(lut_16, maxColorValue = 255)
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # temp_df <- expand.grid(X = iter_X, Y = iter_Y)    
    # temp_df$Z <- c(iter_Z)
    # temp_df <- temp_grid
    seq_X <- sort(unique(temp_df$X))
    seq_Y <- sort(unique(temp_df$Y))
    temp_df$X <- (temp_df$X - min(seq_X))/(max(seq_X) - min(seq_X))
    temp_df$Y <- (temp_df$Y - min(seq_Y))/(max(seq_Y) - min(seq_Y))
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    temp_p <- ggplot(data = temp_df, aes(x = X, y = Y, z = Z)) + 
        geom_raster(data = temp_df, aes(fill = Z)) +
        coord_fixed() +
        ggplot2::labs(title = fig_title, subtitle = fig_subtitle, x = xlab, y = ylab, fill = fill_label) +
        # ggplot2::labs(title = 'fig_title', subtitle = 'fig_subtitle', x = 'xlab', y = 'ylab', fill = 'fill_label') +
        ggplot2::theme_minimal()
    
    # ggplot2::scale_y_discrete(labels = unique(seq_Y)/max(seq_Y), labels = unique(seq_Y))
    # ggplot2::theme_bw()
    if (flag_X) {
        temp_seq_X <- seq_X[seq(from = 1, to = length(seq_X), by = 2)]
        temp_p <- temp_p + 
            ggplot2::scale_x_continuous(
                breaks = (temp_seq_X - min(temp_seq_X))/(max(temp_seq_X) - min(temp_seq_X)),
                labels = temp_seq_X) + 
            ggplot2::scale_y_continuous(
                breaks = (seq_Y - min(seq_Y))/(max(seq_Y) - min(seq_Y)), labels = seq_Y)
    } else {
        temp_p <- temp_p + 
            ggplot2::scale_x_continuous(
                breaks = (seq_X - min(seq_X))/(max(seq_X) - min(seq_X)), 
                labels = seq_X) + 
            ggplot2::scale_y_continuous(
                breaks = (seq_Y - min(seq_Y))/(max(seq_Y) - min(seq_Y)), labels = seq_Y)
    }
    
    if (is.null(scale_range)) {
        temp_p <- temp_p +
            ggplot2::scale_fill_gradient2(
                low = 'blue', high = 'red', mid = 'white', midpoint = 0,
                # colours = colorRamps::matlab.like(300)
                guide = ggplot2::guide_colourbar(
                    raster = T, frame.colour = "black", frame.linewidth = 1
                ), na.value = 'white'
            )
        # ggplot2::scale_fill_gradientn(
        #     colours = lut_16,
        #     # colours = colorRamps::matlab.like(300)
        #     guide = ggplot2::guide_colourbar(
        #         raster = T, frame.colour = "black", frame.linewidth = 1
        #     ), na.value = 'white'
        # )
    }else if (length(scale_range) == 2) {
        temp_p <- temp_p +
            ggplot2::scale_fill_gradient2(
                limits = c(scale_range[1], scale_range[2]),
                low = 'blue', high = 'red', mid = 'white', midpoint = 0,
                # colours = lut_16,
                # colours = colorRamps::matlab.like(300)
                guide = ggplot2::guide_colourbar(
                    raster = T, frame.colour = "black", frame.linewidth = 1
                ), na.value = 'white'
            )
        # ggplot2::scale_fill_gradientn(
        #     limits = c(scale_range[1], scale_range[2]),
        #     colours = lut_16,
        #     # colours = colorRamps::matlab.like(300)
        #     guide = ggplot2::guide_colourbar(
        #         raster = T, frame.colour = "black", frame.linewidth = 1
        #     ), na.value = 'white'
        # )
    }
    
    if (flag_save_plot) {
        cat('Saving', paste0(save_folder, fig_title, '.png'), '...\n')
        ggplot2::ggsave(filename = paste0(save_folder, fig_title, '.png'), plot = temp_p)
        # , width = 7, height = 6.5)
    }
    if (flag_plot) {print(temp_p)}
    return(temp_p)
}


# 1 SpotDiff & sigma_1 ----------------------------------------------------
if (1) {
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    tictoc::tic()
    set.seed(2020)
    # v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v 
    temp_grid <- cbind(expand.grid(seq_SpotDiff, seq_sigma_1), NA)
    # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    colnames(temp_grid) <- c('X', 'Y', 'Z')
    for (iter_i in 1:nrow(temp_grid)) {
        temp_grid$Z[iter_i] <- 
            # v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v 
            Spread_Call_Kirk(S_0_1 = 100 + temp_grid$X[iter_i], 
                             sigma_1 = temp_grid$Y[iter_i]) - 
            Spread_Call(S_0_1 = 100 + temp_grid$X[iter_i], 
                        sigma_1 = temp_grid$Y[iter_i],
                        flag_AV = T, n = 10^6) 
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    }
    p <- ContourPlots(
        temp_grid, 
        fig_subtitle = 'Differences between Kirk & MC for Spread Call',
        scale_range = c(-.6, .6), flag_X = T,
        flag_plot = F, flag_save_plot = T, 
        save_folder = '~/Documents/0_ongoing/fe5222_project2/plots/',
        # v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v 
        xlab = 'Spot Diff ($)', ylab = 'sigma_1', 
        fig_title = paste0('Spot Diff', ' & ', 'sigma_1')
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    )
    print(p)
    print(range(temp_grid$Z))
    tictoc::toc()
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
}



# 2 SpotDiff & sigma_2 ----------------------------------------------------
if (1) {
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    tictoc::tic()
    set.seed(2020)
    # v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v 
    temp_grid <- cbind(expand.grid(seq_SpotDiff, seq_sigma_2), NA)
    # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    colnames(temp_grid) <- c('X', 'Y', 'Z')
    for (iter_i in 1:nrow(temp_grid)) {
        temp_grid$Z[iter_i] <- 
            # v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v 
            Spread_Call_Kirk(S_0_1 = 100 + temp_grid$X[iter_i], 
                             sigma_2 = temp_grid$Y[iter_i]) - 
            Spread_Call(S_0_1 = 100 + temp_grid$X[iter_i], 
                        sigma_2 = temp_grid$Y[iter_i],
                        flag_AV = T, n = 10^6) 
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    }
    p <- ContourPlots(
        temp_grid, 
        fig_subtitle = 'Differences between Kirk & MC for Spread Call',
        scale_range = c(-.6, .6), flag_X = T,
        flag_plot = F, flag_save_plot = T, 
        save_folder = '~/Documents/0_ongoing/fe5222_project2/plots/',
        # v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v 
        xlab = 'Spot Diff ($)', ylab = 'sigma_2', 
        fig_title = paste0('Spot Diff', ' & ', 'sigma_2')
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    )
    print(p)
    print(range(temp_grid$Z))
    tictoc::toc()
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
}



# 3 SpotDiff & rho --------------------------------------------------------
if (1) {
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    tictoc::tic()
    set.seed(2020)
    # v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v 
    temp_grid <- cbind(expand.grid(seq_SpotDiff, seq_rho), NA)
    # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    colnames(temp_grid) <- c('X', 'Y', 'Z')
    for (iter_i in 1:nrow(temp_grid)) {
        temp_grid$Z[iter_i] <- 
            # v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v 
            Spread_Call_Kirk(S_0_1 = 100 + temp_grid$X[iter_i], 
                             rho = temp_grid$Y[iter_i]) - 
            Spread_Call(S_0_1 = 100 + temp_grid$X[iter_i], 
                        rho = temp_grid$Y[iter_i],
                        flag_AV = T, n = 10^6) 
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    }
    p <- ContourPlots(
        temp_grid, 
        fig_subtitle = 'Differences between Kirk & MC for Spread Call',
        scale_range = c(-.1, .1), flag_X = T,
        flag_plot = F, flag_save_plot = T, 
        save_folder = '~/Documents/0_ongoing/fe5222_project2/plots/',
        # v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v 
        xlab = 'Spot Diff ($)', ylab = 'rho', 
        fig_title = paste0('Spot Diff', ' & ', 'rho')
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    )
    print(p)
    print(range(temp_grid$Z))
    tictoc::toc()
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
}



# 4 SpotDiff & StrikePercent ----------------------------------------------
if (1) {
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    tictoc::tic()
    set.seed(2020)
    # v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v 
    temp_grid <- cbind(expand.grid(seq_SpotDiff, seq_StrikePercent), NA)
    # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    colnames(temp_grid) <- c('X', 'Y', 'Z')
    for (iter_i in 1:nrow(temp_grid)) {
        temp_grid$Z[iter_i] <- 
            # v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v 
            Spread_Call_Kirk(
                S_0_1 = 100 + temp_grid$X[iter_i], 
                Strike = temp_grid$X[iter_i] * temp_grid$Y[iter_i]
            ) - 
            Spread_Call(
                S_0_1 = 100 + temp_grid$X[iter_i], 
                Strike = temp_grid$X[iter_i] * temp_grid$Y[iter_i],
                flag_AV = T, n = 10^6) 
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    }
    p <- ContourPlots(
        temp_grid, 
        fig_subtitle = 'Differences between Kirk & MC for Spread Call',
        scale_range = c(-.1, .1), flag_X = T,
        flag_plot = F, flag_save_plot = T, 
        save_folder = '~/Documents/0_ongoing/fe5222_project2/plots/',
        # v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v 
        xlab = 'Spot Diff ($)', ylab = 'Strike in Percentage of SpotDiff', 
        fig_title = paste0('Spot Diff', ' & ', 'Strike Percent')
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    )
    print(p)
    print(range(temp_grid$Z))
    tictoc::toc()
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
}



# 5 sigma_1 & sigma_2 -----------------------------------------------------
if (1) {
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    tictoc::tic()
    set.seed(2020)
    # v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v 
    temp_grid <- cbind(expand.grid(seq_sigma_1, seq_sigma_2), NA)
    # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    colnames(temp_grid) <- c('X', 'Y', 'Z')
    for (iter_i in 1:nrow(temp_grid)) {
        temp_grid$Z[iter_i] <- 
            # v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v 
            Spread_Call_Kirk(
                sigma_1 = temp_grid$X[iter_i], sigma_2 = temp_grid$Y[iter_i] 
            ) - 
            Spread_Call(
                sigma_1 = temp_grid$X[iter_i], sigma_2 = temp_grid$Y[iter_i],
                flag_AV = T, n = 10^6) 
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    }
    p <- ContourPlots(
        temp_grid, 
        fig_subtitle = 'Differences between Kirk & MC for Spread Call',
        scale_range = c(-.3, .3), flag_X = T,
        flag_plot = F, flag_save_plot = T, 
        save_folder = '~/Documents/0_ongoing/fe5222_project2/plots/',
        # v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v 
        xlab = 'sigma_1', ylab = 'sigma_2', 
        fig_title = paste0('sigma_1', ' & ', 'sigma_2')
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    )
    print(p)
    print(range(temp_grid$Z))
    tictoc::toc()
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
}



# 6 sigma_1 & rho ---------------------------------------------------------
if (1) {
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    tictoc::tic()
    set.seed(2020)
    # v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v 
    temp_grid <- cbind(expand.grid(seq_sigma_1, seq_rho), NA)
    # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    colnames(temp_grid) <- c('X', 'Y', 'Z')
    for (iter_i in 1:nrow(temp_grid)) {
        temp_grid$Z[iter_i] <- 
            # v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v 
            Spread_Call_Kirk(
                sigma_1 = temp_grid$X[iter_i], rho = temp_grid$Y[iter_i] 
            ) - 
            Spread_Call(
                sigma_1 = temp_grid$X[iter_i], rho = temp_grid$Y[iter_i],
                flag_AV = T, n = 10^6) 
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    }
    p <- ContourPlots(
        temp_grid, 
        fig_subtitle = 'Differences between Kirk & MC for Spread Call',
        scale_range = c(-.3, .3), flag_X = T,
        flag_plot = F, flag_save_plot = T, 
        save_folder = '~/Documents/0_ongoing/fe5222_project2/plots/',
        # v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v 
        xlab = 'sigma_1', ylab = 'rho', 
        fig_title = paste0('sigma_1', ' & ', 'rho')
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    )
    print(p)
    print(range(temp_grid$Z))
    tictoc::toc()
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
}



# 7 sigma_1 & StrikePercent -----------------------------------------------
if (1) {
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    tictoc::tic()
    set.seed(2020)
    # v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v 
    temp_grid <- cbind(expand.grid(seq_sigma_1, seq_StrikePercent), NA)
    # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    colnames(temp_grid) <- c('X', 'Y', 'Z')
    for (iter_i in 1:nrow(temp_grid)) {
        temp_grid$Z[iter_i] <- 
            # v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v 
            Spread_Call_Kirk(
                sigma_1 = temp_grid$X[iter_i], Strike = 10 * temp_grid$Y[iter_i] 
            ) - 
            Spread_Call(
                sigma_1 = temp_grid$X[iter_i], Strike = 10 * temp_grid$Y[iter_i],
                flag_AV = T, n = 10^6) 
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    }
    p <- ContourPlots(
        temp_grid, 
        fig_subtitle = 'Differences between Kirk & MC for Spread Call',
        scale_range = c(-.3, .3), flag_X = T,
        flag_plot = F, flag_save_plot = T, 
        save_folder = '~/Documents/0_ongoing/fe5222_project2/plots/',
        # v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v 
        xlab = 'sigma_1', ylab = 'Strike in Percentage of SpotDiff ($10)', 
        fig_title = paste0('sigma_1', ' & ', 'Strike Percent')
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    )
    print(p)
    print(range(temp_grid$Z))
    tictoc::toc()
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
}


# 8 sigma_2 & rho ---------------------------------------------------------
if (1) {
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    tictoc::tic()
    set.seed(2020)
    # v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v 
    temp_grid <- cbind(expand.grid(seq_sigma_2, seq_rho), NA)
    # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    colnames(temp_grid) <- c('X', 'Y', 'Z')
    for (iter_i in 1:nrow(temp_grid)) {
        temp_grid$Z[iter_i] <- 
            # v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v 
            Spread_Call_Kirk(
                sigma_2 = temp_grid$X[iter_i], rho = temp_grid$Y[iter_i] 
            ) - 
            Spread_Call(
                sigma_2 = temp_grid$X[iter_i], rho = temp_grid$Y[iter_i],
                flag_AV = T, n = 10^6) 
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    }
    p <- ContourPlots(
        temp_grid, 
        fig_subtitle = 'Differences between Kirk & MC for Spread Call',
        scale_range = c(-.3, .3), flag_X = T,
        flag_plot = F, flag_save_plot = T, 
        save_folder = '~/Documents/0_ongoing/fe5222_project2/plots/',
        # v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v 
        xlab = 'sigma_2', ylab = 'rho', 
        fig_title = paste0('sigma_2', ' & ', 'rho')
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    )
    print(p)
    print(range(temp_grid$Z))
    tictoc::toc()
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
}



# 9 sigma_2 & StrikePercent -----------------------------------------------
if (1) {
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    tictoc::tic()
    set.seed(2020)
    # v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v 
    temp_grid <- cbind(expand.grid(seq_sigma_2, seq_StrikePercent), NA)
    # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    colnames(temp_grid) <- c('X', 'Y', 'Z')
    for (iter_i in 1:nrow(temp_grid)) {
        temp_grid$Z[iter_i] <- 
            # v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v 
            Spread_Call_Kirk(
                sigma_2 = temp_grid$X[iter_i], Strike = 10 * temp_grid$Y[iter_i] 
            ) - 
            Spread_Call(
                sigma_2 = temp_grid$X[iter_i], Strike = 10 * temp_grid$Y[iter_i],
                flag_AV = T, n = 10^6) 
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    }
    p <- ContourPlots(
        temp_grid, 
        fig_subtitle = 'Differences between Kirk & MC for Spread Call',
        scale_range = c(-.7, .7), flag_X = T,
        flag_plot = F, flag_save_plot = T, 
        save_folder = '~/Documents/0_ongoing/fe5222_project2/plots/',
        # v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v 
        xlab = 'sigma_2', ylab = 'Strike in Percentage of SpotDiff ($10)', 
        fig_title = paste0('sigma_2', ' & ', 'Strike Percent')
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    )
    print(p)
    print(range(temp_grid$Z))
    tictoc::toc()
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
}




# 10 rho & StrikePercent --------------------------------------------------
if (1) {
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    tictoc::tic()
    set.seed(2020)
    # v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v 
    temp_grid <- cbind(expand.grid(seq_rho, seq_StrikePercent), NA)
    # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    colnames(temp_grid) <- c('X', 'Y', 'Z')
    for (iter_i in 1:nrow(temp_grid)) {
        temp_grid$Z[iter_i] <- 
            # v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v 
            Spread_Call_Kirk(
                rho = temp_grid$X[iter_i], Strike = 10 * temp_grid$Y[iter_i] 
            ) - 
            Spread_Call(
                rho = temp_grid$X[iter_i], Strike = 10 * temp_grid$Y[iter_i],
                flag_AV = T, n = 10^6) 
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    }
    p <- ContourPlots(
        temp_grid, 
        fig_subtitle = 'Differences between Kirk & MC for Spread Call',
        scale_range = c(-.1, .1), flag_X = T,
        flag_plot = F, flag_save_plot = T, 
        save_folder = '~/Documents/0_ongoing/fe5222_project2/plots/',
        # v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v 
        xlab = 'rho', ylab = 'Strike in Percentage of SpotDiff ($10)', 
        fig_title = paste0('rho', ' & ', 'Strike Percent')
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    )
    print(p)
    print(range(temp_grid$Z))
    tictoc::toc()
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
}
