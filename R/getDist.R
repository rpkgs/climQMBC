#' Get probability distribution function
#'
#' @section Distributions: The available distributions are:
#' 1) Normal distribution;
#' 2) Log-Normal distribution;
#' 3) Gamma 2 parameters distribution;
#' 4) Gamma 3 parameters distribution (Pearson 3 parameters distribution);
#' 5) Log-Gamma 3 parameters distribution (Log-Pearson 3 parameters distribution);
#' 6) Gumbel distribution;
#' 7) Exponential distribution
#' 
#' - For precipitation, only 2), 3) and 5) are considered (1, 4, 6, and 7 are
#' discarded).
#'
#' - For series with negative values, only 1), 3), 4), 6), and 7) are considered
#' (2, 3 and 5 are discarded).
#'
#' @description
#'
#' Get probability distribution function for each month of the period.
#'
#' - For monthly data, it will have 12 rows and each row will represent a month.
#'
#' - For annual data the series will have only one row.
#'
#' This function assigns an independent probability distribution function to
#' each row of the input series by comparing the empirical probability
#' distribution function with seven distributions based on the
#' Kolmogorov-Smirnov (KS) test.
#'
#' Only strictly positive distributions are considered for precipitation and
#' strictly positive distributions are discarded if the series has negative
#' values.
#'
#' @param series A matrix of monthly or annual data (temperature or
#' precipitation). If the series consider monthly data, it will have 12 rows and
#' each row will represent a month. For annual data the series will have only
#' one row.
#'
#' @param coef A named list, with the elements of:
#' - `mu`: A column vector of mean values of the series. `[12,1]` if the series
#' consider monthly data and `[1,1]` if the series consider annual data.
#'
#' - `sigma`: A column vector of standard deviation of the series. `[12,1]` if
#' the series consider monthly data and `[1,1]` if the series consider annual
#' data.
#'
#' - `skew`: A column vector of skewness of the series. `[12,1]` if the series
#' consider monthly data and `[1,1]` if the series consider annual data.
#'
#' - `skewy`: A column vector of skewness of the logarithm of the series.
#' `[12,1]` if the series consider monthly data and `[1,1]` if the series
#' consider annual data.
#'
#' @param var A flag that identifies if data are temperature or precipitation.
#' This flag tells the getDist function if it has to discard distribution
#' functions that allow negative numbers.
#' - Temperature: `var = 0`;
#' - Precipitation: `var = 1`, which not allow negative numbers.
#'
#' @return `PDF`: A column vector with an ID for the resulting distribution from
#' the KS test. `[12,1]` if the series consider monthly data and `[1,1]` if the
#' series consider annual data.
#'
#' The `ID` is related to the numeration of the distribution listed in the
#' description of this function. This `ID` is used in the [getCDF] and
#' [getCDFinv]
#'
#' @example R/example/ex-getDist.R
#'
#' @export
getDist <- function(series, var, coef) {
  # mu, sigma, skew, skewy
  mu    <- coef$mu
  sigma <- coef$sigma
  skew  <- coef$skew
  skewy <- coef$skewy
  
  # 1) Get the number of years to compute the empirical distribution in step
  #    3) and get the number of rows of the input series.
  ny <- dim(series)[2]
  n <- dim(series)[1]

  # 2) Initialize column vectors for the statistics needed for the available
  #    probability distribution functions.
  PDF    <- matrix(0, n, 1)
  sigmay <- matrix(0, n, 1)
  muy    <- matrix(0, n, 1)
  A      <- matrix(0, n, 1)
  B      <- matrix(0, n, 1)
  Alp    <- matrix(0, n, 1)
  Bet    <- matrix(0, n, 1)
  Gam    <- matrix(0, n, 1)
  Alpy   <- matrix(0, n, 1)
  Bety   <- matrix(0, n, 1)
  Gamy   <- matrix(0, n, 1)
  a      <- matrix(0, n, 1)
  u      <- matrix(0, n, 1)

  # 3) Perform the Kolmogorov-Smirnov test for each row.
  for (m in 1:n) {
    # a) Get empirical distribution.
    sortdata <- sort(series[m, ])
    probEmp <- seq(from = 1 / (ny + 1), to = ny / (ny + 1), by = 1 / (ny + 1))

    # b) Compare distributions.
    # i) Normal distribution.
    normal <- pnorm(sortdata, mu[m], sigma[m])
    KSnormal <- max(abs(normal - probEmp))

    # ii) Log Normal distribution.
    if (any(series[m, ] < 0)) {
      KSlognormal <- 1
    } else {
      sigmay[m] <- sqrt(log(1 + (sigma[m] / mu[m])^2))
      muy[m] <- log(mu[m]) - (sigmay[m]^2) / 2
      lognormal <- plnorm(sortdata, muy[m], sigmay[m])
      KSlognormal <- max(abs(lognormal - probEmp))
    }

    # iii) Gamma 2 parameters distribution.
    if (any(series[m, ] < 0)) {
      KSgammaII <- 1
    } else {
      A[m] <- (sigma[m]^2) / mu[m]
      B[m] <- (mu[m] / sigma[m])^2
      GammaII <- pgamma(sortdata, shape = B[m], scale = A[m])
      KSgammaII <- max(abs(GammaII - probEmp))
    }

    # % iv) Gamma 3 parameters distribution.
    #     (Pearson 3 parameters distribution)
    Bet[m] <- (2 / skew[m])^2
    Alp[m] <- sigma[m] / sqrt(Bet[m])
    Gam[m] <- mu[m] - (Alp[m] * Bet[m])
    GammaIII <- pgamma((sortdata - Gam[m]), shape = Bet[m], scale = Alp[m])
    KSgammaIII <- max(abs(GammaIII - probEmp))

    # v) Log-Gamma 3 parameters distribution.
    #    (Log-Pearson 3 parameters distribution)
    if (any(series[m, ] < 0)) {
      KSLpIII <- 1
    } else {
      Bety[m] <- (2 / skewy[m])^2
      Alpy[m] <- sigmay[m] / sqrt(Bety[m])
      Gamy[m] <- muy[m] - (Alpy[m] * Bety[m])
      Lnsortdata <- log(sortdata)
      Lnsortdata[Im(Lnsortdata) != 0] <- 0
      Lnsortdata[!is.finite(Lnsortdata)] <- log(0.01)
      LpIII <- pgamma((Lnsortdata - Gamy[m]), shape = Bety[m], scale = Alpy[m])
      KSLpIII <- max(abs(LpIII - probEmp))
    }

    # vi) Gumbel distribution.
    Sn <- pi / sqrt(6)
    yn <- 0.5772
    a[m] <- Sn / sigma[m]
    u[m] <- mu[m] - (yn / a[m])
    gumbel <- exp(-exp(-a[m] * (sortdata - u[m])))
    KSgumbel <- max(abs(gumbel - probEmp))

    # vii) Exponential distribution.
    gamexp <- mu[m] - sigma[m]
    exponential <- pmax(1 - exp(-1 / sigma[m] * (sortdata - gamexp)), 0)
    KSexponential <- max(abs(exponential - probEmp))

    # c) If variable is precipitation, set KS=1 to distributions that allow
    #    negative values (this will discard those distributions).
    if (var == 1) {
      KSnormal <- 1
      KSgammaIII <- 1
      KSgumbel <- 1
      KSexponential <- 1
    }

    # d) The distribution with lower KS value is considered for each month.
    bestPDF <- which.min(c(KSnormal, KSlognormal, KSgammaII, KSgammaIII, KSLpIII, KSgumbel, KSexponential))

    PDF[m] <- bestPDF
  }
  return(PDF)
}

#' Get probability of a set of values
#'
#' Evaluates each row of the series in the respective cumulative distribution
#' function assigned by the Kolmogorov-Smirnov (KS) test in the [getDist].
#'
#' @inheritSection getDist Distributions
#' @inheritParams getDist
#' @param PDF A column vector with an ID for the resulting distribution from the
#' KS test. `[12,1]` if the series consider monthly data and `[1,1]` if the
#' series consider annual data. The ID is related to the numeration of the
#' distribution listed in the description of this function.
#'
#' @param series A matrix of monthly or annual data (temperature or
#' precipitation). If the series consider monthly data, it will have 12 rows and
#' each row will represent a month. For annual data the series will have only
#' one row.
#'
#' @return `Taot`: A column vector with the non-exceedance probability for each
#' row of the input series. `[12,1]` if the series consider monthly data and
#' `[1,1]` if the series consider annual data.
#'
#' @example R/example/ex-getDist.R
#' @export
getCDF <- function(PDF, series, coef) {
  # mu, sigma, skew, skewy
  mu    <- coef$mu
  sigma <- coef$sigma
  skew  <- coef$skew
  skewy <- coef$skewy

  # 1) Get the number of rows and years of the series.
  n_m <- dim(series)[1]
  n_y <- dim(series)[2]

  # 2) Compute the cumulative distribution function to the values of each
  #    row, based on the distribution assigned in the getDist function of the
  #    climQMBC package.
  Taot <- matrix(0, n_m, n_y)
  for (m in 1:n_m) {
    if (PDF[m] == 1) { # i) Normal distribution.
      Taot[m, ] <- pnorm(series[m, ], mu[m], sigma[m])
    } else if (PDF[m] == 2) { # ii) Log-Normal distribution.
      sigmay <- sqrt(log(1 + (sigma[m] / mu[m])^2))
      muy <- log(mu[m]) - (sigmay^2) / 2
      Taot[m, ] <- plnorm(series[m, ], muy, sigmay)
    } else if (PDF[m] == 3) { # iii) Gamma 2 parameters distribution.
      A <- (sigma[m]^2) / mu[m]
      B <- (mu[m] / sigma[m])^2
      Taot[m, ] <- pgamma(series[m, ], shape = B, scale = A)
    } else if (PDF[m] == 4) { # iv) Gamma 3 parameters distribution.
      Bet <- (2 / skew[m])^2
      Alp <- sigma[m] / sqrt(Bet)
      Gam <- mu[m] - (Alp * Bet)
      Taot[m, ] <- pgamma((series[m, ] - Gam), shape = Bet, scale = Alp)
    } else if (PDF[m] == 5) { # v) Log-Gamma 3 parameters distribution.
      Bety <- (2 / skewy[m])^2
      sigmay <- sqrt(log(1 + (sigma[m] / mu[m])^2))
      Alpy <- sigmay / sqrt(Bety)
      muy <- log(mu[m]) - (sigmay^2) / 2
      Gamy <- muy - (Alpy * Bety)
      Lnsortdata <- log(series[m, ])
      Lnsortdata[Im(Lnsortdata) != 0] <- log(0.01)
      Lnsortdata[!is.finite(Lnsortdata)] <- log(0.01)
      Taot[m, ] <- pgamma((Lnsortdata - Gamy), shape = Bety, scale = Alpy)
    } else if (PDF[m] == 6) { # vi) Gumbel distribution.
      Sn <- pi / sqrt(6)
      yn <- 0.5772
      a <- Sn / sigma[m]
      u <- mu[m] - (yn / a)
      Taot[m, ] <- exp(-exp(-a * (series[m, ] - u)))
    } else if (PDF[m] == 7) { # vii) Exponential distribution.
      gamexp <- mu[m] - sigma[m]
      Taot[m, ] <- pmax(1 - exp(-1 / sigma[m] * (series[m, ] - gamexp)), 0)
    }
  }

  # 3) Tail probabilities are set to (0 +  threshold) and (1 - threshold) to
  #    avoid numerical errors.
  th <- 10^-3
  Taot[Taot > 1 - th] <- 1 - th
  Taot[Taot < th] <- th

  return(Taot)
}


#' Get the value associated to a certain probability
#'
#' Evaluates the probability, Taot, in the respective inverse cumulative
#' distribution function assigned by the Kolmogorov-Smirnov (KS) test in the
#' [getDist].
#'
#' @inheritSection getDist Distributions
#' @inheritParams getDist
#'
#' @param PDF A column vector with an ID for the resulting distribution from the
#' KS test. `[12,1]` if the series consider monthly data and `[1,1]` if the
#' series consider annual data. The ID is related to the numeration of the
#' distribution listed in the description of this function.
#'
#' @param Taot A column vector with the non-exceedance probability for each row
#' of the input series. `[12,1]` if the series consider monthly data and `[1,1]`
#' if the series consider annual data.
#'
#' @return `xhat`: A column vector with the values obtained when the inverse
#' cumulative distribution function is applied. `[12,1]` if the series consider
#' monthly data and `[1,1]` if the series consider annual data.
#'
#' @example R/example/ex-getDist.R
#' @export
getCDFinv <- function(PDF, Taot, coef) {
  # mu, sigma, skew, skewy
  mu    <- coef$mu
  sigma <- coef$sigma
  skew  <- coef$skew
  skewy <- coef$skewy

  # 1) Get the number of rows and years of the series.
  n_m <- dim(Taot)[1]
  n_y <- dim(Taot)[2]

  # 2) Compute the inverse cumulative distribution function to the values of
  #    each row, based on the distribution assigned in the getDist function
  #    of the climQMBC package.
  xhat <- matrix(0, n_m, n_y)
  for (m in 1:n_m) {
    if (PDF[m] == 1) { # i) Normal distribution.
      xhat[m, ] <- qnorm(Taot[m, ], mu[m], sigma[m])
    } else if (PDF[m] == 2) { # ii) Log-Normal distribution.
      sigmay <- sqrt(log(1 + (sigma[m] / mu[m])^2))
      muy <- log(mu[m]) - (sigmay^2) / 2
      xhat[m, ] <- qlnorm(Taot[m, ], muy, sigmay)
    } else if (PDF[m] == 3) { # iii) Gamma 2 parameters distribution.
      A <- (sigma[m]^2) / mu[m]
      B <- (mu[m] / sigma[m])^2
      xhat[m, ] <- qgamma(Taot[m, ], shape = B, scale = A)
    } else if (PDF[m] == 4) { # iv) Gamma 3 parameters distribution.
      Bet <- (2 / skew[m])^2
      Alp <- sigma[m] / sqrt(Bet)
      Gam <- mu[m] - (Alp * Bet)
      xhat[m, ] <- qgamma(Taot[m, ], shape = Bet, scale = Alp) + Gam
    } else if (PDF[m] == 5) { # v) Log-Gamma 3 parameters distribution.
      Bety <- (2 / skewy[m])^2
      sigmay <- sqrt(log(1 + (sigma[m] / mu[m])^2))
      Alpy <- sigmay / sqrt(Bety)
      muy <- log(mu[m]) - (sigmay^2) / 2
      Gamy <- muy - (Alpy * Bety)
      xhat[m, ] <- exp(qgamma(Taot[m, ], shape = Bety, scale = Alpy) + Gamy)
    } else if (PDF[m] == 6) { # vi) Gumbel distribution.
      Sn <- pi / sqrt(6)
      yn <- 0.5772
      a <- Sn / sigma[m]
      u <- mu[m] - (yn / a)
      xhat[m, ] <- u - log(-log(Taot[m, ])) / a
    } else if (PDF[m] == 7) { # vii) Exponential distribution.
      gamexp <- mu[m] - sigma[m]
      xhat[m, ] <- gamexp - (sigma[m] * log(1 - Taot[m, ]))
    }
  }

  return(xhat)
}
