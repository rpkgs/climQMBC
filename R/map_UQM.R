#' @description
#' - Unbiased Quantile Mapping (UQM): Chadwick et al. (2021). 
#' 
#' @details 
#' Correction is performed to monthly or annual precipitation or
#' temperature data in a single location. An independent probability
#' distribution function is assigned to each month and to each projected period
#' based on the Kolmogorov-Smirnov test.
#'
#' _For DQM, QDM, and UQM:_
#'
#'    If annual frequency is specified, this is applied to the complete period.
#' Each projected period considers a window whose length is equal to the number
#' of years of the historical period and ends in the analyzed year.
#'
#'    Correction of the historical period is performed by the Quantile Mapping (QM)
#' method, as described by Cannon et al. (2015).
#'
#' @references 2.Chadwick et al. (2022) `[under revision]`
#'
#' @rdname map_QM
#' @export
map_UQM <- function(obs, mod, var, frq = "M", pp_threshold = 1, pp_factor = 1e-2) {

  l <- formatQM(obs, mod, var, frq, pp_threshold, pp_factor)

  ny_obs <- l$obs$nyear
  ny_mod <- l$mod$nyear
  obs_series <- l$obs$data
  mod_series <- l$mod$data
  nfrq <- dim(mod_series)[1]
  
  coef_obs <- l$obs$coef
  coef_mod <- l$mod$coef

  mu_obs     <- l$obs$coef$mu
  std_obs    <- l$obs$coef$sigma
  skew_obs   <- l$obs$coef$skew
  skewny_obs <- l$obs$coef$skewy

  mu_mod     <- l$mod$coef$mu
  std_mod    <- l$mod$coef$sigma
  skew_mod   <- l$mod$coef$skew
  skewny_mod <- l$mod$coef$skewy

  
  PDF_obs <- getDist(obs_series, var, coef_obs)

  # 3) For each projected period, get the delta factor (delta) and time
  #    dependent (aster) statistics (mean, standard deviation, skewness, and
  #    log skewness). Eq.s X to X in Chadwick et al. (2021).
  xbarmt <- matrix(0, dim(mod_series)[1], ny_mod - ny_obs)
  xhatmt <- matrix(0, dim(xbarmt)[1], dim(xbarmt)[2])
  skwmt <- matrix(0, dim(xbarmt)[1], dim(xbarmt)[2])
  Lskwmt <- matrix(0, dim(xbarmt)[1], dim(xbarmt)[2])

  Dmu <- matrix(0, dim(xbarmt)[1], dim(xbarmt)[2])
  Dsigma <- matrix(0, dim(xbarmt)[1], dim(xbarmt)[2])
  Dskw <- matrix(0, dim(xbarmt)[1], dim(xbarmt)[2])
  DLskw <- matrix(0, dim(xbarmt)[1], dim(xbarmt)[2])
  muAster <- matrix(0, dim(xbarmt)[1], dim(xbarmt)[2])
  sigmaAster <- matrix(0, dim(xbarmt)[1], dim(xbarmt)[2])
  skwAster <- matrix(0, dim(xbarmt)[1], dim(xbarmt)[2])
  LskwAster <- matrix(0, dim(xbarmt)[1], dim(xbarmt)[2])

  if (var == 1) { # Precipitation
    for (i in 1:dim(xbarmt)[1]) {
      for (j in 1:dim(xbarmt)[2]) {
        xbarmt[i, j] <- mean(mod_series[i, (j + 1):(ny_obs + j)], na.rm = TRUE)
        xhatmt[i, j] <- sd(mod_series[i, (j + 1):(ny_obs + j)], na.rm = TRUE)

        skwmt[i, j] <- e1071::skewness(mod_series[i, (j + 1):(ny_obs + j)], na.rm = TRUE, type = 2)
        Lnskwmt <- log(mod_series[i, (j + 1):(ny_obs + j)])
        Lnskwmt[Im(Lnskwmt) != 0] <- 0
        Lnskwmt[!is.finite(Lnskwmt)] <- log(0.01)
        Lskwmt[i, j] <- e1071::skewness(Lnskwmt, na.rm = TRUE, type = 2)

        Dmu[i, j] <- (xbarmt[i, j] - mu_mod[i]) / mu_mod[i]
        Dsigma[i, j] <- (xhatmt[i, j] - std_mod[i]) / std_mod[i]
        Dskw[i, j] <- (skwmt[i, j] - skew_mod[i]) / skew_mod[i]
        DLskw[i, j] <- (Lskwmt[i, j] - skewny_mod[i]) / skewny_mod[i]

        muAster[i, j] <- mu_obs[i] * (1 + Dmu[i, j])
        sigmaAster[i, j] <- std_obs[i] * (1 + Dsigma[i, j])
        skwAster[i, j] <- skew_obs[i] * (1 + Dskw[i, j])
        LskwAster[i, j] <- skewny_obs[i] * (1 + DLskw[i, j])
      }
    }
  } else { # Temperature
    for (i in 1:dim(xbarmt)[1]) {
      for (j in 1:dim(xbarmt)[2]) {
        xbarmt[i, j] <- mean(mod_series[i, (j + 1):(ny_obs + j)], na.rm = TRUE)
        xhatmt[i, j] <- sd(mod_series[i, (j + 1):(ny_obs + j)], na.rm = TRUE)

        skwmt[i, j] <- e1071::skewness(mod_series[i, (j + 1):(ny_obs + j)], na.rm = TRUE, type = 2)
        Lnskwmt <- log(mod_series[i, (j + 1):(ny_obs + j)])
        Lnskwmt[Im(Lnskwmt) != 0] <- 0
        Lnskwmt[!is.finite(Lnskwmt)] <- log(0.01)
        Lskwmt[i, j] <- e1071::skewness(Lnskwmt, na.rm = TRUE, type = 2)

        Dmu[i, j] <- xbarmt[i, j] - mu_mod[i]
        Dsigma[i, j] <- xhatmt[i, j] - std_mod[i]
        Dskw[i, j] <- skwmt[i, j] - skew_mod[i]
        DLskw[i, j] <- Lskwmt[i, j] - skewny_mod[i]

        muAster[i, j] <- mu_obs[i] + (Dmu[i, j])
        sigmaAster[i, j] <- std_obs[i] + (Dsigma[i, j])
        skwAster[i, j] <- skew_obs[i] + (Dskw[i, j])
        LskwAster[i, j] <- skewny_obs[i] + (DLskw[i, j])
      }
    }
  }

  # 4) For each projected period:
  PDF_win <- matrix(0, dim(mod_series)[1], ny_mod - ny_obs)
  Taot <- matrix(0, dim(mod_series)[1], ny_mod - ny_obs)
  UQM <- matrix(0, dim(mod_series)[1], ny_mod - ny_obs)
  for (j in 1:dim(Taot)[2]) {
    win_series <- matrix(mod_series[, (1 + j):(ny_obs + j)], nrow = dim(mod_series)[1])

    mux <- apply(win_series, 1, mean, na.rm = TRUE) # Mean
    sigmax <- apply(win_series, 1, sd, na.rm = TRUE) # Standard deviation
    skewx <- apply(win_series, 1, e1071::skewness, na.rm = TRUE, type = 2) # Skewness
    Ln_win <- log(win_series)
    Ln_win[Im(Ln_win) != 0] <- 0
    Ln_win[!is.finite(Ln_win)] <- log(0.01)
    skewy <- apply(Ln_win, 1, e1071::skewness, na.rm = TRUE, type = 2) # Log-Skewness

    coef <- list(mu = mux, sigma = sigmax, skew = skewx, skewy = skewy)

    # a) Assign a probability distribution function to each month. If
    #    annual frequency is specified, this is applied to the complete
    #    period (getDist).
    PDF_win[, j] <- getDist(win_series, var, coef)

    # b) Apply the CDF of the projected
    #    period, evaluated with the statistics of this period, to the last
    #    data of the period (getCDF).
    #    Eq. X of Chadwick et al. (2021).
    Taot[, j] <- getCDF(PDF_win[, j], matrix(mod_series[, ny_obs + j]), coef)
  }

  # 5) Apply the inverse CDF of the observed
  #    data, evaluated with the time dependent statistics, to the values
  #    obtained in 4b) (getCDFinv). Eq.
  #    X of Chadwick et al. (2021).
  for (yr in 1:dim(Taot)[2]) {
    coef <- list(
      mu    = matrix(muAster[, yr]),
      sigma = matrix(sigmaAster[, yr]), 
      skew  = matrix(skwAster[, yr]),
      skewy = matrix(LskwAster[, yr])
    )
    UQM[, yr] <- getCDFinv(PDF_obs, matrix(Taot[, yr]), coef)
  }

  UQM <- matrix(UQM)

  # 6) Perform QM for the historical period
  mod_h <- mod_series[, 1:ny_obs]
  mod_h <- matrix(mod_h)
  QM_series <- QM(obs, mod_h, var, frq)
  UQM_series <- c(QM_series, UQM)
  if (var == 1) {
    UQM_series[UQM_series < pp_threshold] <- 0
  }
  return(UQM_series)
}
