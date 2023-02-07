#' @description
#' - Scaled Distribution Mapping (SDM): 
#'
#'   Proposed by Switanek et al. (2017). The historical period of the modeled
#' series is bias corrected as a scenario where the complete series is replaced
#' by the bias corrected series. On the other hand, for each year after the
#' historical period the method considers a projected period as a window whose
#' length is equal to the number of years of the historical period and ends in
#' the analyzed year. From the bias corrected series, only the last year is
#' saved and all the others are discarded. This approach aims to build a
#' continuum for the bias corrected series.
#'
#' @references
#'
#' 3.Switanek, B. M., Troch, P., Castro, L. C., Leuprecht, A., Chang, H. I.,
#' Mukherjee, R., and Demaria, M. C. E. (2017) Scaled distribution mapping: A
#' bias correction method that preserves raw climate model projected changes.
#' Hydrology &amp; Earth System Sciences, 21, 2649-2666,
#' https://doi.org/10.5194/hess-21-2649-2017.
#' 
#' @rdname map_QM
#' 
#' @export
map_SDM <- function(obs, mod, var, frq = "M", pp_threshold = 1, pp_factor = 1e-2) {

  lower_lim <- pp_threshold
  CDF_th <- 10^-3

  l <- formatQM(obs, mod, var, frq, pp_threshold, pp_factor)

  ny_obs     <- l$obs$nyear
  ny_mod     <- l$mod$nyear

  obs_series <- l$obs$data
  mod_series <- l$mod$data
  nfrq <- dim(mod_series)[1]

  SDM <- matrix(0, nfrq, ny_mod - ny_obs)
  SDM_h <- matrix(0, dim(obs_series)[1], ny_obs)
  for (m in 1:nfrq) {
    ## 2) Historical period:
    # a) [Switanek et al. (2017), step 1)]
    #     For temperature, get the detrended modeled and observed series in
    #     the historical period.
    #     For precipitation, get the rainday values and its frequency for
    #     the modeled and observed series in the historical period.
    if (var == 0) {
      D_obs <- pracma::detrend(obs_series[m, ])
      D_mod <- pracma::detrend(mod_series[m, 1:ny_obs])

      mu_obs <- mean(obs_series[m, ])
      mu_mod <- mean(mod_series[m, 1:ny_obs])
    } else {
      D_obs <- sort(obs_series[m, obs_series[m, ] > lower_lim])
      D_mod <- sort(mod_series[m, 1:ny_obs][mod_series[m, 1:ny_obs] > lower_lim])

      freq_obs <- length(D_obs) / dim(obs_series)[2]
      freq_mod <- length(D_mod) / length(mod_series[m, 1:ny_obs])
    }
    ## b) [Switanek et al. (2017), step 2)]
    #    For temperature, fit normal probability distribution to the
    #    detrended observed and modeled series in the historical period.
    #    For precipitation, fit gamma distribution parameters by maximum
    #    likelihood to the rainday values of the observed and modeled
    #    series in the historical period.
    #    Using these fitted distributions, get the corresponding cumulative
    #    distribution values. Tail probabilities are set to (0 + threshold)
    #    for temperature and to (1 - threshold) for temperature and
    #    precipitation. Threshold is set in the first lines of this
    #    function. Default is CDF_th = 10^-3.
    if (var == 0) {
      mu_obsD <- mean(D_obs)
      mu_modD <- mean(D_mod)
      sigma_obsD <- sd(D_obs)
      sigma_modD <- sd(D_mod)

      CDF_obs <- pnorm(sort(D_obs), mu_obsD, sigma_obsD)
      CDF_mod <- pnorm(sort(D_mod), mu_modD, sigma_modD)
      CDF_obs[CDF_obs < CDF_th] <- CDF_th
      CDF_mod[CDF_mod < CDF_th] <- CDF_th
      CDF_obs[CDF_obs > 1 - CDF_th] <- 1 - CDF_th
      CDF_mod[CDF_mod > 1 - CDF_th] <- 1 - CDF_th
    } else {
      fit_obs <- fitdistrplus::fitdist(D_obs, distr = "gamma", method = "mle")
      fit_mod <- fitdistrplus::fitdist(D_mod, distr = "gamma", method = "mle")
      CDF_obs <- pgamma(D_obs, shape = fit_obs[[1]][1], rate = fit_obs[[1]][2])
      CDF_mod <- pgamma(D_mod, shape = fit_mod[[1]][1], rate = fit_mod[[1]][2])
      CDF_obs[CDF_obs > 1 - CDF_th] <- 1 - CDF_th
      CDF_mod[CDF_mod > 1 - CDF_th] <- 1 - CDF_th
    }

    # 3) Projected periods:
    for (j in 0:(ny_mod - ny_obs)) {
      # c) Initialize correction array.
      corr_temp <- matrix(0, ny_obs, 1)

      # d) Define projected window.
      win_series <- mod_series[m, (1 + j):(ny_obs + j)]

      # e) [Switanek et al. (2017), step 1)]
      #    For temperature, get the detrended series of the projected
      #    period.
      #    For precipitation, get the rainday values, its frequency, and
      #    expected raindays for the projected period.
      #    Get the index of the sorted detrended or rainday values.
      if (var == 0) {
        D_win <- pracma::detrend(win_series)
        exp_D <- length(win_series)
        win_argsort <- sort(D_win, index.return = TRUE)[[2]]
        mu_win <- mean(win_series)
      } else {
        D_win <- sort(win_series[win_series > lower_lim])
        exp_D <- min(round(length(D_win) * freq_obs / freq_mod), length(win_series))

        win_argsort <- sort(win_series, index.return = TRUE)[[2]]
      }
      # f) [Switanek et al. (2017), step 2)]
      #    For temperature, fit normal probability distribution to the
      #    detrended values of the projected period.
      #    For precipitation, fit gamma distribution parameters by
      #    maximum likelihood to the rainday values of the projected
      #    period.
      #    Using these fitted distributions, get the corresponding
      #    cumulative distribution values. Tail probabilities are set to
      #    (0 + threshold) for temperature and to (1 - threshold) for
      #    temperature and precipitation. Threshold is set in the first
      #    lines of this function. Default is CDF_th = 10^-3.
      if (var == 0) {
        mu_winD <- mean(D_win, na.rm = TRUE)
        sigma_winD <- sd(D_win, na.rm = TRUE)

        diff_win <- win_series - D_win

        CDF_win <- pnorm(sort(D_win), mu_winD, sigma_winD)
        CDF_win[CDF_win < CDF_th] <- CDF_th
        CDF_win[CDF_win > 1 - CDF_th] <- 1 - CDF_th
      } else {
        fit_win <- fitdistrplus::fitdist(D_win, distr = "gamma", method = "mle")
        CDF_win <- pgamma(D_win, shape = fit_win[[1]][1], rate = fit_win[[1]][2])
        CDF_win[CDF_win > 1 - CDF_th] <- 1 - CDF_th
      }

      # g) [Switanek et al. (2017), step 3)]
      #    Get the scaling between the model projected period and
      #    historical period distributions.
      if (var == 0) {
        SF <- (sigma_obsD / sigma_modD) * (qnorm(CDF_win, mu_winD, sigma_winD) - qnorm(CDF_win, mu_modD, sigma_modD))
      } else {
        SF <- qgamma(CDF_win, shape = fit_win[[1]][1], rate = fit_win[[1]][2]) / qgamma(CDF_win, shape = fit_mod[[1]][1], rate = fit_mod[[1]][2])
      }

      # h) Interpolate observed and modeled CDF of the historical period
      #    to the length of the projected period.
      obs_cdf_intpol <- pracma::interp1(pracma::linspace(1, length(D_obs), length(D_obs)), CDF_obs, pracma::linspace(1, length(D_obs), length(D_win)))
      mod_cdf_intpol <- pracma::interp1(pracma::linspace(1, length(D_mod), length(D_mod)), CDF_mod, pracma::linspace(1, length(D_mod), length(D_win)))

      # i) [Switanek et al. (2017), steps 4 and 5)]
      #    Get recurrence intervals and its scaled version for the
      #    observed, historical modeled and projected period modeled
      #    CDFs.
      if (var == 0) {
        RI_obs <- 1 / (0.5 - abs(obs_cdf_intpol - 0.5))
        RI_mod <- 1 / (0.5 - abs(mod_cdf_intpol - 0.5))
        RI_win <- 1 / (0.5 - abs(CDF_win - 0.5))
        RI_scaled <- RI_obs * RI_win / RI_mod

        CDF_scaled <- sign(obs_cdf_intpol - 0.5) * (1 - 1 / RI_scaled)
        CDF_scaled[CDF_scaled < 0] <- CDF_scaled[CDF_scaled < 0] + 1

        CDF_scaled[CDF_scaled < CDF_th] <- CDF_th
        CDF_scaled[CDF_scaled > 1 - CDF_th] <- 1 - CDF_th
      } else {
        RI_obs <- 1 / (1 - obs_cdf_intpol)
        RI_mod <- 1 / (1 - mod_cdf_intpol)
        RI_win <- 1 / (1 - CDF_win)
        RI_scaled <- RI_obs * RI_win / RI_mod
        RI_scaled[RI_scaled < 1] <- 1

        CDF_scaled <- sort(1 - 1 / (RI_scaled))
      }

      # j) [Switanek et al. (2017), step 6)]
      #    Get the initial bias corrected values. For precipitation,
      #    these values are interpolated to the length of the expected
      #    raindays.
      if (var == 0) {
        xvals <- qnorm(sort(CDF_scaled), mu_obsD, sigma_obsD) + SF
        xvals <- xvals - mean(xvals) + mu_obs + (mu_win - mu_mod)
      } else {
        xvals <- qgamma(CDF_scaled, shape = fit_obs[[1]][1], rate = fit_obs[[1]][2]) * SF
        if (length(D_win) > exp_D) {
          xvals <- pracma::interp1(pracma::linspace(1, length(D_win), length(D_win)), xvals, pracma::linspace(1, length(D_win), exp_D))
        } else {
          xvals <- c(matrix(0, 1, exp_D - length(D_win)), xvals)
        }
      }
      # k) [Switanek et al. (2017), step 7)]
      #    Bias corrected values are placed back matching the higher bias
      #    corrected values with the higher rainday or detrended values.
      #    For temperature, the trend of the projected period is added
      #    back.
      corr_temp[win_argsort[(length(win_argsort) - exp_D + 1):length(win_argsort)]] <- matrix(xvals)
      if (var == 0) {
        corr_temp <- corr_temp + diff_win - mu_win
      }
      # l) If the projected period is the historical period (period 0,
      #    j=0) save the complete bias corrected series.
      #    If the projected period is not the historical period (j>0),
      #    save the value of the last year.
      if (j == 0) {
        SDM_h[m, ] <- corr_temp
      } else {
        SDM[m, j] <- corr_temp[length(corr_temp)]
      }
    }
  }

  SDM <- matrix(SDM)
  SDM_h <- matrix(SDM_h)
  SDM_series <- c(SDM_h, SDM)
  if (var == 1) {
    SDM_series[SDM_series < pp_threshold] <- 0
  }
  return(SDM_series)
}
