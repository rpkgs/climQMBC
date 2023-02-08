
#' @description
#' - Quantile Delta Mapping (QDM): Cannon et al. (2015).
#'
#' @param rel_change_th (Optional) A float indicating the maximum scaling factor
#' (Eq. 4 of Cannon et al. (2015)) when the denominator is below inv_mod_th.
#'
#' @param inv_mod_th (Optional) A float indicating the upper threshold of the
#' denominator of the scaling factor (Eq. 4 of Cannon et al. (2015)) to truncate
#' the scaling factor. This parameter is defined as default as the pp_threshold
#' parameter, described above.
#'
#' @rdname map_QM
#'
#' @export
map_QDM <- function(obs, mod, var, frq = "M", pp_threshold = 1, pp_factor = 1e-2,
  rel_change_th = 2, inv_mod_th = pp_threshold) {

  l <- formatQM(obs, mod, var, frq, pp_threshold, pp_factor)

  ny_obs <- l$obs$nyear
  ny_mod <- l$mod$nyear
  ny_exceed <- ny_mod - ny_obs

  obs_series <- l$obs$data
  mod_series <- l$mod$data
  nfrq <- dim(mod_series)[1]

  coef_obs <- l$obs$coef
  coef_mod <- l$mod$coef

  # 2) Assign a probability distribution function to each month for the
  #    observed and modeled data in the historical period.
  PDF_obs <- getDist(obs_series, var, coef_obs)
  PDF_mod <- getDist(mod_series[, 1:ny_obs, drop=FALSE], var, coef_mod)

  # 3) For each projected period:
  PDF_win <- matrix(0, nfrq, ny_exceed)
  Taot <- matrix(0, nfrq, ny_exceed)

  # 这一步是逐年、逐月进行的
  for (j in 1:ny_exceed) {
    win_series <- mod_series[, (1 + j):(ny_obs + j), drop=FALSE]

    mux <- apply(win_series, 1, mean, na.rm = TRUE) # Mean
    sigmax <- apply(win_series, 1, sd, na.rm = TRUE) # Standard deviation
    skewx <- apply(win_series, 1, skewness2) # Skewness
    Ln_win <- log(win_series)
    Ln_win[Im(Ln_win) != 0] <- 0
    Ln_win[!is.finite(Ln_win)] <- log(0.01)
    skewy <- apply(Ln_win, 1, skewness2) # Log-Skewness

    coef <- list(mu = mux, sigma = sigmax, skew = skewx, skewy = skewy)
    
    # a) Assign a probability distribution function to each month.
    PDF_win[, j] <- getDist(win_series, var, coef)

    # b) Apply the CDF of the projected period, evaluated with the statistics of
    # this period, to the last data of the period (getCDF). Eq. 3 of Cannon et al. (2015).
    
    Taot[, j] <- getCDF(PDF_win[, j], matrix(mod_series[, ny_obs + j]), coef)
  }

  # 4) Apply the inverse CDF:
  # a) Of the observed data, evaluated with the statistics of the observed data
  #    in the historical period, to the probabilities obtained from 3b). Eq. 5 of Cannon et al. (2015).
  inv_obs <- getCDFinv(PDF_obs, Taot, coef_obs)

  # b) Of the modeled data, evaluated with the statistics of the observed data
  #    in the historical period, to the probabilities obtained from 3b). Eq.s 4 of Cannon et al. (2015).
  inv_mod <- getCDFinv(PDF_mod, Taot, coef_mod)

  # 5) Get the delta factor or relative change and apply it to the value
  #    obtained in 4b). Eq. 4 and 6 of Cannon et al. (2015).
  if (var == 1) {
    DM <- mod_series[, (ny_obs + 1):ny_mod] / inv_mod
    DM[(DM > rel_change_th) & (inv_mod < inv_mod_th)] <- rel_change_th
    QDM <- inv_obs * DM
  } else {
    DM <- mod_series[, (ny_obs + 1):ny_mod] - inv_mod
    QDM <- inv_obs + DM
  }
  QDM <- matrix(QDM)

  # 6) Perform QM for the historical period.
  mod_h <- mod_series[, 1:ny_obs]
  mod_h <- matrix(mod_h)
  QM_series <- QM(obs, mod_h, var, frq)
  QDM_series <- c(QM_series, QDM)
  if (var == 1) {
    QDM_series[QDM_series < pp_threshold] <- 0
  }
  return(QDM_series)
}
