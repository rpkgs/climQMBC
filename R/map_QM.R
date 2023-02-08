#' Quantile Mapping functions
#' 
#' @description
#' - Quantile Mapping (QM): Cannon et al. (2015).
#' 
#' @inheritParams formatQM
#'
#' @return A column vector of monthly or annual modeled data (temperature or
#' precipitation) corrected by the QM method. 
#' 
#' - If monthly frequency is specified, the length of this vector is 12 times
#' the number of observed years `[12 x nyear_mod, 1]`. 
#' 
#' - If annual frequency is specified, the length of this vector is the number
#' of observed  years `[nyear_mod, 1]`.
#' 
#' @references
#' 
#' 1.Cannon, A. J., S. R. Sobie, and T. Q. Murdock, 2015: Bias correction of GCM
#' precipitation by quantile mapping: How well do methods preserve changes in
#' quantiles and extremes? J. Climate, 28(17), 6938-6959,
#' https://doi.org/10.1175/JCLI-D-14-00754.1
#'
#' @export
map_QM <- function(obs, mod, var, frq = "M", pp_threshold = 1, pp_factor = 1e-2) {
  l <- formatQM(obs, mod, var, frq, pp_threshold, pp_factor)

  ny_obs     <- l$obs$nyear
  obs_series <- l$obs$data
  mod_series <- l$mod$data

  coef_obs <- l$obs$coef
  coef_mod <- l$mod$coef

  # 2) Assign a probability distribution function to each month for the
  #    observed and modeled data in the historical period. If annual
  #    frequency is specified, this is applied to the complete historical
  #    period (getDist).
  PDF_obs <- getDist(obs_series, var, coef_obs)
  PDF_mod <- getDist(mod_series[, 1:ny_obs, drop = FALSE], var, coef_mod)

  # 3) Apply the CDF of the modeled data,
  #    evaluated with the statistics of the modeled data in the historical
  #    period, to the modeled data (getCDF).
  #    Eq. 1 of Cannon et al. (2015).
  Taot <- getCDF(PDF_mod, mod_series, coef_mod)

  # 4) Apply the inverse CDF of the observed
  #    data, evaluated with the statistics of the observed data in the
  #    historical period, to the probabilities obtained from 3) (getCDFinv).
  #    Eq. 1 of Cannon et al. (2015).
  QM_series <- getCDFinv(PDF_obs, Taot, coef_obs)
  QM_series <- matrix(QM_series)
  if (var == 1) {
    QM_series[QM_series < pp_threshold] <- 0
  }
  return(QM_series)
}
