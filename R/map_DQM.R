#' @description
#' - Detrended Quantile Mapping (DQM): Cannon et al. (2015).
#' @rdname map_QM
#' 
#' @export
map_DQM <- function(obs, mod, var, frq = "M", pp_threshold = 1, pp_factor = 1e-2) {
  l <- formatQM(obs, mod, var, frq, pp_threshold, pp_factor)

  ny_obs <- l$obs$nyear
  ny_mod <- l$mod$nyear
  obs_series <- l$obs$data
  mod_series <- l$mod$data
  nfrq <- dim(mod_series)[1]
  
  mu_obs     <- l$obs$coef$mu
  std_obs    <- l$obs$coef$sigma
  skew_obs   <- l$obs$coef$skew
  skewny_obs <- l$obs$coef$skewy

  mu_mod     <- l$mod$coef$mu
  std_mod    <- l$mod$coef$sigma
  skew_mod   <- l$mod$coef$skew
  skewny_mod <- l$mod$coef$skewy

  ## 2) Assign a probability distribution function to each month for the observed
  #    and modeled data in the historical period. If annual frequency is
  #    specified, this is applied to the complete historical period (getDist).
  PDF_obs <- getDist(obs_series, mu_obs, std_obs, skew_obs, skewny_obs, var)
  PDF_mod <- getDist(mod_series[, 1:ny_obs, drop=FALSE], mu_mod, std_mod, skew_mod, skewny_mod, var)

  ## 3) Extract the long-term trend from the modeled data: 

  # a) Get the monthly mean of the historical period. If annually frequency is
  # specified, this is applied to the complete period).
  xbarmh <- apply(mod_series[, 1:ny_obs, drop = FALSE], 1, mean, na.rm = TRUE)

  # b) Get the monthly mean of each projected period. If annually frequency is
  # specified, this is applied to the complete period).
  xbarmp <- matrix(0, nfrq, (ny_mod - ny_obs))

  for (j in 1:(ny_mod - ny_obs)) {
    xbarmp[, j] <- apply(mod_series[, (j + 1):(ny_obs + j), drop = FALSE], 1, mean)
  }

  # 4)Compute the linear scaled values (value in square brackets in Eq. 2
  # of Cannon et al. (2015)).
  LS <- matrix(0, nfrq, (ny_mod - ny_obs))
  
  scale_mult <- function(x_mp, xbar_mp, xbar_mh) {
    x_mp * xbar_mh / xbar_mp # adjust mean to his
  }
  scale_add <- function(x_mp, xbar_mp, xbar_mh) {
    x_mp + xbar_mh - xbar_mp # adjust mean to his
  }
  
  fun_scale = ifelse(var == 1, scale_mult, scale_add)
  for (m in 1:nfrq) {
    LS[m, ] <- fun_scale(mod_series[m, (ny_obs + 1):ny_mod], xbarmp[m, ], xbarmh[m])
  }

  # 5) Apply the CDF of the modeled data, evaluated with the statistics of the
  #    modeled period, to the future modeled data (getCDF). Eq.2, Cannon 2015.
  Taot <- getCDF(PDF_mod, LS, mu_mod, std_mod, skew_mod, skewny_mod)

  # 6) Apply the inverse CDF of the observed data, evaluated with the statistics
  #    of the observed data in the historical period, to the probabilities
  #    obtained from 5) (getCDFinv). Eq.2, Cannon 2015.
  DQM_LS <- getCDFinv(PDF_obs, Taot, mu_obs, std_obs, skew_obs, skewny_obs)

  # 7) Reimpose the trend to the values obtained in 6). Eq.2, Cannon 2015.
  ## ReverseScale
  xbar_mh = matrix(rep(xbarmh, ny_mod - ny_obs), ncol = ny_mod - ny_obs)
  if (var == 1) { 
    DQM <- DQM_LS * (xbarmp / xbar_mh)
  } else { 
    DQM <- DQM_LS + (xbarmp - xbar_mh)
  }
  DQM <- matrix(DQM)

  # 8) Perform QM for the historical period.
  mod_h <- mod_series[, 1:ny_obs, drop=FALSE]
  QM_series <- QM(obs, mod_h, var, frq)
  DQM_series <- c(QM_series, DQM)
  if (var == 1) {
    DQM_series[DQM_series < pp_threshold] <- 0
  }
  return(DQM_series)
}
