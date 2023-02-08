#' formatQM
#'
#' This function formats the inputs and gets basic statistics for the different
#' Quantile Mapping (QM, DQM, QDM, UQM and SDM) methods available in the
#' climQMBC package. If monthly data is specified, the input series will be
#' reshaped to a matrix of 12 rows and several columns equal to the number of
#' years of each series. If annual data is specified, the input is reshaped to a
#' row vector with same entries as the input series. For precipitation,
#' physically null values (values below `pp_threshold`) are replaced by random
#' positive values below pp_factor.
#'
#' @param obs,mod A column vector of monthly or annual observed (modeled) data
#' (temperature or precipitation).
#'
#' If monthly frequency is specified, the length of this vector is 12 times the
#' number of years `[12 x ny]`.
#'
#' If annual frequency is specified, the length of this vector is equal to the
#' number of years `[ny]`.
#'
#' @param var A flag that identifies if data are temperature or precipitation.
#' This flag tells the getDist function if it has to discard distribution
#' functions that allow negative numbers, and if the terms in the correction
#' equations are multiplied/divided or added/subtracted. 
#' 
#' - `Temperature`  : `var = 0`;
#' - `Precipitation`: `var = 1`.
#' 
#' @param frq A string specifying if the input is annual or monthly data. If not
#' specified, it is set monthly as default. Monthly: `frq = 'M'`; Annual: `frq =
#' 'A'`
#'
#' @param pp_threshold A float indicating the threshold to consider physically
#' null precipitation values.
#'
#' @param pp_factor A float indicating the maximum value of the random values
#' that replace physically null precipitation values.

#' @return `list(obs, mod)`, `obs` and `mod` have the same structure:
#'
#' #### `data`:
#'
#' A column vector of monthly or annual observed data (temperature or
#'   precipitation). If monthly frequency is specified, the length of this
#'   vector is 12 times the number of observed years `[12, y_obs]`. If annual
#'   frequency is specified, the length of this vector is equal to the number of
#'   observed years `[1, y_obs]`.
#'
#' #### `coef`:
#'
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
#' #### `nyear`
#' 
#' length of year
#' 
#' @export
formatQM <- function(obs, mod, var, frq = "M", pp_threshold = 1, pp_factor = 1e-2) {
  # 0) Check if annually or monthly data is specified.
  # matrix: n_month, n_year
  nrow <- ifelse(frq == "A", 1, 12)
  nyear <- floor(length(obs) / nrow) # ncol
  n <- nyear * nrow
  obs <- obs[1:n]

  # 让obs和mod具有相同的长度
  nyear_mod <- floor(length(mod) / nrow) # ncol
  n_mod <- nyear_mod * nrow
  mod_bak <- mod[1:n_mod]
  mod <- mod_bak[1:n]
  
  process <- function(obs) {
    # 1. If variable is precipitation, replace low values with random values.
    if (var == 1) {
      bool_low_obs <- obs < pp_threshold
      obs[bool_low_obs] <- runif(sum(bool_low_obs)) * pp_factor
    }
    # 2. If monthly data is specified, reshape the input series to a matrix of
    #    12 rows and several columns equal to the number of years of each
    #    series. If annually data is specified, reshape the input to a row
    #    vector with same entries as the input series.
    obs_series <- matrix(obs, nrow = nrow, ncol = nyear)
    # 3. If monthly data is specified, get monthly mean, standard deviation,
    #   skewness, and log-skewness for the historical period of the observed
    #   and modeled series. If annually data is specified, get monthly mean,
    #   standard deviation, skewness, and log-skewness for the historical
    #   period of the observed and modeled series.
    mu <- apply(obs_series, 1, mean, na.rm = TRUE) # Mean
    sigma <- apply(obs_series, 1, sd, na.rm = TRUE) # Standard deviation
    skew <- apply(obs_series, 1, skewness2) # Skewness
    Ln_obs <- log(obs_series)
    Ln_obs[Im(Ln_obs) != 0] <- 0
    Ln_obs[!is.finite(Ln_obs)] <- log(0.01)
    skewy <- apply(Ln_obs, 1, skewness2) # Log-Skewness

    coef <- listk(mu, sigma, skew, skewy) %>% lapply(as.matrix)
    listk(data = obs_series, coef, nyear) #
  }

  l_obs <- process(obs)

  l_mod <- process(mod)
  l_mod$data = matrix(mod_bak, ncol = nyear_mod) # retrieve original data
  l_mod$nyear = nyear_mod

  list(obs = l_obs, mod = l_mod)
}
