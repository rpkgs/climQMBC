#' Univariate bias correction via quantile delta mapping
#'
#' Univariate bias correction based on the quantile delta mapping `QDM`
#' version of the quantile mapping algorithm from Cannon et al. (2015).
#' `QDM` constrains model-projected changes in quantiles to be preserved
#' following bias correction by quantile mapping.
#'
#' @param o.c vector of observed samples during the calibration period.
#' @param m.c vector of model outputs during the calibration period.
#' @param m.p vector of model outputs during the projected period.
#' @param ratio logical value indicating if samples are of a ratio quantity
#' (e.g., precipitation).
#' @param trace numeric value indicating the threshold below which values of a
#' ratio quantity (e.g., `ratio=TRUE`) should be considered exact zeros.
#' @param trace.calc numeric value of a threshold used internally when handling
#' of exact zeros; defaults to one half of `trace`.
#' @param jitter.factor optional strength of jittering to be applied when
#' quantities are quantized.
#' @param n.tau number of quantiles used in the quantile mapping; `NULL`
#' equals the length of the `m.p` series.
#' @param ratio.max numeric value indicating the maximum proportional change
#' allowed for ratio quantities below the `ratio.max.trace` threshold.
#' @param ratio.max.trace numeric value of a trace threshold used to constrain
#' the proportional change in ratio quantities to `ratio.max`; defaults to
#' ten times `trace`.
#' @param ECBC logical value indicating whether `mhat.p` outputs should be
#' ordered according to `o.c` ranks, i.e., as in the empirical copula-bias
#' correction (ECBC) algorithm.
#' @param ties method used to handle ties when calculating ordinal ranks.
#' @param subsample use `subsample` draws of size `n.tau` to
#' calculate empirical quantiles; if `NULL`, calculate normally.
#' @param pp.type type of plotting position used in `quantile`.
#' 
#' 
#' @return a list of with elements consisting of:
#'
#' - `mhat.c`: vector of bias corrected `m.c` values for the calibration period.
#'
#' - `mhat.p`: vector of bias corrected `m.p` values for the projection period.

#' @references
#' 1. Cannon, A.J., S.R. Sobie, and T.Q. Murdock, 2015. Bias correction of
#' simulated precipitation by quantile mapping: How well do methods preserve
#' relative changes in quantiles and extremes? Journal of Climate, 28:6938-6959.
#' doi:10.1175/JCLI-D-14-00754.1
#'
#' @seealso `[MBCp], [MBCr], [MRS], [escore]`
#' @export
MBC_QDM <- function(
    o.c, m.c, m.p,
    ratio = FALSE, trace = 0.05, trace.calc = 0.5 * trace,
    jitter.factor = 0, n.tau = length(m.p),
    ratio.max = 2, ratio.max.trace = 10 * trace,
    ECBC = FALSE, ties = "first", subsample = NULL, pp.type = 7) {
  if (jitter.factor == 0 &&
    (length(unique(o.c)) == 1 ||
      length(unique(m.c)) == 1 ||
      length(unique(m.p)) == 1)) {
    jitter.factor <- sqrt(.Machine$double.eps)
  }
  if (jitter.factor > 0) {
    o.c <- jitter(o.c, jitter.factor)
    m.c <- jitter(m.c, jitter.factor)
    m.p <- jitter(m.p, jitter.factor)
  }
  # For ratio data, treat exact zeros as left censored values less than trace.calc
  if (ratio) {
    o.c %<>% tidy_censored(trace.calc)
    m.c %<>% tidy_censored(trace.calc)
    m.p %<>% tidy_censored(trace.calc)
  }
  # Calculate empirical quantiles
  n <- length(m.p)
  # if (is.null(n.tau)) n.tau <- n
  tau <- seq(0, 1, length = n.tau)
  if (!is.null(subsample)) {
    quant.o.c <- rowMeans(apply(
      replicate(subsample, sample(o.c, size = length(tau))),
      2, quantile,
      probs = tau, type = pp.type
    ))
    quant.m.c <- rowMeans(apply(
      replicate(subsample, sample(m.c, size = length(tau))),
      2, quantile,
      probs = tau, type = pp.type
    ))
    quant.m.p <- rowMeans(apply(
      replicate(subsample, sample(m.p, size = length(tau))),
      2, quantile,
      probs = tau, type = pp.type
    ))
  } else {
    quant.o.c <- quantile(o.c, tau, type = pp.type)
    quant.m.c <- quantile(m.c, tau, type = pp.type)
    quant.m.p <- quantile(m.p, tau, type = pp.type)
  }

  # Apply quantile delta mapping bias correction
  tau.m.p <- approx(quant.m.p, tau, m.p, rule = 2)$y # tau of all values, ecdf(m.p)(m.p)

  # 这里划分成30年的滑动时间段，更合理
  m.p_c <- inverse_cdf(tau, quant.m.c, tau.m.p) # for increasing ratio
  m.p_o <- inverse_cdf(tau, quant.o.c, pout = tau.m.p)

  if (ratio) {
    # 对微量降水进行处理，设定最大变化倍数为`ratio.max`
    delta.m <- m.p / m.p_c # change ratio, Cannon 2015, Eq. 4
    delta.m[(delta.m > ratio.max) & (m.p_c < ratio.max.trace)] <- ratio.max

    mhat.p <- m.p_o * delta.m
  } else {
    delta.m <- m.p - m.p_c
    mhat.p <- m.p_o + delta.m
  }

  # 对训练期进行处理，这里采用的分段线性插值
  mhat.c <- approx(quant.m.c, quant.o.c, m.c, rule = 2)$y

  if (ratio) {
    # For ratio data, set values less than trace to zero
    # 对微量降水进行处理
    mhat.c[mhat.c < trace] <- 0
    mhat.p[mhat.p < trace] <- 0
  }

  if (ECBC) {
    # empirical copula coupling/Schaake shuffle
    if (length(mhat.p) == length(o.c)) {
      mhat.p <- sort(mhat.p)[rank(o.c, ties.method = ties)]
    } else {
      stop("Schaake shuffle failed due to incompatible lengths")
    }
  }
  list(mhat.c = mhat.c, mhat.p = mhat.p)
}

inverse_cdf <- function(p, x, pout) {
  approx(p, x, pout, rule = 2)$y
}

tidy_censored <- function(x, trace.calc = 0.025) {
  x[x < trace.calc] <- runif(
    sum(x < trace.calc),
    min = .Machine$double.eps, max = trace.calc
  )
  x
}
