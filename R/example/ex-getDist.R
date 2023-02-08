# temperature data
var = 0
input <- clim_Temp

l <- with(input, formatQM(obs, mod, var))
print(str(l))

l_obs <- l$obs
l_mod <- l$mod

DIST_mod <- with(l_mod, getDist(data, var, coef))

cdf_mod <- with(l_mod, getCDF(DIST_mod, data, coef))

DIST_obs <- with(l_obs, getDist(data, var, coef))
QM_series <- getCDFinv(DIST_obs, cdf_mod, l_obs$coef)

r <- list(dist = DIST_mod, cdf = cdf_mod)
print(str(r))
