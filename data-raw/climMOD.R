# devtools::load_all()
library(climQMBC2)
library(glue)
library(data.table)
library(magrittr)

# Load observed and model data.
# - `temperature`   : `var = 0`
# - `precipitation` : `var = 1`

get_data <- function(var = 0) {
  if (var == 1) {
    var_txt <- "pp"
  } else {
    var_txt <- "tmax"
  }

  f <- glue("data-raw/obs_{var_txt}.csv") %>% system.file(package = "climQMBC2")
  obs <- fread(f)$V1 %>% matrix()

  f <- glue("data-raw/mod_{var_txt}.csv") %>% system.file(package = "climQMBC2")
  mod <- fread(f)$V1 %>% matrix()

  listk(
    frq = "M",
    var,
    var_txt,
    obs, mod
  )
}

clim_Prcp <- get_data(var = 1)
clim_Temp <- get_data(var = 0)

library(usethis)
# > monthly data
use_data(clim_Prcp, overwrite = TRUE)
use_data(clim_Temp, overwrite = TRUE)
