# devtools::load_all()
library(climQMBC2)
library(glue)
library(data.table)
library(magrittr)

# Load observed and model data.
# - `temperature`   : `var = 0`
# - `precipitation` : `var = 1`
var <- 0

if (var == 1) {
  agg <- sum
  var_txt <- "pp"
} else {
  var_txt <- "tmax"
  agg <- mean
}

f <- glue("data-raw/obs_{var_txt}.csv") %>% system.file(package = "climQMBC2")
obs <- fread(f)$V1 %>% matrix()

f <- glue("data-raw/mod_{var_txt}.csv") %>% system.file(package = "climQMBC2")
mod <- fread(f)$V1 %>% matrix()


climMOD <- listk(
  frq = "M",
  var,
  var_txt,
  obs, mod
)
use_data(climMOD, overwrite = TRUE)
# > monthly data
