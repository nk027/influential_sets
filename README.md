# Influence – Sensitivity to Influential Sets

This repository contains:

1. An R package, `influence`, to study sensitivity of linear regression (OLS and 2SLS/IV) to *influential sets* of observations.
2. A replication bundle for the paper *Hidden in plain sight: influential sets in linear regression* ([Kuschnig, Zens, and Crespo Cuaresma](https://www.kuschnig.eu/files/wp_influential-sets_wip.pdf)) in `paper/`.

## Usage

Install the package from this repository and load it:
```r
devtools::install_github("nk027/influential_sets")
```

### Quick Start
For a quick start, try out the following:
```
library("influence") # Load the package

# 1. reproduce the toy example -----
set.seed(0463) # Simulate some random (seeded) data
N <- 54 # We'll add 3 and 3 outliers to the end
x <- c(rnorm(N), rnorm(3, 6, 0.25), rnorm(3, 8, 0.25))
y <- c( # The outliers differ in effect
  x[seq(N)] * -0.5 + rnorm(N, 0, 1),
  x[seq(N + 1, N + 3)] * 0.1 + rnorm(3, 0, 0.1),
  x[seq(N + 4, N + 6)] * 0.4 + rnorm(3, 0, 0.1)
)
mdl <- lm(y ~ x - 1) # Fit a linear model (w/o intercept) ---
plot(x, y); abline(mdl, col = "darkgray", lwd = 2)
# Find the influential sets ---
mdl_sens <- sens(mdl) # Turns the first coef negative by default
mdl_sens$influence$id[1:7] # Adaptive results
mdl_sens$initial$id[1:7] # Initial approximation

# 2. reproduce the ruggedness application -----
data <- read.csv("https://short.wu.ac.at/rugged_data") # Read from GitHub
data$diamonds <- with(data, gemstones / (land_area / 100))
data <- data[!is.na(data$rgdppc_2000), ] # Remove NAs beforehand
# Fit the baseline model ---
mdl <- lm(log(rgdppc_2000) ~ rugged * cont_africa +
  diamonds * cont_africa + soil * cont_africa +
  tropical * cont_africa + dist_coast * cont_africa, data = data)
summary(mdl) # Summary and non-robust SE
# lmtest::coeftest(mdl, vcov = sandwich::vcovHC(mdl, "HC1")) # Robust SE
# Assess influential sets ---
mdl_sens <- sens(mdl, # Target t of rugged:cont_africa at position 8
  lambda = set_lambda("tstat", pos = 8, sign = sign(coef(mdl)[8])),
  options = set_options("all"), # No approximations
  cluster = seq_along(mdl$fitted) # For HC1 SE
)
plot(mdl_sens, threshold = qnorm(.975)) # Plot the path
summary(mdl_sens) # Get the summary of set sizes
data[mdl_sens$influence$id[seq(5)], ] # Check out the top 5 set
```

### Functionality

The core functions are:
- `infl()` computes influence diagnostics for `lm` and `ivreg` objects.
- `sens()` computes sensitivity paths under sequential removals (adaptive and non-adaptive).
- `init()` computes an initial approximation used as a baseline/comparator.
- `goal()` searches for the smallest set size needed to hit a target (experimental).
- `set_lambda()`, `set_target()`, `set_compute()`, `set_options()` define targets and computational settings.

## Replication Package

Please see [`papers/README.md`](papers/README.md).
