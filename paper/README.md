
# Hidden in plain sight: influential sets in linear regression

The results of the paper depend on this package (nk027/influential_sets) as well as some additional data (provided) and packages. To install them, run 

```r
# Reproduce results
install.packages(c("ivreg", "AER", "sandwich", "lmtest"))
# Robust estimation
install.packages(c("MASS", "robustbase", "robust"))
# Tooling
install.packages(c("dplyr", "tidyr", "foreign"))
# Plots and tables
install.packages(c("ggplot2", "viridisLite", "kableExtra", "memisc"))
# For the spatial data and map
install.packages(c("sf", "tmap", "rworldmap"))
# For comparisons to Broderick, Giordano, and Meager's AMIP
# v0.0.0.9000 / commit 43b53f7 and recently 2d29fbb
devtools::install_github("rgiordan/zaminfluence", subdir = "zaminfluence")

## Notes

- Note that the scripts use `cairo_pdf` and some figures use custom fonts.
- The code and results are made available under te package's GPL3 license; data provided for replication is subject to original sources' licenses/terms; see the paper and source references.

## Structure

- `0_masking.R` produces outputs for the masking demonstration.
- `1_ruggedness.R` analyses the results of Nunn and Puga (2012).
- `2_mistrust.R` analyses the results of Nunn and Wantchekon (2011).
- `3_tsetse.R` analyses the results of Alsan (2015).
- `8_migration.R` analyses the results of Droller (2018).
- `8_poverty.R` analyses the results of Ravallion (2012) and Crespo Cuaresma et al. (2022).
- `9_amip.R` investigates properties of AMIP (Broderick et al., 2021)
- `9_microcredit.R` analyses the results of Broderick et al. (2021).
- `9_sim.R` conducts two simulation exercises to investigate algorithm performance in the presence of influential sets, and the technical sensitivities of ordinary and two-stage least squares.
