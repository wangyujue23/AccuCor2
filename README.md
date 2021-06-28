# AccuCor2 - Isotope Natural Abundance Correction for dual-isotope data

<!-- badges: start -->
[![R-CMD-check](https://github.com/wangyujue23/AccuCor2/workflows/R-CMD-check/badge.svg)](https://github.com/wangyujue23/AccuCor2/actions)
<!-- badges: end -->

AccuCor2 is an isotope natural abundance correction algorithm for dual-isotope experimental data. AccuCor2 supports 13C-2H and 13C-15N isotope correction.

AccuCor2 accepts Excel (.xls and .xlsx) files.

## Installation

Download [R](https://www.r-project.org/) to install accucor2.

```R
# Please make sure all the required packages are updated.
install.packages("devtools")
library(devtools)
devtools::install_github("wangyujue23/AccuCor2",INSTALL_opts = "--no-multiarch")
```

## Usage

```R
library(accucor2)

# Example for 13C-15N tracer data:
# Input file and the Sheet name (example file included)
# Or use your own: input_file <- "/path/to/my/datafile.xlsx"
# You may use the example files in the "inst/extdata" folder as a template to format your own data.

input_file <- system.file("extdata", "CN_Compound_Test.xlsx", package = "accucor2")
sheet_name <- "Sheet1"
metabolite_list <-system.file("extdata", "metabolite_formula_and_charge_info.csv", package = "accucor2")

# Be sure to specify the appropriate resolution.
# For Exactive, the resolution is 70000, defined at 200 Mw
# Output is written into "corrected_[time].xlsx" under the directory folder by default (a message will show you where it is)
# The results are also returned as a named list of dataframes for further processing in R
# "Original", "Corrected", "Normalized", "Pool size"

corrected <- dual_correction(input_file,sheet_name,metabolite_list,"CN",Resolution = 70000)

corrected


# Example for 13C-2H tracer data:


input_file <- system.file("extdata", "CH_Compound_Test.xlsx", package = "accucor2")
sheet_name <- "Sheet1"
metabolite_list <-system.file("extdata", "Metabolite Formula and Charge Info.csv", package = "accucor2")

# The default value for tracer C/H/N isotopic purity is 0.99, you could change them by specify the value of C13Purity and H2N15Purity.
# You could also specify the natual abundance value. In this example, the NitrogenNaturalAbundace (14N, 15N) and SulfurNaturalAbundace (32S, 33S, 34S) are specified.

corrected <- dual_correction(input_file,sheet_name,metabolite_list,"CH",C13Purity = 1, H2N15Purity = 1, Resolution = 750000,NitrogenNaturalAbundance = c(0.99632, 0.00368),SulfurNaturalAbundance = c(0.9493, 0.0076, 0.0431))

corrected

```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License

## Potential problem and solution
Error:
```R
Warning in install.packages :
  download of package ‘utf8’ failed
```
Solution:
```R
install.packages("utf8")
```

Error:
```R
Loading required package: usethis
```
Solution:
```R
install.packages("usethis")
library(usethis)
```
