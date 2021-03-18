# AccuCor2 - Isotope Natural Abundance Correction for dual-isotope data

AccuCor2 is an isotope natural abundance correction algorithm for dual-isotope experimental data. AccuCor2 supports 13C-2H and 13C-15N isotope correction.

AccuCor2 accepts Excel (.xls and .xlsx) files.

## Installation

Download [R](https://www.r-project.org/) to install AccuCor2.

```R
install.packages("devtools")
library(devtools)
devtools::install_github("wangyujue23/AccuCor2")
```

## Usage

```R
library(accucor2)

# Input file (example file included)
# Or use your own: input_file <- "/path/to/my/datafile.xlsx"

input_file <- system.file("extdata", "CH_Compound_Test.xlsx", package = "accucor2")


# Output is written to [input_file]_corrected.xlsx by default
# Be sure to specify the appropriate resolution.
# For Exactive, the resolution is 100000, defined at 200 Mw

corrected <- natural_abundance_correction(
  path = input_file,
  resolution = 100000)


# The results are also returned as a named list of dataframes for further processing in R
# "Original", "Corrected", "Normalized", "PoolBeforeDF", "PoolAfterDF"

corrected


# Purity is set to 0.99 for C and N
# Be sure to specify purity if your samples differ
carbon_corrected <- natural_abundance_correction(
  path = input_file,
  resolution = 100000,
  C_purity = 0.97, 
  H_purity = 0.97)

```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[MIT](https://choosealicense.com/licenses/mit/)
