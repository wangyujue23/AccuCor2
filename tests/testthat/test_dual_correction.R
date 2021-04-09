context("Dual isotope abundance correction")
library(accucor2)

read_expected <- function(file, sheet) {
  expected <- readxl::read_excel(path = file, sheet = sheet)
  expected <- dplyr::mutate_at(expected,
                               dplyr::vars(dplyr::ends_with("_Label")),
                               as.integer)
}

test_that("CN correction", {
  resolution <- 70000
  input_file <- system.file("extdata", 
                            "CN_Compound_Test.xlsx", 
                            package = "accucor2")
  sheet_name <- "Sheet1"
  metabolite_list <-system.file("extdata", 
                                "metabolite_formula_and_charge_info.csv", 
                                package = "accucor2")
  
  corrected <- dual_correction(
    input_file,
    sheet_name,
    metabolite_list,
    "CN",
    Resolution = resolution)
  
  expected_output <- list(
    "Original" = read_expected(
      system.file("extdata", "CN_Compound_Test_corrected.xlsx", package = "accucor2"),
      sheet = "Original"),
    "Corrected" = read_expected(
      system.file("extdata", "CN_Compound_Test_corrected.xlsx", package = "accucor2"),
      sheet = "Corrected"),
    "Normalized" = read_expected(
      system.file("extdata", "CN_Compound_Test_corrected.xlsx", package = "accucor2"),
      sheet = "Normalized"),
    "Pool size" = read_expected(
      system.file("extdata", "CN_Compound_Test_corrected.xlsx", package = "accucor2"),
      sheet = "Pool size")
  )

  expect_equal(corrected, expected_output)
})

test_that("CH correction", {
  resolution <- 750000
  input_file <- system.file("extdata",
                            "CH_Compound_Test.xlsx",
                            package = "accucor2")
  sheet_name <- "Sheet1"
  metabolite_list <-system.file("extdata",
                                "metabolite_formula_and_charge_info.csv",
                                package = "accucor2")

  corrected <- dual_correction(
    input_file,
    sheet_name,
    metabolite_list,
    "CH",
    Resolution = resolution)

  expected_output <- list(
    "Original" = read_expected(
      system.file("extdata", "CH_Compound_Test_corrected.xlsx", package = "accucor2"),
      sheet = "Original"),
    "Corrected" = read_expected(
      system.file("extdata", "CH_Compound_Test_corrected.xlsx", package = "accucor2"),
      sheet = "Corrected"),
    "Normalized" = read_expected(
      system.file("extdata", "CH_Compound_Test_corrected.xlsx", package = "accucor2"),
      sheet = "Normalized"),
    "Pool size" = read_expected(
      system.file("extdata", "CH_Compound_Test_corrected.xlsx", package = "accucor2"),
      sheet = "Pool size")
  )

  expect_equal(corrected, expected_output)
})

# Other possible tests
# dual_correction("Compare with real data3.xlsx","raw","KNOWN_2018May.csv","CN",Resolution = 70000)
# dual_correction("CH_Compound_Test.xlsx","Sheet1","KNOWN_2018May.csv","CH",C13Purity = 1, H2N15Purity = 1, Resolution = 750000,NitrogenNaturalAbundace = c(0.99632, 0.00368),SulfurNaturalAbundace = c(0.9493, 0.0076, 0.0431))