context("Dual isotope abundance correction")
library(accucor2)

read_expected <- function(file, sheet) {
  expected <- readxl::read_excel(path = file, sheet = sheet)
  expected <- dplyr::mutate_at(expected,
                               dplyr::vars(dplyr::ends_with("_Label")),
                               as.integer)
}

test_that("Carbon correction (Excel, simple format)", {
  resolution <- 100000
  input_file <- system.file("extdata", 
                            "CN_Compound_Test.xlsx", 
                            package = "accucor2")
  sheet_name <- "Sheet1"
  metabolite_list <-system.file("extdata", 
                                "Metabolite Formula and Charge Info.csv", 
                                package = "accucor2")
  
  corrected <- dual_correction(
    input_file,
    sheet_name,
    metabolite_list,
    "CN",
    Resolution = 70000)
  
  # Does not test contents, just that we got something of length 4
  # TODO Update to test contents of `corrected`
  expect_equal(length(corrected), 4)
  
  # 
  # expected_output <- list(
  #   "Original" = read_expected(
  #     system.file("extdata", "C_Sample_Input_Simple.xlsx", package = "accucor2"),
  #     sheet = 1),
  #   "Corrected" = read_expected(
  #     system.file("extdata", "C_Sample_Input_Simple_corrected.xlsx", package = "accucor2"),
  #     sheet = "Corrected"),
  #   "Normalized" = read_expected(
  #     system.file("extdata", "C_Sample_Input_Simple_corrected.xlsx", package = "accucor2"),
  #     sheet = "Normalized"),
  #   "PoolAfterDF" = read_expected(
  #     system.file("extdata", "C_Sample_Input_Simple_corrected.xlsx", package = "accucor2"),
  #     sheet = "PoolAfterDF")
  # )
  # 
  # expect_equal(corrected, expected_output)
})
