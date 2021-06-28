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
    Resolution = resolution,
    output_base = FALSE)
  
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

test_that("CN correction XLSX output", {
  resolution <- 70000
  input_file <- system.file("extdata", 
                            "CN_Compound_Test.xlsx", 
                            package = "accucor2")
  sheet_name <- "Sheet1"
  metabolite_list <-system.file("extdata", 
                                "metabolite_formula_and_charge_info.csv", 
                                package = "accucor2")
  
  output_filetype = "xlsx"
  output_file_base <- tempfile(fileext = paste0(".", output_filetype))
  output_file <- paste0(tools::file_path_sans_ext(output_file_base),
                        "_corrected.",
                        output_filetype)
  
  
  expect_message(
    dual_correction(
      input_file,
      sheet_name,
      metabolite_list,
      "CN",
      Resolution = resolution,
      output_base = output_file_base
    ),
    "Output written to:.*_corrected\\.xlsx"
  )
  
  expected_output_file <- system.file("extdata", "CN_Compound_Test_corrected.xlsx", package = "accucor2")
  expected_output <- list(
    "Original" = read_expected(expected_output_file, sheet = "Original"),
    "Corrected" = read_expected(expected_output_file, sheet = "Corrected"),
    "Normalized" = read_expected(expected_output_file, sheet = "Normalized"),
    "Pool size" = read_expected(expected_output_file, sheet = "Pool size")
  )
  
  corrected <- list(
    "Original" = read_expected(output_file, sheet = "Original"),
    "Corrected" = read_expected(output_file, sheet = "Corrected"),
    "Normalized" = read_expected(output_file, sheet = "Normalized"),
    "Pool size" = read_expected(output_file, sheet = "Pool size")
  )
  
  expect_equal(corrected, expected_output)
})


test_that("CN correction csv output", {
  resolution <- 70000
  input_file <- system.file("extdata", 
                            "CN_Compound_Test.xlsx", 
                            package = "accucor2")
  sheet_name <- "Sheet1"
  metabolite_list <-system.file("extdata", 
                                "metabolite_formula_and_charge_info.csv", 
                                package = "accucor2")
  
  output_filetype = "csv"
  output_file_base <- tempfile(fileext = paste0(".", output_filetype))
  output_file <- paste0(tools::file_path_sans_ext(output_file_base),
                        "_corrected.",
                        output_filetype)
  
  
  expect_message(expect_message(
    expect_message(
      dual_correction(
        input_file,
        sheet_name,
        metabolite_list,
        "CN",
        Resolution = resolution,
        output_base = output_file_base,
        output_filetype = output_filetype
      ),
      "Output written to:.*_corrected\\.csv"
    ),
    "Output written to:.*_normalized\\.csv"
  ),
  "Output written to:.*_pool size\\.csv")
  
  for (sheet in list("Corrected", "Normalized", "Pool size")) {
    output_file_suffix <- paste0("_", tolower(sheet), ".", output_filetype)
    output_file <- paste0(tools::file_path_sans_ext(output_file_base), output_file_suffix)
    expect_snapshot_file(output_file, paste0("CN_Compound_Test", output_file_suffix), compare = compare_file_text)
  }
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
    Resolution = resolution,
    output_base = FALSE)

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