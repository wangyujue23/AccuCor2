write_output <- function(dataframe_list, path, filetype = NULL, ...) {
  if(identical(FALSE, path)) {
    return(path)
  }
  
  if (is.null(path)) {
    return(path)
  }
  
  if (is.null(filetype)) {
    filetype = tools::file_ext(path)
  }
  
  if (filetype %in% c("xls", "xlsx")) {
    output_path = paste(tools::file_path_sans_ext(path), "_corrected.", filetype, sep="")
    message(paste("Output written to: '", output_path, "'", sep = ""))
    writexl::write_xlsx(dataframe_list, output_path, ...)
  } else if (filetype %in% c("csv", "txt", "tsv")) {
    if (filetype %in% c("csv", "txt")) {
      write_func <- readr::write_csv
    } else if (filetype %in% c("tsv")) {
      write_func <- readr::write_tsv
    }
    for (sheet in names(dataframe_list)) {
      if (tolower(sheet) == "original") {
        next()
      }
      sheet_path <- paste(tools::file_path_sans_ext(path), "_", tolower(sheet),
                          ".", filetype, sep="")
      message(paste("Output written to: '", sheet_path, "'", sep = ""))
      write_func(dataframe_list[[sheet]], sheet_path, ...)
    }
  } else {
    stop(paste("Unsupported output filetype: '", filetype, "'", sep = ""))
  }
  return(path)
}