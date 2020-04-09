library('shiny')

runCHARTS <- function(data_dir) {
  .GlobalEnv$DATA_DIR <- data_dir
  runApp("src/shiny")
}

Args = commandArgs(trailingOnly=TRUE)

data_dir <- Args[1]
runCHARTS(data_dir)
