NULL
#' Get path to `SGDr` example
#'
#' SGDr comes bundled with a number of sample files in its `inst/extdata`
#' directory. This function make them easy to access
#'
#' @param file Name of file. If `NULL`, the example files will be listed.
#' @examples
#' SGDr_example()
#' SGDr_example("test_data.csv")
#' @export
SGDr_example <- function(file = NULL) {
  if (is.null(file)) {
    dir(system.file("extdata", package = "SGDr"))
  } else {
    system.file("extdata", file, package = "SGDr", mustWork = TRUE)
  }
}