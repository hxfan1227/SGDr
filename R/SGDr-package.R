#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @useDynLib SGDr, .registration = TRUE
NULL

# Expose the Bucket class
#' @export Bucket
Rcpp::loadModule(module = "BucketModule", TRUE)