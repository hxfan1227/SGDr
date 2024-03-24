#' Create a Parameters Object from the Parameters C++ Class
#'
#' Allows for the creation of a Parameters Object in _C++_ from _R_
#' using the _C++_ Parameters class.
#'
#' @param name Name of Parameters
#' @param age  Age of Parameters
#' @param male Is Parameters a Male?
#'
#' @return
#' A `Parameters` object from the _C++_ Parameters Class.
#'
#' @examples
#' ##################
#' ## Constructor
#'
#' # Construct new Parameters object called "ben"
#' ben = new(Parameters, name = "Ben", age = 26, male = TRUE)
#'
#' ##################
#' ## Getters
#'
#' ben$LikesBlue()
#'
#' ben$GetAge()
#'
#' ben$IsMale()
#'
#' ben$GetName()
#'
#' ben$GetFavoriteNumbers()
#' @name Parameters
#' @export Parameters

# ^^^^^^^^^^^^^^^^
# Export the "Parameters" C++ class by explicitly requesting Parameters be
# exported via roxygen2's export tag.
# Also, provide a name for the Rd file.


# Load the Rcpp module exposed with RCPP_MODULE( ... ) macro.
loadModule(module = "ParametersModule", TRUE)
loadModule(module = "ModelModule", TRUE)