
#' Create a Bucket Object from the Bucket C++ Class
#'
#' Allows for the creation of a Bucket Object in _C++_ from _R_
#' using the _C++_ Bucket class.
#'
#' @param layer Soil layer number in bucket (-). In this model we only have 2 layers (1 and 2).
#' 


#' @return
#' A `Bucket` object from the _C++_ Bucket Class.
#'
#' @examples
#' ##################
#' ## Constructor
#'
#' # Construct new Bucket object called "bucket1"
#' 
#' bucketParams <- list(layer = 1, 
#'                      WP = 100, 
#'                      z = 200, 
#'                      FCmm = 300, 
#'                      SAT = 400, 
#'                      Swres = 500, 
#'                      Ksat = 600, 
#'                      n = 0.5, 
#'                      phi_soil = 0.3, 
#'                      Y = 0.2, 
#'                      Z = 0.1)
#' bucket1 <- new(Bucket, bucketParams)
#'
#' ##################
#' ## Getters
#'
#' bucket1$get_WP()
#' bucket1$get_z()
#' bucket1$get_FCmm()
#' bucket1$get_SAT()
#' bucket1$get_Swres()
#' bucket1$get_Ksat()
#' bucket1$get_n()
#' bucket1$get_phi_soil()
#' bucket1$get_Y()
#' bucket1$get_Z()
#' @name Parameters
#' @export Parameters

# ^^^^^^^^^^^^^^^^
# Export the "Bucket" C++ class by explicitly requesting Bucket be
# exported via roxygen2's export tag.
# Also, provide a name for the Rd file.


# Load the Rcpp module exposed with RCPP_MODULE( ... ) macro.
Rcpp::loadModule(module = "ParametersModule", TRUE)

NULL

#' @name Model
#' @export Model
Rcpp::loadModule(module = "ModelModule", TRUE)
