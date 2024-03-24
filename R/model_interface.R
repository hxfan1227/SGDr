#' @import jsonlite
NULL
#' @name create_sgd_model
#' @aliases R_model_interface
#' @title Gernerate a SGD model in R
#' @rdname R_model_interface
#' @description This function allows for the creation of a `Model` object in C++ from R.
#' @param input_data A data frame containing the input data for the model
#' @param calibratable_parameters A list of calibratable parameters
#' @param const_paramters A list of constant parameters
#' @return A `Model` object from the _C++_ Model Class.
#' @details This function allows for the creation of a \code{Model} object in C++ from _R.
#' * The \code{input_data} should be a data.frame containing the required columns for the model (R, E0, Pumping, GWL, H2O_SB1, H2O_SB2).
#'      * \code{R} is a numeric vector of the daily precipitation (mm). 
#'      * \code{E0} is a numeric vector of the daily potential evapotranspiration (mm). 
#'      * \code{Pumping} is a numeric vector of the daily pumping rate (m3/d).
#'      * \code{GWL} is a numeric vector of the initial groundwater level (mm), Only the first value is used as the initial value for the SB1. 
#'      * \code{H2O_SB1} is a numeric vector of the initial water level in the first soil bucket (mm). Only the first value is used as the initial value for the SB1. 
#'      * \code{H2O_SB2} is a numeric vector of the inital water level in the second soil bucket (mm). Only the first value is used as the initial value for the SB2. 
#' * \code{calibratable_parameters} should be a list of calibratable parameters. It's recommended to use \code{\link{json_to_paramter_list}} to create this list.
#' * \code{const_paramters} should be a list of constant parameters. It's recommended to use \code{\link{json_to_paramter_list}} to create this list.
#' @export
create_sgd_model <- function(input_data, calibratable_parameters, const_paramters) {
    # stop if input_data does not contain the required columns for the model (R, E0, Pumping, GWL, H2O_SB1, H2O_SB2)
    if (!all(c("R", "E0", "Pumping", "GWL", "H2O_SB1", "H2O_SB2") %in% colnames(input_data))) {
        stop("Input data does not contain the required columns for the model (R, E0, Pumping, GWL, H2O_SB1, H2O_SB2)")
    }
    model <- new(Model, input_data, calibratable_parameters, const_paramters)
    return(model)
}

NULL 
#' @name json_to_paramter_list
#' @aliases json_to_paramter_list
#' @title Convert a JSON file to a list of parameters
#' @rdname json_to_paramter_list
#' @description This function converts a JSON file to a list of parameters.
#' @param json_file A JSON file containing the parameters
#' @return A list of parameters (invisible)
#' @export

json_to_paramter_list <- function(json_file) {
    invisible(as.relistable(jsonlite::fromJSON(json_file)))
}

NULL 
#' @name paramter_vec_to_list
#' @aliases paramter_vec_to_list
#' @title Convert a vector of parameters to a list of parameters (useful for calibrating parameters)
#' @rdname paramter_vec_to_list
#' @description This function converts a vector of parameters to a list of parameters.
#' @param parameter_vec A vector of parameters
#' @param prameter_skeleton A list of parameters to use as a skeleton for the new list
#' @return A list of parameters (invisible)
#' @export

paramter_vec_to_list <- function(parameter_vec, prameter_skeleton) {
    invisible(relist(parameter_vec, skeleton = prameter_skeleton))
}
