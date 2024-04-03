#include <Rcpp.h>
#include "Model.h"

//' Estimate the SGD volume.
//' @param inputData A data frame containing the input data for the model.
//' @param calibratableParams A list of calibratable parameters.
//' @param constParams A list of constant parameters.
//' @param windowSize A integer indicating the period you want to average the recharge.
//' @description Estimate the SGD volume based on Strack (1976) analytical solution.
//' @details Estimate the SGD volume based on Strack (1976) analytical solution.
//' * The \code{inputData} should be a data.frame containing the required columns for the model (R, E0, Pumping, GWL, H2O_SB1, H2O_SB2).
//'   * \code{R} is a numeric vector of the daily precipitation (mm).
//'   * \code{E0} is a numeric vector of the daily potential evapotranspiration (mm).
//'   * \code{Pumping} is a numeric vector of the daily pumping rate (m3/d).
//'   * \code{GWL} is a numeric vector of the initial groundwater level (mm), Only the first value is used as the initial value for the SB1.
//'   * \code{H2O_SB1} is a numeric vector of the initial water level in the first soil bucket (mm). Only the first value is used as the initial value for the SB1.
//'   * \code{H2O_SB2} is a numeric vector of the inital water level in the second soil bucket (mm). Only the first value is used as the initial value for the SB2.
//' * \code{calibratableParams} should be a list of calibratable parameters. It's recommended to use \code{\link{json_to_paramter_list}} to create this list.
//' * \code{windowSize} should be a list of constant parameters. It's recommended to use \code{\link{json_to_paramter_list}} to create this list.
//' @return A data frame containing the estimated SGD volume.
//' @export 
// [[Rcpp::export]]
Rcpp::DataFrame estimate_sgd(const Rcpp::DataFrame& inputData, const Rcpp::List& calibratableParams, const Rcpp::List& constParams, int windowSize)
{
    Model model(inputData, calibratableParams, constParams);
    model.calc_recharge();
    model.calc_sgd(windowSize);
    return model.get_sgd_output();
}

//' Estimate the SGD volume with fixed hn and xn.
//' @inheritParams estimate_sgd  
//' @param hn A double indicating the prefered value for hn.
//' @param xn A double indicating the prefered value for xn.
//' @description Estimate the SGD volume based on Strack (1976) analytical solution with fixed hn and xn.
//' @details Estimate the SGD volume based on Strack (1976) analytical solution.
//' * The \code{inputData} should be a data.frame containing the required columns for the model (R, E0, Pumping, GWL, H2O_SB1, H2O_SB2).
//'   * \code{R} is a numeric vector of the daily precipitation (mm).
//'   * \code{E0} is a numeric vector of the daily potential evapotranspiration (mm).
//'   * \code{Pumping} is a numeric vector of the daily pumping rate (m3/d).
//'   * \code{GWL} is a numeric vector of the initial groundwater level (mm), Only the first value is used as the initial value for the SB1.
//'   * \code{H2O_SB1} is a numeric vector of the initial water level in the first soil bucket (mm). Only the first value is used as the initial value for the SB1.
//'   * \code{H2O_SB2} is a numeric vector of the inital water level in the second soil bucket (mm). Only the first value is used as the initial value for the SB2.
//' * \code{calibratableParams} should be a list of calibratable parameters. It's recommended to use \code{\link{json_to_paramter_list}} to create this list.
//' * \code{windowSize} should be a list of constant parameters. It's recommended to use \code{\link{json_to_paramter_list}} to create this list.
//' @return A data frame containing the estimated SGD volume.
//' @export 
// [[Rcpp::export]]
Rcpp::DataFrame estimate_sgd_fixed_hxn(const Rcpp::DataFrame& inputData, const Rcpp::List& calibratableParams, 
                                       const Rcpp::List& constParams, int windowSize, const double hn, const double xn)
{
  Model model(inputData, calibratableParams, constParams);
  model.calc_recharge();
  model.calc_sgd(windowSize, hn, xn);
  return model.get_sgd_output();
}
