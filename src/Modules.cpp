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

RCPP_MODULE(ModelModule)
{
    using namespace Rcpp;
    function("estimate_sgd", &estimate_sgd);

    // Rcpp::class_<Model>("Model")
    //     .constructor<Rcpp::DataFrame, Rcpp::List, Rcpp::List>()
    //     .method("calc_recharge", &Model::calc_recharge)
    //     .method("get_recharge_output", &Model::get_recharge_output)
    //     .method("get_sgd_output", &Model::get_sgd_output)
    //     .method("get_inputData", &Model::get_inputData)
    //     .method("get_SW", &Model::get_SW)
    //     .method("get_S", &Model::get_S)
    //     .method("get_Ia", &Model::get_Ia)
    //     .method("get_CN", &Model::get_CN)
    //     .method("get_Qsurf", &Model::get_Qsurf)
    //     .method("get_MC_SB1", &Model::get_MC_SB1)
    //     .method("get_DoS_SB1", &Model::get_DoS_SB1)
    //     .method("get_finf_SB1", &Model::get_finf_SB1)
    //     .method("get_finfla_SB1", &Model::get_finfla_SB1)
    //     .method("get_H2O1_SB1", &Model::get_H2O1_SB1)
    //     .method("get_IntercH2O", &Model::get_IntercH2O)
    //     .method("get_CanopyH2O", &Model::get_CanopyH2O)
    //     .method("get_E0_Int", &Model::get_E0_Int)
    //     .method("get_YE", &Model::get_YE)
    //     .method("get_E", &Model::get_E)
    //     .method("get_H2O2_SB1", &Model::get_H2O2_SB1)
    //     .method("get_QSE", &Model::get_QSE)
    //     .method("get_SWexcess", &Model::get_SWexcess)
    //     .method("get_Kunsat", &Model::get_Kunsat)
    //     .method("get_Kunsat2", &Model::get_Kunsat2)
    //     .method("get_TT", &Model::get_TT)
    //     .method("get_TT2", &Model::get_TT2)
    //     .method("get_Wdperc", &Model::get_Wdperc)
    //     .method("get_Wperc", &Model::get_Wperc)
    //     .method("get_Wperc2", &Model::get_Wperc2)
    //     .method("get_H2O3_SB1", &Model::get_H2O3_SB1)
    //     .method("get_H2O3_SB2", &Model::get_H2O3_SB2)
    //     .method("get_MC_SB2", &Model::get_MC_SB2)
    //     .method("get_DoS_SB2", &Model::get_DoS_SB2)
    //     .method("get_H2O1_SB2", &Model::get_H2O1_SB2)
    //     .method("get_YE2", &Model::get_YE2)
    //     .method("get_Edd", &Model::get_Edd)
    //     .method("get_H2O2_SB2", &Model::get_H2O2_SB2)
    //     .method("get_Sw", &Model::get_Sw)
    //     .method("get_Sw2", &Model::get_Sw2)
    //     .method("get_SWexcess2", &Model::get_SWexcess2)
    //     .method("get_Wrechg", &Model::get_Wrechg)
    //     .method("get_H2O1_AQ", &Model::get_H2O1_AQ)
    //     .method("get_WrechgAve", &Model::get_WrechgAve)
    //     .method("get_H2O2_AQ", &Model::get_H2O2_AQ)
    //     .method("get_H2O3_AQ", &Model::get_H2O3_AQ)
    //     .method("calc_sgd", &Model::calc_sgd)
    //     .method("get_SGD1", &Model::get_SGD1)
    //     .method("get_SGD2", &Model::get_SGD2)
    //     .method("get_xn", &Model::get_xn)
    //     .method("get_xn1", &Model::get_xn1)
    //     .method("get_xn2", &Model::get_xn2)
    //     .method("get_hn", &Model::get_hn)
    //     .method("get_hn1", &Model::get_hn1)
    //     .method("get_hn2", &Model::get_hn2)
    //     .method("get_M1", &Model::get_M1)
    //     .method("get_M2", &Model::get_M2)
    //     .method("get_xT1", &Model::get_xT1)
    //     .method("get_xT2", &Model::get_xT2)
    //     .method("get_SGD", &Model::get_SGD)
    //     .method("get_SGDdrop", &Model::get_SGDdrop)
    //     .method("get_PumpingDrop", &Model::get_PumpingDrop)
    //     .method("get_dxT", &Model::get_dxT)
    //     .method("get_xT3", &Model::get_xT3)
    //     .method("get_SWvol", &Model::get_SWvol)
    //     .method("get_FWLdrop", &Model::get_FWLdrop)
    //     .method("update_parameters", &Model::update_parameters)
    //     .method("get_all_params_list", &Model::get_all_params_list);
}