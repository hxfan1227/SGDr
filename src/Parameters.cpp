#include "Paramters.h"

Parameters::Parameters(Rcpp::List calibratableParams, Rcpp::List constParams) : bucket1(Rcpp::as<Rcpp::List>(calibratableParams["bucket1"]), constParams),
                                                                                bucket2(Rcpp::as<Rcpp::List>(calibratableParams["bucket2"]), constParams),
                                                                                curveNumber(Rcpp::as<Rcpp::List>(calibratableParams["curveNumber"]), bucket1, bucket2),
                                                                                aquifer(Rcpp::as<Rcpp::List>(calibratableParams["aquifer"])),
                                                                                constParameter(constParams)
{
}

Bucket Parameters::get_bucket1() { return bucket1; }

Bucket Parameters::get_bucket2() { return bucket2; }

CurveNumber Parameters::get_curveNumber() { return curveNumber; }

Aquifer Parameters::get_aquifer() { return aquifer; }

ConstParameter Parameters::get_constParameter() { return constParameter; }

Rcpp::List Parameters::get_all_params_list()
{
    return Rcpp::List::create(
        Rcpp::Named("bucket1") = bucket1.get_all_params_list(),
        Rcpp::Named("bucket2") = bucket2.get_all_params_list(),
        Rcpp::Named("curveNumber") = curveNumber.get_all_params_list(),
        Rcpp::Named("aquifer") = aquifer.get_all_params_list());
}

Rcpp::List Parameters::get_const_params_list()
{
    return constParameter.get_all_params_list();
}

Rcpp::List Parameters::get_calibratable_params_list()
{
    return Rcpp::List::create(
        Rcpp::Named("bucket1") = bucket1.get_calibratable_params_list(),
        Rcpp::Named("bucket2") = bucket2.get_calibratable_params_list(),
        Rcpp::Named("curveNumber") = curveNumber.get_calibratable_params_list(),
        Rcpp::Named("aquifer") = aquifer.get_calibratable_params_list());
}

RCPP_MODULE(ParametersModule)
{
    Rcpp::class_<Parameters>("Parameters")
        .constructor<Rcpp::List, Rcpp::List>()
        .method("get_all_params_list", &Parameters::get_all_params_list)
        .method("get_const_params_list", &Parameters::get_const_params_list)
        .method("get_calibratable_params_list", &Parameters::get_calibratable_params_list);
}