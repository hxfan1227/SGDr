#include "ConstParameter.h"
ConstParameter::ConstParameter(Rcpp::List constParams)
{
    if (!constParams.containsElementNamed("x"))
    {
        Rcpp::stop("Const parameter x not found");
    }
    x = Rcpp::as<double>(constParams["x"]);
    if (!constParams.containsElementNamed("W"))
    {
        Rcpp::stop("Const parameter W not found");
    }
    W = Rcpp::as<double>(constParams["W"]);
    if (!constParams.containsElementNamed("Area"))
    {
        Rcpp::stop("Const parameter Area not found");
    }
    Area = Rcpp::as<double>(constParams["Area"]);
    if (!constParams.containsElementNamed("L"))
    {
        Rcpp::stop("Const parameter L not found");
    }
    L = Rcpp::as<double>(constParams["L"]);
    if (!constParams.containsElementNamed("waterContent"))
    {
        Rcpp::stop("Const parameter waterContent not found");
    }
    waterContent = Rcpp::as<Rcpp::List>(constParams["waterContent"]);
}

double ConstParameter::get_x()
{
    return {x};
}

double ConstParameter::get_W()
{
    return {W};
}

double ConstParameter::get_Area()
{
    return {Area};
}

double ConstParameter::get_L()
{
    return {L};
}

Rcpp::List ConstParameter::get_all_params_list()
{
    return Rcpp::List::create(Rcpp::Named("x") = x,
                              Rcpp::Named("W") = W,
                              Rcpp::Named("Area") = Area,
                              Rcpp::Named("L") = L,
                              Rcpp::Named("waterContent") = waterContent);
}