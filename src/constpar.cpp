#include "constpar.h"

ConstParameter::ConstParameter(Rcpp::List pars)
{
    if (!pars.containsElementNamed("x"))
    {
        Rcpp::stop("Const parameter x not found");
    }
    x_ = Rcpp::as<double>(pars["x"]);
    if (!pars.containsElementNamed("W"))
    {
        Rcpp::stop("Const parameter W not found");
    }
    W_ = Rcpp::as<double>(pars["W"]);
    if (!pars.containsElementNamed("A"))
    {
        Rcpp::stop("Const parameter A not found");
    }
    A_ = Rcpp::as<double>(pars["A"]);
    if (!pars.containsElementNamed("L"))
    {
        Rcpp::stop("Const parameter L not found");
    }
    L_ = Rcpp::as<double>(pars["L"]);
    if (!pars.containsElementNamed("WC"))
    {
        Rcpp::stop("Const parameter WC not found");
    }
    WC_ = Rcpp::as<Rcpp::List>(pars["WC"]);
}

double ConstParameter::x() { return x_; }
double ConstParameter::W() { return W_; }
double ConstParameter::A() { return A_; }
double ConstParameter::L() { return L_; }
Rcpp::List ConstParameter::const_pars()
{
    return Rcpp::List::create(Rcpp::Named("x") = x_,
                              Rcpp::Named("W") = W_,
                              Rcpp::Named("Area") = A_,
                              Rcpp::Named("L") = L_,
                              Rcpp::Named("WC") = WC_);
}