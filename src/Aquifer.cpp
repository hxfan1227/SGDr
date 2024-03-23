#include "Aquifer.h"

Aquifer::Aquifer(const Rcpp::List aquiferParams)
{
    if (!aquiferParams.containsElementNamed("delta"))
    {
        Rcpp::stop("Aquifer parameter delta not found");
    }
    delta = Rcpp::as<double>(aquiferParams["delta"]);
    if (!aquiferParams.containsElementNamed("rho_s"))
    {
        Rcpp::stop("Aquifer parameter rho_s not found");
    }
    rho_s = Rcpp::as<double>(aquiferParams["rho_s"]);
    if (!aquiferParams.containsElementNamed("rho_f"))
    {
        Rcpp::stop("Aquifer parameter rho_f not found");
    }
    rho_f = Rcpp::as<double>(aquiferParams["rho_f"]);
    if (!aquiferParams.containsElementNamed("K"))
    {
        Rcpp::stop("Aquifer parameter K not found");
    }
    K = Rcpp::as<double>(aquiferParams["K"]);
    if (!aquiferParams.containsElementNamed("z0"))
    {
        Rcpp::stop("Aquifer parameter z0 not found");
    }
    z0 = Rcpp::as<double>(aquiferParams["z0"]);
    if (!aquiferParams.containsElementNamed("Sy"))
    {
        Rcpp::stop("Aquifer parameter Sy not found");
    }
    Sy = Rcpp::as<double>(aquiferParams["Sy"]);
    if (!aquiferParams.containsElementNamed("xT"))
    {
        Rcpp::stop("Aquifer parameter xT not found");
    }
    xT = Rcpp::as<double>(aquiferParams["xT"]);
    if (!aquiferParams.containsElementNamed("dxT"))
    {
        Rcpp::stop("Aquifer parameter dxT not found");
    }
    dxT = Rcpp::as<double>(aquiferParams["dxT"]);
    a = (rho_s - rho_f) / rho_f;
    he = z0 * rho_s / rho_f - z0;
}

double Aquifer::get_delta()
{
    return {delta};
}

double Aquifer::get_rho_s()
{
    return {rho_s};
}

double Aquifer::get_rho_f()
{
    return {rho_f};
}

double Aquifer::get_K()
{
    return {K};
}

double Aquifer::get_z0()
{
    return {z0};
}

double Aquifer::get_Sy()
{
    return {Sy};
}

double Aquifer::get_a()
{
    return {a};
}

double Aquifer::get_he()
{
    return {he};
}

double Aquifer::get_xT()
{
    return {xT};
}

double Aquifer::get_dxT()
{
    return {dxT};
}

Rcpp::List Aquifer::get_all_params_list()
{
    return Rcpp::List::create(Rcpp::Named("delta") = delta,
                              Rcpp::Named("rho_s") = rho_s,
                              Rcpp::Named("rho_f") = rho_f,
                              Rcpp::Named("K") = K,
                              Rcpp::Named("z0") = z0,
                              Rcpp::Named("Sy") = Sy,
                              Rcpp::Named("a") = a,
                              Rcpp::Named("he") = he,
                              Rcpp::Named("xT") = xT,
                              Rcpp::Named("dxT") = dxT);
}

Rcpp::List Aquifer::get_calibratable_params_list()
{
    return Rcpp::List::create(Rcpp::Named("delta") = delta,
                              Rcpp::Named("K") = K,
                              Rcpp::Named("z0") = z0,
                              Rcpp::Named("Sy") = Sy,
                              Rcpp::Named("xT") = xT,
                              Rcpp::Named("dxT") = dxT);
}

