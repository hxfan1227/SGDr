#ifndef ConstParameter_H
#define ConstParameter_H

#include <Rcpp.h>

class ConstParameter
{
private:
    double x_; // distance to representative well [m]
    double W_; // shoreline width of region [m]
    double A_; // area of the region [m2]
    double L_; // length to the inland boundary
    Rcpp::List WC_; // water content table based on the SWAT manual.

public:
    ConstParameter(Rcpp::List pars);
    double x();
    double W();
    double A();
    double L();
    Rcpp::List const_pars();
};
#endif // ConstParameter_H