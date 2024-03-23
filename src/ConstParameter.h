#ifndef PARAMETER_H
#define PARAMETER_H
#include <Rcpp.h>

class ConstParameter
{
private:
    double x;
    double W;
    double Area;
    double L;
    Rcpp::List waterContent;
public:
    ConstParameter( Rcpp::List Params);
    double get_x();
    double get_W();
    double get_Area();
    double get_L();
    Rcpp::List get_all_params_list();

};
#endif // PARAMETER_H