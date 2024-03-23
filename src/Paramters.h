#ifndef Parameters_H
#define Parameters_H
#include <Rcpp.h>
#include "Bucket.h"
#include "Aquifer.h"
#include "ConstParameter.h"
#include "CurveNumber.h"

class Parameters
{
private:
    Bucket bucket1;
    Bucket bucket2;    
    CurveNumber curveNumber;
    Aquifer aquifer;
    ConstParameter constParameter;
public:
    Parameters(Rcpp::List calibratableParams, Rcpp::List constParams);
    Bucket get_bucket1();
    Bucket get_bucket2();
    CurveNumber get_curveNumber();
    Aquifer get_aquifer();
    ConstParameter get_constParameter();
    Rcpp::List get_calibratable_params_list();
    Rcpp::List get_const_params_list();
    Rcpp::List get_all_params_list();
}; 

#endif // Parameters_H