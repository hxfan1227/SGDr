#ifndef MODEL_H
#define MODEL_H
#include "Parameters.h"
class Model
{
private:
    Rcpp::DataFrame inputData;
    Parameters parameters;
    int simLength;
    int warmUp;
    Bucket &bucket1 = parameters.get_bucket1();
    Bucket &bucket2 = parameters.get_bucket2();
    CurveNumber &cn = parameters.get_curveNumber();
    Aquifer &aquifer = parameters.get_aquifer();
    ConstParameter &constpar = parameters.get_constParameter();
    Rcpp::NumericVector R;       // precipitation (mm)
    Rcpp::NumericVector E0;      // potential evapotranspiration (mm)
    Rcpp::NumericVector q1; // pumping rate (m3/day)
    Rcpp::NumericVector hf;     // groundwater level (m)
    Rcpp::NumericVector phi1; // head in soil bucket 1 (mm)
    Rcpp::NumericVector phi2; // head in soil bucket 2 (mm)

public:
    // constructor
    Model(Rcpp::DataFrame inputData,
          const Rcpp::List &calibratableParams,
          const Rcpp::List &constParams, int warmUp = 0);
    void calc(int nw, int period);
    void run(int windowSize, int warmUp);
    Rcpp::List estimate_sgd(const Rcpp::DataFrame &inputData, const Rcpp::List &calibratableParams, const Rcpp::List &constParams, int windowSize, int warmUp);
    Rcpp::DataFrame get_recharge_output();
    Rcpp::DataFrame get_sgd_output();
    Rcpp::List get_all_params_list();
    Rcpp::DataFrame get_inputData();
    Rcpp::NumericVector get_Qsurf();

};

#endif // MODEL_H