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
    Bucket &bucket1 = parameters.bucket1();
    Bucket &bucket2 = parameters.bucket2();
    CurveNumber &cn = parameters.cn();
    Aquifer &aquifer = parameters.aquifer();
    ConstParameter &constpar = parameters.constpar();
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
    Rcpp::DataFrame get_recharge_output();
    Rcpp::DataFrame q0();
    Rcpp::List pars();
    Rcpp::DataFrame input();
    Rcpp::NumericVector Q();

};

Rcpp::List estimate_sgd(const Rcpp::DataFrame& inputData, const Rcpp::List& calibratableParams, const Rcpp::List& constParams, int nw, int warmup);

#endif // MODEL_H