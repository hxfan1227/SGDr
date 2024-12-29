#ifndef CURVENUMBER_H
#define CURVENUMBER_H

#include <Rcpp.h>

class Bucket; // forward declaration

class CurveNumber
{
private:
    double CN2_;  // curve number II [-]
    double pF_;   // percentage of Ia that is pre-runoff infiltration [%]
    double CN1_;  // curve number I [-]
    double CN3_;  // curve number III [-]
    double Smax_; // daily maximum retention parameter [mm H2O]
    double S3_;   // daily retention parameter at curve number III [-]
    double w1_;   // first shape coefficient [-]
    double w2_;   // second shape coefficient [-]

    // the following variables are used in the calculation
    Rcpp::NumericVector Phi_; // amount of water held in the soil (summing the soil water amounts in soil buckets 1 and 2) excluding the water held in the soil profile at wilting point [mm H2O]
    Rcpp::NumericVector S_;   // retention parameter [mm H2O]
    Rcpp::NumericVector Ia_;  // initial abstraction that includes surface storage, interception and infiltration prior to runoff [mm H2O]
    Rcpp::NumericVector CN_;  // curve number [-]
    Rcpp::NumericVector Q_;   // runoff [mm H2O]

public:
    // constructor
    CurveNumber(const Rcpp::List CNParams, Bucket &SB1, Bucket &SB2);
    // initialize the curve number based on the simulation length
    void initialize(int sim_length);
    // calculate the runoff for a given day
    void calc_runoff(int i, const Rcpp::NumericVector &R, Bucket &bucket1, Bucket &bucket2);
    // amount of water held in the soil (summing the soil water amounts in soil buckets 1 and 2) excluding the water held in the soil profile at wilting point [mm H2O]
    double Phi(int i);
    // retention parameter [mm H2O]
    double S(int i);
    // initial abstraction that includes surface storage, interception and infiltration prior to runoff [mm H2O]
    double Ia(int i);
    // curve number [-]
    double CN(int i);
    // runoff [mm H2O]
    double Q(int i);
    //[timeseries] amount of water held in the soil (summing the soil water amounts in soil buckets 1 and 2) excluding the water held in the soil profile at wilting point [mm H2O]
    Rcpp::NumericVector Phi();
    //[timeseries] retention parameter [mm H2O]
    Rcpp::NumericVector S();
    //[timeseries] initial abstraction that includes surface storage, interception and infiltration prior to runoff [mm H2O]
    Rcpp::NumericVector Ia();
    //[timeseries] curve number [-]
    Rcpp::NumericVector CN();
    //[timeseries] runoff [mm H2O]
    Rcpp::NumericVector Q();
    // curve number II [-]
    double CN2();
    // percentage of Ia that is pre-runoff infiltration [%]
    double pF();
    // curve number I [-]
    double CN1();
    // curve number III [-]
    double CN3();
    // daily maximum retention parameter [mm H2O]
    double Smax();
    // daily retention parameter at curve number III [-]
    double S3();
    // first shape coefficient [-]
    double w1();
    // second shape coefficient [-]
    double w2();
    void reset(int warmup);
    Rcpp::List get_all_params_list();
    Rcpp::List get_calibratable_params_list();
};

#endif // CURVENUMBER_H