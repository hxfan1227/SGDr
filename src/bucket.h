#ifndef BUCKET_H
#define BUCKET_H

#include <Rcpp.h>

class CurveNumber; // forward declaration

class Bucket
{
private:
    int layer_;      // layer number
    double z_;       // soil depth [mm]
    double rho_b_;   // bulk density [Mg/m3]
    double m_;       // clay percentage [%]
    double rho_p_;   // particle density [Mg/m3]
    double phi_s_;   // amount of water in the soil bucket when completely saturated [mm H2O]
    double n_;       // empirical moisture retention parameter [-]
    double phi_w_;   // amount of water in the soil bucket at wilting point [mm H2O]
    double theta_w_; // wilting point of the soil bucket [-]
    double theta_f_; // field capacity of the soil bucket [-]
    double phi_f_;   // amount of water in the soil bucket at field capacity [mm H2O]
    double theta_r_; // residual saturation of the soil bucket [-]
    double Ks_;      // saturated hydraulic conductivity of the soil bucket [mm/hr]
    // double Ksatmd;
    // double TTperc;
    double theta_s_; // saturation of the soil bucket [-]
    // double AWClyr;
    double pE_; // percentage of evaporation from soil bucket 1 in the total evaporation from the soil buckets [%]
    double pZ_; // percentage of capacity between field capacity and saturation that is available for percolation [%]
    double interp1(Rcpp::NumericVector x, Rcpp::NumericVector y, double xi);
    void calculate_bucket_params(Rcpp::List Params);

    // the following variables are used in the calculation

    Rcpp::NumericVector phi_;   // amount of water held in the soil excluding the water held in the soil profile at wilting point at the beginning of the day [mm H2O]
    Rcpp::NumericVector _phi_;  // amount of water held in the soil excluding the water held in the soil profile at wilting point at the end of the day [mm H2O]
    Rcpp::NumericVector F_;     // infiltration without considering the pre-runoff infiltration [mm H2O]
    Rcpp::NumericVector FI_;    // pre-runoff infiltration [mm H2O]
    Rcpp::NumericVector I_;     // amount of water intercepted by the vegetation canopy [mm H2O]
    Rcpp::NumericVector Ec_;    // canopy evaporation [mm H2O]
    Rcpp::NumericVector Es_;    // soil evaporation [mm H2O]
    Rcpp::NumericVector Esi_;   // soil evaporation in the current bucket [mm H2O]
    Rcpp::NumericVector Ea_;    // actual evaporation from soil bucket i [mm H2O]
    Rcpp::NumericVector phi_e_; // drainable volume of water in soil bucket i [mm H2O]
    Rcpp::NumericVector theta_; // dimensionless saturation of soil bucket i [-]
    Rcpp::NumericVector K_;     // unsaturated hydraulic conductivity of soil bucket i [mm/hr]
    Rcpp::NumericVector Tp_;    // percolation travel time in soil bucket i [hr]
    Rcpp::NumericVector P0_;    // potential percolation from soil bucket i [mm H2O]
    Rcpp::NumericVector P_;     // actual percolation from soil bucket i to the underlying bucket limited to field capacity and available space of the underlying bucket [mm H2O]

public:
    // constructor
    Bucket(Rcpp::List bucket_pars, Rcpp::List const_pars);
    // initialize the bucket based on the simulation length and initial soil moisture content
    void initialize(int sim_length, double init_phi);
    // calculate the amount of water held in the soil profile at the end of the day
    void calc_phi(int i, Rcpp::NumericVector &R, Rcpp::NumericVector &E0, CurveNumber &cn, Bucket &bucket);
    // update the amount of water held in the soil excluding the water held in the soil profile at wilting point at the beginning of the day
    void update(int i);
    // amount of water held in the soil bucket at the beginning of the day [mm H2O]
    double phi(int i);
    // amount of water held in the soil bucket at the end of the day [mm H2O]
    double _phi(int i);
    // infiltration without considering the pre-runoff infiltration [mm H2O]
    double F(int i);
    // pre-runoff infiltration [mm H2O]
    double FI(int i);
    // amount of water intercepted by the vegetation canopy [mm H2O]
    double I(int i);
    // canopy evaporation [mm H2O]
    double Ec(int i);
    // soil evaporation [mm H2O]
    double Es(int i);
    // soil evaporation in the current bucket [mm H2O]
    double Esi(int i);
    // actual evaporation from soil bucket i [mm H2O]
    double Ea(int i);
    // drainable volume of water in soil bucket i [mm H2O]
    double phi_e(int i);
    // dimensionless saturation of soil bucket i [-]
    double theta(int i);
    // unsaturated hydraulic conductivity of soil bucket i [mm/hr]
    double K(int i);
    // percolation travel time in soil bucket i [hr]
    double Tp(int i);
    // potential percolation from soil bucket i [mm H2O]
    double P0(int i);
    // actual percolation from soil bucket i to the underlying bucket limited to field capacity and available space of the underlying bucket [mm H2O]
    double P(int i);

    //[timeseries] amount of water held in the soil profile at the beginning of the day [mm H2O]
    Rcpp::NumericVector phi();
    //[timeseries] amount of water held in the soil profile at the end of the day [mm H2O]
    Rcpp::NumericVector _phi();
    //[timeseries] infiltration without considering the pre-runoff infiltration [mm H2O]
    Rcpp::NumericVector F();
    //[timeseries] pre-runoff infiltration [mm H2O]
    Rcpp::NumericVector FI();
    //[timeseries] amount of water intercepted by the vegetation canopy [mm H2O]
    Rcpp::NumericVector I();
    //[timeseries] canopy evaporation [mm H2O]
    Rcpp::NumericVector Ec();
    //[timeseries] soil evaporation [mm H2O]
    Rcpp::NumericVector Es();
    //[timeseries] soil evaporation in the current bucket [mm H2O]
    Rcpp::NumericVector Esi();
    //[timeseries] actual evaporation from soil bucket i [mm H2O]
    Rcpp::NumericVector Ea();
    //[timeseries] drainable volume of water in soil bucket i [mm H2O]
    Rcpp::NumericVector phi_e();
    //[timeseries] dimensionless saturation of soil bucket i [-]
    Rcpp::NumericVector theta();
    //[timeseries] unsaturated hydraulic conductivity of soil bucket i [mm/hr]
    Rcpp::NumericVector K();
    //[timeseries] percolation travel time in soil bucket i [hr]
    Rcpp::NumericVector Tp();
    //[timeseries] potential percolation from soil bucket i [mm H2O]
    Rcpp::NumericVector P0();
    //[timeseries] actual percolation from soil bucket i to the underlying bucket limited to field capacity and available space of the underlying bucket [mm H2O]
    Rcpp::NumericVector P();    

    // layer number
    int layer();
    // soil depth [mm]
    double z();
    // bulk density [Mg/m3]
    double rho_b();
    // clay percentage [%]
    double m();
    // particle density [Mg/m3]
    double rho_p();
    // amount of water in the soil bucket when completely saturated [mm H2O]
    double phi_s();
    // empirical moisture retention parameter [-]
    double n();
    // amount of water in the soil bucket at wilting point [mm H2O]
    double phi_w();
    // wilting point of the soil bucket [-]
    double theta_w();
    // field capacity of the soil bucket [-]
    double theta_f();
    // amount of water in the soil bucket at field capacity [mm H2O]
    double phi_f();
    // residual saturation of the soil bucket [-]
    double theta_r();
    // saturated hydraulic conductivity of the soil bucket [mm/hr]
    double Ks();
    // saturation of the soil bucket [-]
    double theta_s();
    // percentage of evaporation from soil bucket 1 in the total evaporation from the soil buckets [%]
    double pE();
    // percentage of capacity between field capacity and saturation that is available for percolation [%]
    double pZ();
    void reset(int warmup);
    Rcpp::List all_pars();
    Rcpp::List cali_pars();
};

#endif // BUCKET_H