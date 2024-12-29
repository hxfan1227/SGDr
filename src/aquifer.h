#ifndef AQIFER_H
#define AQIFER_H

#include "bucket.h"
#include "constpar.h"
#include <Rcpp.h>

class Aquifer
{
private:
    double Td_;      // groundwater delay time [d]
    double rho_s_;   // seawater density [kg/m3]
    double rho_f_;   // freshwater density [kg/m3]
    double Ka_;      // hydraulic conductivity of the aquifer layer [m/d]
    double z0_;      // aquifer depth [m]
    double Sy_;      // specific yield [-]
    double a_;       // dimensionless density difference between freshwater and seawater [-]
    double he_;      // equivalent freshwater head at aquifer base [m AHD]
    double XT_;      // initial seawater wedge toe location [m]
    double dxT_max_; // maximum seawater wedge toe migration rate [m/d]

    // the following variables are used in the calculation
    Rcpp::NumericVector Wnet_; // recharge entering the aquifer bucket [mm H2O]
    Rcpp::NumericVector _Wnet_; // average recharge [mm H2O] over period nw [d]
    Rcpp::NumericVector hf_; // hydraulic head at the beginning of the day [m AHD]
    Rcpp::NumericVector _hf_; // hydraulic head at the end of the day [m AHD]
    Rcpp::NumericVector q0_; // SFGD depending on the location of the saltwater wedge [m2/d]
    Rcpp::NumericVector q01_; // SFGD when the distance to the representative well (x) is less than the distance to the toe of the saltwater wedge (xT) [m2/d]
    Rcpp::NumericVector q02_; // SFGD when the distance to the representative well (x) is greater than the distance to the toe of the saltwater wedge (xT) [m2/d]
    Rcpp::NumericVector xn_; // the distance from the coast to the peak of the groundwater mound [m]
    Rcpp::NumericVector xn1_; // the distance from the coast to the peak of the groundwater mound when x < xT [m]
    Rcpp::NumericVector xn2_; // the distance from the coast to the peak of the groundwater mound when x > xT [m]
    Rcpp::NumericVector hn_; // the elevation of the groundwater mount (0m AHD) [m]
    Rcpp::NumericVector hn1_; // the elevation of the groundwater mount (0m AHD) when x < xT [m]
    Rcpp::NumericVector hn2_; // the elevation of the groundwater mount (0m AHD) when x > xT [m]
    Rcpp::NumericVector M1_; // mixed convection ratio when x < xT [-]
    Rcpp::NumericVector M2_; // mixed convection ratio when x > xT [-]
    Rcpp::NumericVector xT1_; // positioning of the saltwater toe when x < xT [m]
    Rcpp::NumericVector xT2_; // positioning of the saltwater toe when x > xT [m]
    Rcpp::NumericVector dh0_; // head drop in aquifer layer due to SFGD [m H2O]
    Rcpp::NumericVector dh1_; // head drop in aquifer layer due to pumping extraction [m H2O]
    Rcpp::NumericVector dh2_; // head drop in aquifer layer due to saltwater wedge migration [m H2O]
    Rcpp::NumericVector dxT_; // saltwater toe moving speed [m/d]
    Rcpp::NumericVector _XT_; // saltwater toe location at the end of the day[m]
    Rcpp::NumericVector V_; // volume of salt water in the aquifer bucket [m3]
    
    // private methods

    // calculate SFGD depending on the location of the saltwater wedge [m2/d]
    // void calc_sfgds(int i, ConstParameter &constpar, Rcpp::NumericVector &q1);
    // SFGD when the distance to the representative well (x) is less than the distance to the toe of the saltwater wedge (xT) [m2/d]
    double __q01(int i, ConstParameter &constpar);
    // SFGD when the distance to the representative well (x) is greater than the distance to the toe of the saltwater wedge (xT) [m2/d]
    double __q02(int i, ConstParameter &constpar);
    // calculate the moving average of the recharge entering the aquifer bucket [mm H2O] over period nw [d]
    Rcpp::NumericVector moving_average(Rcpp::NumericVector &x, int nw);
    // calculate the distance from the coast to the peak of the groundwater mound when x < xT [m]
    double __xn1(int i, ConstParameter &constpar);
    // calculate the distance from the coast to the peak of the groundwater mound when x < xT [m]
    double __xn2(int i, ConstParameter &constpar);
    // calculate the elevation of the groundwater mount (0m AHD) when x < xT [m]
    double __hn1(int i, ConstParameter &constpar);
    // calculate the elevation of the groundwater mount (0m AHD) when x < xT [m]
    double __hn2(int i, ConstParameter &constpar);
    // calculate the mixed convection ratio when x < xT [-]
    double __M1(int i, ConstParameter &constpar);
    // calculate the mixed convection ratio when x > xT [-]
    double __M2(int i, ConstParameter &constpar);
    // calculate the positioning of the saltwater toe when x < xT [m]
    double __xT1(int i, ConstParameter &constpar);
    // calculate the positioning of the saltwater toe when x > xT [m]
    double __xT2(int i, ConstParameter &constpar);

public:
    // constructor
    Aquifer(const Rcpp::List aquiferParams);
    // initialize the aquifer object
    void initialize(int sim_length, double init_hf);
    // groundwater delay time [d]
    double Td();
    // seawater density [kg/m3]
    double rho_s();
    // freshwater density [kg/m3]
    double rho_f();
    // hydraulic conductivity of the aquifer layer [m/d]
    double Ka();
    // aquifer depth [m]
    double z0();
    // specific yield [-]
    double Sy();
    // dimensionless density difference between freshwater and seawater [-]
    double get_a();
    // equivalent freshwater head at aquifer base [m AHD]
    double get_he();
    // initial seawater wedge toe location [m]
    double XT();
    // maximum seawater wedge toe migration rate [m/d]
    double dxT_max();
    // recharge entering the aquifer bucket [mm H2O]
    double Wnet(int i);
    // average recharge [mm H2O] over period nw [d]
    double _Wnet(int i);
    // hydraulic head at the beginning of the day [m AHD]
    double hf(int i);
    // hydraulic head at the end of the day [m AHD]
    double _hf(int i);
    // SFGD depending on the location of the saltwater wedge [m2/d]
    double q0(int i);
    // distance from the coast to the peak of the groundwater mound [m]
    double xn(int i);
    // the elevation of the groundwater mount (0m AHD) [m]
    double hn(int i);
    // head drop in aquifer layer due to SFGD [m H2O]
    double dh0(int i);
    // head drop in aquifer layer due to pumping extraction [m H2O]
    double dh1(int i);
    // head drop in aquifer layer due to saltwater wedge migration [m H2O]
    double dh2(int i);
    // saltwater toe moving speed [m/d]
    double dxT(int i);
    // saltwater toe location at the end of the day[m]
    double _XT(int i);
    // calculate the recharge entering the aquifer bucket [mm H2O]
    void calc_recharge(int i, Bucket &bucket);
    // average the recharge entering the aquifer bucket [mm H2O] over period nw [d]
    void average_recharge(int nw);
    // update the aquifer object
    void update(int i);
    // calculate SFGD depending on the location of the saltwater wedge [m2/d]
    void calc_sfgd(int sim_length, ConstParameter &constpar, Rcpp::NumericVector &q1);
    // reset the aquifer object
    void reset(int warmup);
    Rcpp::List get_all_params_list();
    Rcpp::List get_calibratable_params_list();
    //[timeseries] recharge entering the aquifer bucket [mm H2O] 
    Rcpp::NumericVector Wnet(); 
    //[timeseries] average recharge [mm H2O] over period nw [d]
    Rcpp::NumericVector _Wnet();
    //[timeseries] hydraulic head at the beginning of the day [m AHD]
    Rcpp::NumericVector hf();
    //[timeseries] hydraulic head at the end of the day [m AHD]
    Rcpp::NumericVector _hf(); 
    //[timeseries] SFGD depending on the location of the saltwater wedge [m2/d]
    Rcpp::NumericVector q0();
    //[timeseries] distance from the coast to the peak of the groundwater mound [m]
    Rcpp::NumericVector xn();
    //[timeseries] the elevation of the groundwater mount (0m AHD) [m]
    Rcpp::NumericVector hn();
    //[timeseries] head drop in aquifer layer due to SFGD [m H2O]
    Rcpp::NumericVector dh0();
    //[timeseries] head drop in aquifer layer due to pumping extraction [m H2O]
    Rcpp::NumericVector dh1();
    //[timeseries] head drop in aquifer layer due to saltwater wedge migration [m H2O]
    Rcpp::NumericVector dh2();
    //[timeseries] saltwater toe moving speed [m/d]
    Rcpp::NumericVector dxT();
    //[timeseries] saltwater toe location at the end of the day[m]
    Rcpp::NumericVector _XT();
    //[timeseries] volume of salt water in the aquifer bucket [m3]
    Rcpp::NumericVector V();
};
#endif // AQIFER_H