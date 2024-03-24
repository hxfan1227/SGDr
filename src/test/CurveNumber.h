// #ifndef CURVENUMBER_H
// #define CURVENUMBER_H

// #include <Rcpp.h>
// #include "Bucket.h"

// class CurveNumber
// {
// private:
//     double CN2;  // Curve number 2
//     double PIa;  // Proportion of Ia that is pre-runoff infiltration         
//     double CN1;  // Curve number 1
//     double CN3;  // Curve number 3
//     double Smax; // Maximum retention
//     double S3;   // Retention at CN3
//     double w1;   // First shape coefficient
//     double w2;   // Second shape coefficient

// public:
//     CurveNumber(const Rcpp::List CNParams, Bucket& SB1, Bucket& SB2);
//     double get_CN2();
//     double get_PIa();
//     double get_CN1();
//     double get_CN3();
//     double get_Smax();
//     double get_S3();
//     double get_w1();
//     double get_w2();
//     Rcpp::List get_all_params_list();
//     Rcpp::List get_calibratable_params_list();
// };

// #endif // CURVENUMBER_H