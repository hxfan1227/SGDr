#ifndef SGD_H_
#define SGD_H_
#include <Rcpp.h>

Rcpp::NumericVector moving_average(Rcpp::NumericVector x, int windowSize);


// class CurveNumber
// {
// private:
//   double PIa;  // proportion of Ia that is pre-runoff infiltration (-)
//   double Smax; // maximum retention parameter (mm H20)
//   double w1;   // first shape coefficient (-)
//   double w2;   // second shape coefficient (-)

// public:
//   // Constructor that initializes the class with an Rcpp::List
//   CurveNumber(const Rcpp::List curveNumberParams)
//   {
//     if (!curveNumberParams.containsElementNamed("PIa"))
//     {
//       Rcpp::stop("Curve number parameter PIa not found");
//     }
//     PIa = Rcpp::as<double>(curveNumberParams["PIa"]);
//     if (!curveNumberParams.containsElementNamed("Smax"))
//     {
//       Rcpp::stop("Curve number parameter Smax not found");
//     }
//     Smax = Rcpp::as<double>(curveNumberParams["Smax"]);
//     if (!curveNumberParams.containsElementNamed("w1"))
//     {
//       Rcpp::stop("Curve number parameter w1 not found");
//     }
//     w1 = Rcpp::as<double>(curveNumberParams["w1"]);
//     if (!curveNumberParams.containsElementNamed("w2"))
//     {
//       Rcpp::stop("Curve number parameter w2 not found");
//     }
//     w2 = Rcpp::as<double>(curveNumberParams["w2"]);
//   }

//   // Getter methods for the class
//   double get_PIa() { return PIa; }
//   double get_Smax() { return Smax; }
//   double get_w1() { return w1; }
//   double get_w2() { return w2; }
// }; // End of CurveNumber class

// // Expose the class to R with the specific methods instantiated for expected types


// RCPP_MODULE(CurveNumberModule)
// {
//   Rcpp::class_<CurveNumber>("CurveNumber")
//       .constructor<Rcpp::List>()
//       .method("get_PIa", &CurveNumber::get_PIa)
//       .method("get_Smax", &CurveNumber::get_Smax)
//       .method("get_w1", &CurveNumber::get_w1)
//       .method("get_w2", &CurveNumber::get_w2);
// }
#endif // SGD_H_
