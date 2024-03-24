// #include "CurveNumber.h"
// #include "Bucket.h"

// CurveNumber::CurveNumber(const Rcpp::List cnParams, Bucket &SB1, Bucket &SB2)
// {
//     if (!cnParams.containsElementNamed("CN2"))
//     {
//         Rcpp::stop("Curve number parameter CN2 not found");
//     }
//     CN2 = Rcpp::as<double>(cnParams["CN2"]);
//     if (!cnParams.containsElementNamed("PIa"))
//     {
//         Rcpp::stop("Curve number parameter PIa not found");
//     }
//     PIa = Rcpp::as<double>(cnParams["PIa"]);
//     CN1 = CN2 - 20 * (100 - CN2) / (100 - CN2 + exp(2.533 - 0.0636 * (100 - CN2)));
//     CN3 = CN2 * exp(0.00673 * (100 - CN2));
//     Smax = 25.4 * ((1000 / CN1) - 10);
//     S3 = 25.4 * ((1000 / CN3) - 10);
//     w2 = (log(((SB1.get_FCmm() + SB2.get_FCmm()) / (1 - S3 / Smax)) - (SB1.get_FCmm() + SB2.get_FCmm())) - log((SB1.get_SAT() + SB2.get_SAT()) / (1 - 2.54 / Smax) - (SB1.get_SAT() + SB2.get_SAT()))) / ((SB1.get_SAT() + SB2.get_SAT()) - (SB1.get_FCmm() + SB2.get_FCmm()));
//     w1 = log((SB1.get_FCmm() + SB2.get_FCmm()) / (1 - S3 / Smax) - (SB1.get_FCmm() + SB2.get_FCmm())) + w2 * (SB1.get_FCmm() + SB2.get_FCmm());
// }

// double CurveNumber::get_CN2()
// {
//     return CN2;
// }

// double CurveNumber::get_PIa()
// {
//     return PIa;
// }

// double CurveNumber::get_CN1()
// {
//     return CN1;
// }

// double CurveNumber::get_CN3()
// {
//     return CN3;
// }

// double CurveNumber::get_Smax()
// {
//     return Smax;
// }

// double CurveNumber::get_S3()
// {
//     return S3;
// }

// double CurveNumber::get_w1()
// {
//     return w1;
// }

// double CurveNumber::get_w2()
// {
//     return w2;
// }

// Rcpp::List CurveNumber::get_all_params_list()
// {

//     return Rcpp::List::create(Rcpp::Named("CN2") = CN2,
//                               Rcpp::Named("PIa") = PIa,
//                               Rcpp::Named("CN1") = CN1,
//                               Rcpp::Named("CN3") = CN3,
//                               Rcpp::Named("Smax") = Smax,
//                               Rcpp::Named("S3") = S3,
//                               Rcpp::Named("w1") = w1,
//                               Rcpp::Named("w2") = w2);
// }

// Rcpp::List CurveNumber::get_calibratable_params_list()
// {
//     return Rcpp::List::create(Rcpp::Named("CN2") = CN2,
//                               Rcpp::Named("PIa") = PIa);
// }
