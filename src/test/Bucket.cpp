// #include "Bucket.h"

// double Bucket::interp1(Rcpp::NumericVector x, Rcpp::NumericVector y, double xi)
// {
//   Rcpp::Environment pracmaEnv = Rcpp::Environment::namespace_env("pracma");
//   Rcpp::Function interp1_r = pracmaEnv["interp1"];
//   return Rcpp::as<double>(interp1_r(x, y, xi));
// }

// Bucket::Bucket(Rcpp::List bucketParams, Rcpp::List Params)
// {
//   // ructor implementation remains the same
//   if (!bucketParams.containsElementNamed("layer"))
//   {
//     Rcpp::stop("Bucket parameter layer not found");
//   }
//   layer = Rcpp::as<int>(bucketParams["layer"]);
//   if (!bucketParams.containsElementNamed("z"))
//   {
//     Rcpp::stop("Bucket parameter z not found");
//   }
//   z = Rcpp::as<double>(bucketParams["z"]);
//   if (!bucketParams.containsElementNamed("rho_b"))
//   {
//     Rcpp::stop("Bucket parameter rho_b not found");
//   }
//   rho_b = Rcpp::as<double>(bucketParams["rho_b"]);
//   if (!bucketParams.containsElementNamed("mc"))
//   {
//     Rcpp::stop("Bucket parameter mc not found");
//   }
//   mc = Rcpp::as<double>(bucketParams["mc"]);
//   if (!bucketParams.containsElementNamed("rho_s"))
//   {
//     Rcpp::stop("Bucket parameter rho_s not found");
//   }
//   rho_s = Rcpp::as<double>(bucketParams["rho_s"]);
//   if (!bucketParams.containsElementNamed("Ksat"))
//   {
//     Rcpp::stop("Bucket parameter Ksat not found");
//   }
//   Ksat = Rcpp::as<double>(bucketParams["Ksat"]);
//   if (!bucketParams.containsElementNamed("n"))
//   {
//     Rcpp::stop("Bucket parameter n not found");
//   }
//   n = Rcpp::as<double>(bucketParams["n"]);
//   if (!bucketParams.containsElementNamed("Y"))
//   {
//     if (layer == 1)
//     {
//       Rcpp::stop("Bucket parameter Y not found");
//     }
//     else
//     {
//       Y = 0;
//     }
//   }
//   else
//   {
//     Y = Rcpp::as<double>(bucketParams["Y"]);
//   }
//   if (!bucketParams.containsElementNamed("Z"))
//   {
//     if (layer > 1)
//     {
//       Rcpp::stop("Bucket parameter Z not found");
//     }
//     else
//     {
//       Z = 0;
//     }
//   }
//   else
//   {
//     Z = Rcpp::as<double>(bucketParams["Z"]);
//   }
//   calculate_bucket_params(Params);
// }

// int Bucket::get_layer() { return layer; }
// double Bucket::get_z() { return z; }
// double Bucket::get_rho_b() { return rho_b; }
// double Bucket::get_mc() { return mc; }
// double Bucket::get_rho_s() { return rho_s; }
// double Bucket::get_SAT() { return SAT; }
// double Bucket::get_n() { return n; }
// double Bucket::get_WPmm() { return WPmm; }
// double Bucket::get_WP() { return WP; }
// double Bucket::get_FC() { return FC; }
// double Bucket::get_FCmm() { return FCmm; }
// double Bucket::get_Swres() { return Swres; }
// double Bucket::get_Ksat() { return Ksat; }
// double Bucket::get_phi_soil() { return phi_soil; }
// double Bucket::get_Y() { return Y; }
// double Bucket::get_Z() { return Z; }
// double Bucket::get_Ksatmd() { return Ksatmd; }
// double Bucket::get_TTperc() { return TTperc; }
// double Bucket::get_AWClyr() { return AWClyr; }
// Rcpp::List Bucket::get_all_params_list()
// {
//   return Rcpp::List::create(Rcpp::Named("layer") = layer,
//                             Rcpp::Named("z") = z,
//                             Rcpp::Named("rho_b") = rho_b,
//                             Rcpp::Named("mc") = mc,
//                             Rcpp::Named("rho_s") = rho_s,
//                             Rcpp::Named("SAT") = SAT,
//                             Rcpp::Named("n") = n,
//                             Rcpp::Named("WPmm") = WPmm,
//                             Rcpp::Named("WP") = WP,
//                             Rcpp::Named("FC") = FC,
//                             Rcpp::Named("FCmm") = FCmm,
//                             Rcpp::Named("Swres") = Swres,
//                             Rcpp::Named("Ksat") = Ksat,
//                             Rcpp::Named("phi_soil") = phi_soil,
//                             Rcpp::Named("Y") = Y,
//                             Rcpp::Named("Z") = Z,
//                             Rcpp::Named("Ksatmd") = Ksatmd,
//                             Rcpp::Named("TTperc") = TTperc,
//                             Rcpp::Named("AWClyr") = AWClyr);
// }

// void Bucket::calculate_bucket_params(Rcpp::List Params)
// {
//   if (!Params.containsElementNamed("waterContent"))
//   {
//     Rcpp::stop(" parameter waterContent not found");
//   }
//   Rcpp::List waterContent = Rcpp::as<Rcpp::List>(Params["waterContent"]);
//   phi_soil = 1 - rho_b / rho_s;
//   SAT = phi_soil * z;
//   WP = 0.4 * ((mc * rho_b) / 100);
//   WPmm = WP * z;
//   Swres = WPmm / SAT;
//   FC = interp1(waterContent["ClayContent"], waterContent["FC"], mc);
//   FCmm = FC * z;
//   AWClyr = FC - WP;
//   Ksatmd = Ksat / 1000 * 24;
//   TTperc = (SAT - FCmm) / Ksat;
// }

// Rcpp::List Bucket::get_calibratable_params_list()
// {
//   return Rcpp::List::create(Rcpp::Named("layer") = layer,
//                             Rcpp::Named("z") = z,
//                             Rcpp::Named("rho_b") = rho_b,
//                             Rcpp::Named("mc") = mc,
//                             Rcpp::Named("rho_s") = rho_s,
//                             Rcpp::Named("Ksat") = Ksat,
//                             Rcpp::Named("n") = n,
//                             Rcpp::Named("Y") = Y,
//                             Rcpp::Named("Z") = Z);
// }