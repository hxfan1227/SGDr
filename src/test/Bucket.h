// #ifndef BUCKET_H
// #define BUCKET_H

// #include <Rcpp.h>

// class Bucket
// {
// private:
//     int layer; 
//     double z;
//     double rho_b;
//     double mc;
//     double rho_s;
//     double SAT;
//     double n;
//     double WPmm;
//     double WP;
//     double FC;
//     double FCmm;
//     double Swres;
//     double Ksat;
//     double Ksatmd;
//     double TTperc;
//     double phi_soil;
//     double AWClyr;
//     double Y;
//     double Z;
//     double interp1(Rcpp::NumericVector x, Rcpp::NumericVector y, double xi);
//     void calculate_bucket_params( Rcpp::List Params);

// public:
//     Bucket( Rcpp::List bucketParams,  Rcpp::List Params);
//     int get_layer();
//     double get_z();
//     double get_rho_b();
//     double get_mc();
//     double get_rho_s();
//     double get_SAT();
//     double get_n();
//     double get_WPmm();
//     double get_WP();
//     double get_FC();
//     double get_FCmm();
//     double get_Swres();
//     double get_Ksat();
//     double get_phi_soil();
//     double get_Y();
//     double get_Z();
//     double get_Ksatmd();
//     double get_TTperc();
//     double get_AWClyr();
//     Rcpp::List get_all_params_list();
//     Rcpp::List get_calibratable_params_list();
// };

// #endif // BUCKET_H
