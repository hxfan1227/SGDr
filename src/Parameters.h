#ifndef Parameters_H
#define Parameters_H
#include <Rcpp.h>

#ifndef BUCKET_H
#define BUCKET_H

class Bucket
{
private:
    int layer;
    double z;
    double rho_b;
    double mc;
    double rho_s;
    double SAT;
    double n;
    double WPmm;
    double WP;
    double FC;
    double FCmm;
    double Swres;
    double Ksat;
    double Ksatmd;
    double TTperc;
    double phi_soil;
    double AWClyr;
    double Y;
    double Z;
    double interp1(Rcpp::NumericVector x, Rcpp::NumericVector y, double xi);
    void calculate_bucket_params(Rcpp::List Params);

public:
    Bucket(Rcpp::List bucketParams, Rcpp::List Params);
    int get_layer();
    double get_z();
    double get_rho_b();
    double get_mc();
    double get_rho_s();
    double get_SAT();
    double get_n();
    double get_WPmm();
    double get_WP();
    double get_FC();
    double get_FCmm();
    double get_Swres();
    double get_Ksat();
    double get_phi_soil();
    double get_Y();
    double get_Z();
    double get_Ksatmd();
    double get_TTperc();
    double get_AWClyr();
    Rcpp::List get_all_params_list();
    Rcpp::List get_calibratable_params_list();
};

#endif // BUCKET_H

#ifndef AQIFER_H
#define AQIFER_H

class Aquifer
{
private:
    double delta;
    double rho_s;
    double rho_f;
    double K;
    double z0;
    double Sy;
    double a;
    double he;
    double xT;
    double dxT;

public:
    Aquifer(const Rcpp::List aquiferParams);
    double get_delta();
    double get_rho_s();
    double get_rho_f();
    double get_K();
    double get_z0();
    double get_Sy();
    double get_a();
    double get_he();
    double get_xT();
    double get_dxT();
    Rcpp::List get_all_params_list();
    Rcpp::List get_calibratable_params_list();
};
#endif // AQIFER_H

#ifndef CURVENUMBER_H
#define CURVENUMBER_H

class CurveNumber
{
private:
    double CN2;  // Curve number 2
    double PIa;  // Proportion of Ia that is pre-runoff infiltration
    double CN1;  // Curve number 1
    double CN3;  // Curve number 3
    double Smax; // Maximum retention
    double S3;   // Retention at CN3
    double w1;   // First shape coefficient
    double w2;   // Second shape coefficient

public:
    CurveNumber(const Rcpp::List CNParams, Bucket &SB1, Bucket &SB2);
    double get_CN2();
    double get_PIa();
    double get_CN1();
    double get_CN3();
    double get_Smax();
    double get_S3();
    double get_w1();
    double get_w2();
    Rcpp::List get_all_params_list();
    Rcpp::List get_calibratable_params_list();
};

#endif // CURVENUMBER_H

#ifndef ConstParameter_H
#define ConstParameter_H

class ConstParameter
{
private:
    double x;
    double W;
    double Area;
    double L;
    Rcpp::List waterContent;

public:
    ConstParameter(Rcpp::List Params);
    double get_x();
    double get_W();
    double get_Area();
    double get_L();
    Rcpp::List get_all_params_list();
};
#endif // ConstParameter_H

class Parameters
{
private:
    Bucket bucket1;
    Bucket bucket2;
    CurveNumber curveNumber;
    Aquifer aquifer;
    ConstParameter constParameter;

public:
    Parameters(Rcpp::List calibratableParams, Rcpp::List constParams);
    Bucket get_bucket1();
    Bucket get_bucket2();
    CurveNumber get_curveNumber();
    Aquifer get_aquifer();
    ConstParameter get_constParameter();
    Rcpp::List get_calibratable_params_list();
    Rcpp::List get_const_params_list();
    Rcpp::List get_all_params_list();
    void update(const Rcpp::List& newCalibratableParams, const Rcpp::List& newConstParams);
};
#endif // Parameters_H