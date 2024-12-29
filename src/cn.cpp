#include "cn.h"
#include "bucket.h"

CurveNumber::CurveNumber(const Rcpp::List cnParams, Bucket &SB1, Bucket &SB2)
{
    if (!cnParams.containsElementNamed("CN2"))
    {
        Rcpp::stop("Curve number parameter CN2 not found");
    }
    CN2_ = Rcpp::as<double>(cnParams["CN2"]);
    if (!cnParams.containsElementNamed("PIa"))
    {
        Rcpp::stop("Curve number parameter PIa not found");
    }
    pF_ = Rcpp::as<double>(cnParams["PIa"]);
    CN1_ = CN2_ - 20 * (100 - CN2_) / (100 - CN2_ + exp(2.533 - 0.0636 * (100 - CN2_)));
    CN3_ = CN2_ * exp(0.00673 * (100 - CN2_));
    Smax_ = 25.4 * ((1000 / CN1_) - 10);
    S3_ = 25.4 * ((1000 / CN3_) - 10);
    w2_ = (log(((SB1.phi_f() + SB2.phi_f()) / (1 - S3_ / Smax_)) - (SB1.phi_f() + SB2.phi_f())) - log((SB1.phi_s() + SB2.phi_s()) / (1 - 2.54 / Smax_) - (SB1.phi_s() + SB2.phi_s()))) / ((SB1.phi_s() + SB2.phi_s()) - (SB1.phi_f() + SB2.phi_f()));
    w1_ = log((SB1.phi_f() + SB2.phi_f()) / (1 - S3_ / Smax_) - (SB1.phi_f() + SB2.phi_f())) + w2_ * (SB1.phi_f() + SB2.phi_f());
}

void CurveNumber::initialize(int sim_length)
{
    Phi_ = Rcpp::NumericVector(sim_length, 0.0);
    S_ = Rcpp::NumericVector(sim_length, 0.0);
    Ia_ = Rcpp::NumericVector(sim_length, 0.0);
    CN_ = Rcpp::NumericVector(sim_length, 0.0);
    Q_ = Rcpp::NumericVector(sim_length, 0.0);
}

void CurveNumber::calc_runoff(int i, const Rcpp::NumericVector &R, Bucket &bucket1, Bucket &bucket2)
{

    Phi_[i] = bucket1.phi(i) + bucket2.phi(i) - (bucket1.phi_w() + bucket2.phi_w());
    S_[i] = Smax_ * (1 - Phi(i) / (Phi(i) + exp(w1_ - w2_ * Phi(i))));
    Ia_[i] = S(i) * 0.2;
    CN_[i] = 25400 / (S(i) + 254);
    Q_[i] = (R[i] < Ia(i)) ? 0 : pow((R[i] - Ia(i)), 2) / (R[i] + 0.8 * S(i));
}

double CurveNumber::Phi(int i) { return Phi_[i]; }
double CurveNumber::S(int i) { return S_[i]; }
double CurveNumber::Ia(int i) { return Ia_[i]; }
double CurveNumber::CN(int i) { return CN_[i]; }
double CurveNumber::Q(int i) { return Q_[i]; }
Rcpp::NumericVector CurveNumber::Phi() { return Phi_; }
Rcpp::NumericVector CurveNumber::S() { return S_; }
Rcpp::NumericVector CurveNumber::Ia() { return Ia_; }
Rcpp::NumericVector CurveNumber::CN() { return CN_; }
Rcpp::NumericVector CurveNumber::Q() { return Q_; }
double CurveNumber::CN2() { return CN2_; }
double CurveNumber::pF() { return pF_; }
double CurveNumber::CN1() { return CN1_; }
double CurveNumber::CN3() { return CN3_; }
double CurveNumber::Smax() { return Smax_; }
double CurveNumber::S3() { return S3_; }
double CurveNumber::w1() { return w1_; }
double CurveNumber::w2() { return w2_; }
void CurveNumber::reset(int warmup)
{
    Phi_[0] = Phi_[warmup - 1];
    S_[0] = S_[warmup - 1];
    Ia_[0] = Ia_[warmup - 1];
    CN_[0] = CN_[warmup - 1];
    Q_[0] = Q_[warmup - 1];
}

Rcpp::List CurveNumber::get_all_params_list()
{

    return Rcpp::List::create(Rcpp::Named("CN2") = CN2_,
                              Rcpp::Named("PIa") = pF_,
                              Rcpp::Named("CN1") = CN1_,
                              Rcpp::Named("CN3") = CN3_,
                              Rcpp::Named("Smax") = Smax_,
                              Rcpp::Named("S3") = S3_,
                              Rcpp::Named("w1") = w1_,
                              Rcpp::Named("w2") = w2_);
}

Rcpp::List CurveNumber::get_calibratable_params_list()
{
    return Rcpp::List::create(Rcpp::Named("CN2") = CN2_,
                              Rcpp::Named("PIa") = pF_);
}