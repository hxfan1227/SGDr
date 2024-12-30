#include "bucket.h"
#include "cn.h"

double Bucket::interp1(Rcpp::NumericVector x, Rcpp::NumericVector y, double xi)
{
    Rcpp::Environment pracmaEnv = Rcpp::Environment::namespace_env("pracma");
    Rcpp::Function interp1_r = pracmaEnv["interp1"];
    return Rcpp::as<double>(interp1_r(x, y, xi));
}

Bucket::Bucket(Rcpp::List bucket_pars, Rcpp::List const_pars)
{
    if (!bucket_pars.containsElementNamed("layer"))
    {
        Rcpp::stop("Bucket parameter layer not found");
    }
    layer_ = Rcpp::as<int>(bucket_pars["layer"]);

    if (!bucket_pars.containsElementNamed("z"))
    {
        Rcpp::stop("Bucket parameter z not found");
    }
    z_ = Rcpp::as<double>(bucket_pars["z"]);

    if (!bucket_pars.containsElementNamed("rho_b"))
    {
        Rcpp::stop("Bucket parameter rho_b not found");
    }
    rho_b_ = Rcpp::as<double>(bucket_pars["rho_b"]);

    if (!bucket_pars.containsElementNamed("m"))
    {
        Rcpp::stop("Bucket parameter m not found");
    }
    m_ = Rcpp::as<double>(bucket_pars["m"]);
    
    if (!bucket_pars.containsElementNamed("rho_p"))
    {
        Rcpp::stop("Bucket parameter rho_p not found");
    }
    rho_p_ = Rcpp::as<double>(bucket_pars["rho_p"]);

    if (!bucket_pars.containsElementNamed("Ks"))
    {
        Rcpp::stop("Bucket parameter Ks not found");
    }
    Ks_ = Rcpp::as<double>(bucket_pars["Ks"]);

    if (!bucket_pars.containsElementNamed("n"))
    {
        Rcpp::stop("Bucket parameter n not found");
    }
    n_ = Rcpp::as<double>(bucket_pars["n"]);

    if (!bucket_pars.containsElementNamed("pE"))
    {
        if (layer_ == 1)
        {
            Rcpp::stop("Bucket parameter pE not found");
        }
        else
        {
            pE_ = 0;
        }
    }
    else
    {
        pE_ = Rcpp::as<double>(bucket_pars["pE"]);
    }
    if (!bucket_pars.containsElementNamed("pZ"))
    {
        if (layer_ > 1)
        {
            Rcpp::stop("Bucket parameter pZ not found");
        }
        else
        {
            pZ_ = 0;
        }
    }
    else
    {
        pZ_ = Rcpp::as<double>(bucket_pars["pZ"]);
    }
    calculate_bucket_params(const_pars);
}

void Bucket::initialize(int sim_length, double init_phi)
{
    phi_ = Rcpp::NumericVector(sim_length, init_phi);
    _phi_ = Rcpp::NumericVector(sim_length, 0.0);
    F_ = Rcpp::NumericVector(sim_length, 0.0);
    FI_ = Rcpp::NumericVector(sim_length, 0.0);
    I_ = Rcpp::NumericVector(sim_length, 0.0);
    Ec_ = Rcpp::NumericVector(sim_length, 0.0);
    Es_ = Rcpp::NumericVector(sim_length, 0.0);
    Esi_ = Rcpp::NumericVector(sim_length, 0.0);
    Ea_ = Rcpp::NumericVector(sim_length, 0.0);
    phi_e_ = Rcpp::NumericVector(sim_length, 0.0);
    theta_ = Rcpp::NumericVector(sim_length, 0.0);
    K_ = Rcpp::NumericVector(sim_length, 0.0);
    Tp_ = Rcpp::NumericVector(sim_length, 0.0);
    P0_ = Rcpp::NumericVector(sim_length, 0.0);
    P_ = Rcpp::NumericVector(sim_length, 0.0);
}

void Bucket::update(int i)
{
    if (i > 0)
    {
        phi_[i] = _phi_[i - 1];
    }
}

void Bucket::calc_phi(int i, Rcpp::NumericVector &R, Rcpp::NumericVector &E0, CurveNumber &cn, Bucket &bucket)
{
    if (layer_ == 1)
    {
        F_[i] = (R[i] < cn.Ia(i)) ? 0 : R[i] - cn.Q(i) - cn.Ia(i);
        FI_[i] = R[i] > cn.Ia(i) ? cn.pF() * cn.Ia(i) : (R[i] > (cn.Ia(i) * (1 - cn.pF())) ? R[i] - (1 - cn.pF()) * cn.Ia(i) : 0);
        _phi_[i] = phi_[i] + (F_[i] + FI_[i]);
        I_[i] = R[i] > ((1 - cn.pF()) * cn.Ia(i)) ? (1 - cn.pF()) * cn.Ia(i) : R[i];
        Ec_[i] = i == 0 ? 0.0 : std::max(Ec_[i - 1] + I_[i] - E0[i], 0.0);
        Es_[i] = i == 0 ? std::max(0.0, E0[i] - I_[i] - Ec_[i]): std::max(0.0, E0[i] - I_[i] - (Ec_[i - 1] - Ec_[i]));
        Esi_[i] = Es_[i] * pE_ / 100;
        Ea_[i] = std::min(_phi_[i] < phi_f_ ? Esi_[i] * exp(2.5 * (_phi_[i] - phi_f_) / (phi_f_ - phi_w_)) : Esi_[i], 0.8 * (_phi_[i] - phi_w_));
        _phi_[i] -= Ea_[i];
        phi_e_[i] = std::max(0.0, _phi_[i] - phi_f_);
        theta_[i] = std::min(1.0, (_phi_[i] / phi_s_ - theta_r_) / (1 - theta_r_));
        K_[i] = std::max(0.0, pow(theta_[i], 0.5) * pow((1 - pow(1 - pow(theta_[i], n_ / (n_ - 1)), (n_ - 1) / n_)), 2) * Ks_);
        Tp_[i] = std::max((phi_s_ - phi_f_) / K_[i], 72.0);
        P0_[i] = phi_e_[i] * (1 - exp(-24 / Tp_[i]));
        P_[i] = bucket.phi(i) >= bucket.phi_f() + bucket.pZ() / 100 * (bucket.phi_s() - bucket.phi_f()) ? 0.0 : (P0_[i] <= (bucket.phi_s() - bucket.phi(i)) ? P0_[i] : (bucket.phi_s() - bucket.phi(i)));
        _phi_[i] -= P_[i];
    }
    if (layer_ > 1)
    {
        _phi_[i] = phi_[i] + bucket.P(i);
        Esi_[i] = bucket.Es(i) * (100 - bucket.pE()) / 100;
        Ea_[i] = std::min((_phi_[i] < phi_f_ ? Esi_[i] * exp(2.5 * (_phi_[i] - phi_f_) / (phi_f_ - phi_w_)) : Esi_[i]),
                          0.8 * (_phi_[i] - phi_w_));
        _phi_[i] -= Ea_[i];
        phi_e_[i] = std::max(0.0, _phi_[i] - phi_f_);
        theta_[i] = std::min(1.0, (_phi_[i] / phi_s_ - theta_r_) / (1 - theta_r_));
        K_[i] = std::max(0.0, pow(theta_[i], 0.5) * pow((1 - pow(1 - pow(theta_[i], n_ / (n_ - 1)), (n_ - 1) / n_)), 2) * Ks_);
        Tp_[i] = (phi_s_ - phi_f_) / K_[i];
        P0_[i] = phi_e_[i] * (1 - exp(-24 / Tp_[i]));
        P_[i] = P0_[i];
        _phi_[i] -= P_[i];
    }
}

double Bucket::phi(int i) { return phi_[i]; }
double Bucket::_phi(int i) { return _phi_[i]; }
double Bucket::F(int i) { return F_[i]; }
double Bucket::FI(int i) { return FI_[i]; }
double Bucket::I(int i) { return I_[i]; }
double Bucket::Ec(int i) { return Ec_[i]; }
double Bucket::Es(int i) { return Es_[i]; }
double Bucket::Esi(int i) { return Esi_[i]; }
double Bucket::Ea(int i) { return Ea_[i]; }
double Bucket::phi_e(int i) { return phi_e_[i]; }
double Bucket::theta(int i) { return theta_[i]; }
double Bucket::K(int i) { return K_[i]; }
double Bucket::Tp(int i) { return Tp_[i]; }
double Bucket::P0(int i) { return P0_[i]; }
double Bucket::P(int i) { return P_[i]; }

Rcpp::NumericVector Bucket::phi() { return phi_; }
Rcpp::NumericVector Bucket::_phi() { return _phi_; }
Rcpp::NumericVector Bucket::F() { return F_; }
Rcpp::NumericVector Bucket::FI() { return FI_; }
Rcpp::NumericVector Bucket::I() { return I_; }
Rcpp::NumericVector Bucket::Ec() { return Ec_; }
Rcpp::NumericVector Bucket::Es() { return Es_; }
Rcpp::NumericVector Bucket::Esi() { return Esi_; }
Rcpp::NumericVector Bucket::Ea() { return Ea_; }
Rcpp::NumericVector Bucket::phi_e() { return phi_e_; }
Rcpp::NumericVector Bucket::theta() { return theta_; }
Rcpp::NumericVector Bucket::K() { return K_; }
Rcpp::NumericVector Bucket::Tp() { return Tp_; }
Rcpp::NumericVector Bucket::P0() { return P0_; }
Rcpp::NumericVector Bucket::P() { return P_; }

int Bucket::layer() { return layer_; }
double Bucket::z() { return z_; }
double Bucket::rho_b() { return rho_b_; }
double Bucket::m() { return m_; }
double Bucket::rho_p() { return rho_p_; }
double Bucket::phi_s() { return phi_s_; }
double Bucket::n() { return n_; }
double Bucket::phi_w() { return phi_w_; }
double Bucket::theta_w() { return theta_w_; }
double Bucket::theta_f() { return theta_f_; }
double Bucket::phi_f() { return phi_f_; }
double Bucket::theta_r() { return theta_r_; }
double Bucket::Ks() { return Ks_; }
double Bucket::theta_s() { return theta_s_; }
double Bucket::pE() { return pE_; }
double Bucket::pZ() { return pZ_; }

void Bucket::reset(int warmup)
{
    phi_[0] = phi_[warmup - 1];
    _phi_[0] = _phi_[warmup - 1];
    F_[0] = F_[warmup - 1];
    FI_[0] = FI_[warmup - 1];
    I_[0] = I_[warmup - 1];
    Ec_[0] = Ec_[warmup - 1];
    Es_[0] = Es_[warmup - 1];
    Esi_[0] = Esi_[warmup - 1];
    Ea_[0] = Ea_[warmup - 1];
    phi_e_[0] = phi_e_[warmup - 1];
    theta_[0] = theta_[warmup - 1];
    K_[0] = K_[warmup - 1];
    Tp_[0] = Tp_[warmup - 1];
    P0_[0] = P0_[warmup - 1];
    P_[0] = P_[warmup - 1];
}

Rcpp::List Bucket::all_pars()
{
    return Rcpp::List::create(Rcpp::Named("layer") = layer_,
                              Rcpp::Named("z") = z_,
                              Rcpp::Named("rho_b") = rho_b_,
                              Rcpp::Named("m") = m_,
                              Rcpp::Named("rho_p") = rho_p_,
                              Rcpp::Named("phi_s") = phi_s_,
                              Rcpp::Named("n") = n_,
                              Rcpp::Named("phi_w") = phi_w_,
                              Rcpp::Named("theta_w") = theta_w_,
                              Rcpp::Named("theta_f") = theta_f_,
                              Rcpp::Named("phi_f") = phi_f_,
                              Rcpp::Named("theta_r") = theta_r_,
                              Rcpp::Named("Ks") = Ks_,
                              Rcpp::Named("theta_s") = theta_s_,
                              Rcpp::Named("pE") = pE_,
                              Rcpp::Named("pZ") = pZ_
    );
}

void Bucket::calculate_bucket_params(Rcpp::List pars)
{
    if (!pars.containsElementNamed("WC"))
    {
        Rcpp::stop(" parameter WC not found");
    }
    Rcpp::List WC = Rcpp::as<Rcpp::List>(pars["WC"]);
    theta_s_ = 1 - rho_b_ / rho_p_;
    phi_s_ = theta_s_ * z_;
    theta_w_ = 0.4 * ((m_ * rho_b_) / 100);
    phi_w_ = theta_w_ * z_;
    theta_r_ = phi_w_ / phi_s_;
    theta_f_ = interp1(WC["ClayContent"], WC["FC"], m_);
    phi_f_ = theta_f_ * z_;
}

Rcpp::List Bucket::cali_pars()
{
    return Rcpp::List::create(Rcpp::Named("layer") = layer_,
                              Rcpp::Named("z") = z_,
                              Rcpp::Named("rho_b") = rho_b_,
                              Rcpp::Named("m") = m_,
                              Rcpp::Named("rho_p") = rho_p_,
                              Rcpp::Named("Ks") = Ks_,
                              Rcpp::Named("n") = n_,
                              Rcpp::Named("pE") = pE_,
                              Rcpp::Named("pZ") = pZ_);
}