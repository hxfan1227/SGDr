#include "aquifer.h"
#include "constpar.h"

class Bucket;

void Aquifer::initialize(int sim_length, double init_hf)
{
    Wnet_ = Rcpp::NumericVector(sim_length, 0.0);
    _Wnet_ = Rcpp::NumericVector(sim_length, 0.0);
    hf_ = Rcpp::NumericVector(sim_length, init_hf);
    _hf_ = Rcpp::NumericVector(sim_length, 0.0);
    q0_ = Rcpp::NumericVector(sim_length, 0.0);
    q01_ = Rcpp::NumericVector(sim_length, 0.0);
    q02_ = Rcpp::NumericVector(sim_length, 0.0);
    xn_ = Rcpp::NumericVector(sim_length, 0.0);
    xn1_ = Rcpp::NumericVector(sim_length, 0.0);
    xn2_ = Rcpp::NumericVector(sim_length, 0.0);
    hn_ = Rcpp::NumericVector(sim_length, 0.0);
    hn1_ = Rcpp::NumericVector(sim_length, 0.0);
    hn2_ = Rcpp::NumericVector(sim_length, 0.0);
    M1_ = Rcpp::NumericVector(sim_length, 0.0);
    M2_ = Rcpp::NumericVector(sim_length, 0.0);
    xT1_ = Rcpp::NumericVector(sim_length, 0.0);
    xT2_ = Rcpp::NumericVector(sim_length, 0.0);
    dh0_ = Rcpp::NumericVector(sim_length, 0.0);
    dh1_ = Rcpp::NumericVector(sim_length, 0.0);
    dh2_ = Rcpp::NumericVector(sim_length, 0.0);
    dxT_ = Rcpp::NumericVector(sim_length, 0.0);
    _XT_ = Rcpp::NumericVector(sim_length, 0.0);
    V_ = Rcpp::NumericVector(sim_length, 0.0);
}

void Aquifer::reset(int warmup)
{
    Wnet_[0] = Wnet_[warmup - 1];
    _Wnet_[0] = _Wnet_[warmup - 1];
    hf_[0] = hf_[warmup - 1];
    _hf_[0] = _hf_[warmup - 1];
    q0_[0] = q0_[warmup - 1];
    q01_[0] = q01_[warmup - 1];
    q02_[0] = q02_[warmup - 1];
    xn_[0] = xn_[warmup - 1];
    xn1_[0] = xn1_[warmup - 1];
    xn2_[0] = xn2_[warmup - 1];
    hn_[0] = hn_[warmup - 1];
    hn1_[0] = hn1_[warmup - 1];
    hn2_[0] = hn2_[warmup - 1];
    M1_[0] = M1_[warmup - 1];
    M2_[0] = M2_[warmup - 1];
    xT1_[0] = xT1_[warmup - 1];
    xT2_[0] = xT2_[warmup - 1];
    dh0_[0] = dh0_[warmup - 1];
    dh1_[0] = dh1_[warmup - 1];
    dh2_[0] = dh2_[warmup - 1];
    dxT_[0] = dxT_[warmup - 1];
    _XT_[0] = _XT_[warmup - 1];
    V_[0] = V_[warmup - 1];
}

Aquifer::Aquifer(const Rcpp::List aquiferParams)
{
    if (!aquiferParams.containsElementNamed("delta"))
    {
        Rcpp::stop("Aquifer parameter delta not found");
    }
    Td_ = Rcpp::as<double>(aquiferParams["delta"]);
    if (!aquiferParams.containsElementNamed("rho_s"))
    {
        Rcpp::stop("Aquifer parameter rho_s not found");
    }
    rho_s_ = Rcpp::as<double>(aquiferParams["rho_s"]);
    if (!aquiferParams.containsElementNamed("rho_f"))
    {
        Rcpp::stop("Aquifer parameter rho_f not found");
    }
    rho_f_ = Rcpp::as<double>(aquiferParams["rho_f"]);
    if (!aquiferParams.containsElementNamed("K"))
    {
        Rcpp::stop("Aquifer parameter K not found");
    }
    Ka_ = Rcpp::as<double>(aquiferParams["K"]);
    if (!aquiferParams.containsElementNamed("z0"))
    {
        Rcpp::stop("Aquifer parameter z0 not found");
    }
    z0_ = Rcpp::as<double>(aquiferParams["z0"]);
    if (!aquiferParams.containsElementNamed("Sy"))
    {
        Rcpp::stop("Aquifer parameter Sy not found");
    }
    Sy_ = Rcpp::as<double>(aquiferParams["Sy"]);
    if (!aquiferParams.containsElementNamed("xT"))
    {
        Rcpp::stop("Aquifer parameter xT not found");
    }
    XT_ = Rcpp::as<double>(aquiferParams["xT"]);
    if (!aquiferParams.containsElementNamed("dxT"))
    {
        Rcpp::stop("Aquifer parameter dxT not found");
    }
    dxT_max_ = Rcpp::as<double>(aquiferParams["dxT"]);
    a_ = (rho_s_ - rho_f_) / rho_f_;
    he_ = z0_ * rho_s_ / rho_f_ - z0_;
}
double Aquifer::Td() { return Td_; }
double Aquifer::rho_s() { return rho_s_; }
double Aquifer::rho_f() { return rho_f_; }
double Aquifer::Ka() { return Ka_; }
double Aquifer::z0() { return z0_; }
double Aquifer::Sy() { return Sy_; }
double Aquifer::get_a() { return a_; }
double Aquifer::get_he() { return he_; }
double Aquifer::XT() { return XT_; }
double Aquifer::dxT_max() { return dxT_max_; }

double Aquifer::Wnet(int i) { return Wnet_[i]; }
double Aquifer::_Wnet(int i) { return _Wnet_[i]; }
double Aquifer::hf(int i) { return hf_[i]; }
double Aquifer::_hf(int i) { return _hf_[i]; }
double Aquifer::q0(int i) { return q0_[i]; }

double Aquifer::__q01(int i, ConstParameter &constpar)
{
    q01_[i] = ((Ka_ / 2 / constpar.get_x() * rho_s_ / (rho_s_ - rho_f_) * std::pow(hf_[i], 2)) + (_Wnet_[i] / 1000) * constpar.get_x() / 2);
    return q01_[i];
}

double Aquifer::__q02(int i, ConstParameter &constpar)
{
    q02_[i] = (Ka_ * ((std::pow((hf_[i] + z0_), 2) - rho_s_ / rho_f_ * std::pow(z0_, 2))) + (_Wnet_[i] / 1000) * std::pow(constpar.get_x(), 2)) / 2 / constpar.get_x();
    return q02_[i];
}

double Aquifer::xn(int i) { return xn_[i]; }

double Aquifer::__xn1(int i, ConstParameter &constpar)
{
    if (_Wnet_[i] == 0)
    {
        xn1_[i] = 1000000000;
    }
    else
    {
        xn1_[i] = q01_[i] / (_Wnet_[i] / 1000);
    }
    return xn1_[i];
}
double Aquifer::__xn2(int i, ConstParameter &constpar)
{
    if (_Wnet_[i] == 0)
    {
        xn2_[i] = 1000000000;
    }
    else
    {
        xn2_[i] = q02_[i] / (_Wnet_[i] / 1000);
    }
    return xn2_[i];
}

double Aquifer::hn(int i) { return hn_[i]; }
double Aquifer::__hn1(int i, ConstParameter &constpar)
{
    if (_Wnet_[i] == 0)
    {
        hn1_[i] = 1000;
    }
    else
    {
        hn1_[i] = sqrt(std::pow(q01_[i], 2) / (_Wnet_[i] / 1000) / Ka_ + rho_s_ / rho_f_ * z0_ * z0_) - z0_;
    }
    return hn1_[i];
}
double Aquifer::__hn2(int i, ConstParameter &constpar)
{
    if (_Wnet_[i] == 0)
    {
        hn2_[i] = 1000;
    }
    else
    {
        hn2_[i] = sqrt(std::pow(q02_[i], 2) / (_Wnet_[i] / 1000) / Ka_ + rho_s_ / rho_f_ * z0_ * z0_) - z0_;
    }
    return hn2_[i];
}
double Aquifer::__M1(int i, ConstParameter &constpar)
{
    if (_Wnet_[i] == 0)
    {
        M1_[i] = 0;
    }
    else
    {
        M1_[i] = Ka_ * a_ * (1 + a_) * z0_ * z0_ / (_Wnet_[i] / 1000) / (xn1_[i] * xn1_[i]);
    }
    return M1_[i];
}
double Aquifer::__M2(int i, ConstParameter &constpar)
{
    if (_Wnet_[i] == 0)
    {
        M2_[i] = 0;
    }
    else
    {
        M2_[i] = Ka_ * a_ * (1 + a_) * z0_ * z0_ / (_Wnet_[i] / 1000) / (xn2_[i] * xn2_[i]);
    }
    return M2_[i];
}
double Aquifer::__xT1(int i, ConstParameter &constpar)
{
    if (_Wnet_[i] == 0)
    {
        xT1_[i] = Ka_ * a_ * (1 + a_) * z0_ * z0_ / 2 / (q01_[i]);
    }
    else
    {
        xT1_[i] = (q01_[i]) / (_Wnet_[i] / 1000) - sqrt(pow((q01_[i]) / (_Wnet_[i] / 1000), 2) - Ka_ * a_ * (1 + a_) * (z0_ * z0_) / (_Wnet_[i] / 1000));
    }
    return xT1_[i];
}
double Aquifer::__xT2(int i, ConstParameter &constpar)
{
    if (_Wnet_[i] == 0)
    {
        xT2_[i] = Ka_ * a_ * (1 + a_) * z0_ * z0_ / 2 / (q02_[i]);
    }
    else
    {
        xT2_[i] = (q02_[i]) / (_Wnet_[i] / 1000) - sqrt(pow((q02_[i]) / (_Wnet_[i] / 1000), 2) - Ka_ * a_ * (1 + a_) * (z0_ * z0_) / (_Wnet_[i] / 1000));
    }
    return xT2_[i];
}
double Aquifer::dh0(int i) { return dh0_[i]; }
double Aquifer::dh1(int i) { return dh1_[i]; }
double Aquifer::dh2(int i) { return dh2_[i]; }
double Aquifer::dxT(int i) { return dxT_[i]; }
double Aquifer::_XT(int i) { return _XT_[i]; }

Rcpp::List Aquifer::get_all_params_list()
{
    return Rcpp::List::create(Rcpp::Named("delta") = Td_,
                              Rcpp::Named("rho_s") = rho_s_,
                              Rcpp::Named("rho_f") = rho_f_,
                              Rcpp::Named("K") = Ka_,
                              Rcpp::Named("z0") = z0_,
                              Rcpp::Named("Sy") = Sy_,
                              Rcpp::Named("a") = a_,
                              Rcpp::Named("he") = he_,
                              Rcpp::Named("xT") = XT_,
                              Rcpp::Named("dxT") = dxT_max_);
}

Rcpp::List Aquifer::get_calibratable_params_list()
{
    return Rcpp::List::create(Rcpp::Named("delta") = Td_,
                              Rcpp::Named("K") = Ka_,
                              Rcpp::Named("z0") = z0_,
                              Rcpp::Named("Sy") = Sy_,
                              Rcpp::Named("xT") = XT_,
                              Rcpp::Named("dxT") = dxT_max_);
}

void Aquifer::calc_recharge(int i, Bucket &bucket)
{
    Wnet_[i] = i == 0 ? (1 - exp(-1 / Td_)) * bucket.P(i) : (1 - exp(-1 / Td_)) * bucket.P(i) + exp(-1 / Td_) * Wnet_[i - 1];
}

void Aquifer::update(int i)
{
    if (i > 0)
    {
        hf_[i] = _hf_[i - 1];
    }
}

Rcpp::NumericVector Aquifer::moving_average(Rcpp::NumericVector &x, int nw)
{
    int n = x.size();
    Rcpp::NumericVector avg(n); // Initialize the result vector

    // Initial pass for moving average calculation
    for (int i = 0; i < n; ++i)
    {
        if (i < nw - 1)
        {
            avg[i] = NA_REAL; // NA for elements before enough data points are available
        }
        else
        {
            double sum = 0;
            for (int j = i - nw + 1; j <= i; ++j)
            {
                sum += x[j];
            }
            avg[i] = sum / nw;
        }
    }

    // Fill initial NAs with the first available average value
    if (n >= nw)
    {
        double firstAvg = avg[nw - 1];
        for (int i = 0; i < nw - 1; ++i)
        {
            avg[i] = firstAvg;
        }
    }

    return avg;
}

void Aquifer::average_recharge(int nw)
{
    _Wnet_ = moving_average(Wnet_, nw);
}

// void Aquifer::calc_sfgd(int sim_length, int nw, ConstParameter &constpar, Rcpp::NumericVector &q1)
// {
//     _Wnet_ = moving_average(Wnet_, nw);
//     for (int i = 0; i < sim_length; ++i)
//     {
//         update(i);
//         calc_sfgds(i, constpar, q1);
//     }
// }

void Aquifer::calc_sfgd(int i, ConstParameter &constpar, Rcpp::NumericVector &q1)
{
    __q01(i, constpar);
    __q02(i, constpar);
    __xn1(i, constpar);
    __xn2(i, constpar);
    __hn1(i, constpar);
    __hn2(i, constpar);
    __M1(i, constpar);
    __M2(i, constpar);
    if (M1_[i] < 1)
    {
        __xT1(i, constpar);
    }
    if (M2_[i] < 1)
    {
        __xT2(i, constpar);
    }

    if (xT2_[i] < 0)
    {
        q0_[i] = q01_[i];
        xn_[i] = xn1_[i];
        hn_[i] = hn1_[i];
    }
    else
    {
        if (M1_[i] >= 1)
        {
            q0_[i] = q01_[i];
            xn_[i] = xn1_[i];
            hn_[i] = hn1_[i];
        }
        else
        {
            if (constpar.get_W() <= xT2_[i])
            {
                q0_[i] = q01_[i];
                xn_[i] = xn1_[i];
                hn_[i] = hn1_[i];
            }
            else
            {
                if (constpar.get_x() > xT2_[i])
                {
                    q0_[i] = q02_[i];
                    xn_[i] = xn2_[i];
                    hn_[i] = hn2_[i];
                }
                else
                {
                    q0_[i] = q01_[i];
                    xn_[i] = xn1_[i];
                    hn_[i] = hn1_[i];
                }
            }
        }
    }
    // update heads
    _hf_[i] = hf_[i] + Wnet_[i] / 1000 / Sy_;
    dh0_[i] = q0_[i] * constpar.get_W() / constpar.get_Area() * 1000 / Sy_;
    dh1_[i] = q1[i] / constpar.get_Area() * 1000 / Sy_;
    _hf_[i] = _hf_[i] - dh0_[i] / 1000 - dh1_[i] / 1000;
    if (i == 0)
    {
        dxT_[i] = 0;
        _XT_[i] = XT_;
    }
    else
    {
        if (M1_[i] > 1)
        {
            dxT_[i] = dxT_max_;
        }
        else
        {
            if (xT1_[i] >= _XT_[i - 1])
            {
                dxT_[i] = std::min(dxT_max_, xT1_[i] - _XT_[i - 1]);
            }
            else
            {
                dxT_[i] = std::max(-dxT_max_, xT1_[i] - _XT_[i - 1]);
            }
        }
        _XT_[i] = _XT_[i - 1] + dxT_[i];
    }
    V_[i] = Sy_ * z0_ * _XT_[i] / 3;
    dh2_[i] = i == 0 ? 0 : (V_[i - 1] - V_[i]) * constpar.get_W() / constpar.get_Area() * 1000;
    _hf_[i] = i == 0 ? _hf_[i] : _hf_[i] - dh2_[i] / 1000;
}

Rcpp::NumericVector Aquifer::Wnet() { return Wnet_; }
Rcpp::NumericVector Aquifer::_Wnet() { return _Wnet_; }
Rcpp::NumericVector Aquifer::hf() { return hf_; }
Rcpp::NumericVector Aquifer::_hf() { return _hf_; }
Rcpp::NumericVector Aquifer::q0() { return q0_; }
Rcpp::NumericVector Aquifer::xn() { return xn_; }
Rcpp::NumericVector Aquifer::hn() { return hn_; }
Rcpp::NumericVector Aquifer::dh0() { return dh0_; }
Rcpp::NumericVector Aquifer::dh1() { return dh1_; }
Rcpp::NumericVector Aquifer::dh2() { return dh2_; }
Rcpp::NumericVector Aquifer::dxT() { return dxT_; }
Rcpp::NumericVector Aquifer::_XT() { return _XT_; }
Rcpp::NumericVector Aquifer::V() { return V_; }