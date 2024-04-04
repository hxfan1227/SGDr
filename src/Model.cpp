#include "Model.h"

Model::Model(Rcpp::DataFrame inputData,
             const Rcpp::List &calibratableParams,
             const Rcpp::List &constParams) : inputData(inputData),
                                              parameters(calibratableParams, constParams),
                                              simLength(inputData.nrows()),
                                              R(Rcpp::as<Rcpp::NumericVector>(inputData["R"])),
                                              E0(Rcpp::as<Rcpp::NumericVector>(inputData["E0"])),
                                              Pumping(Rcpp::as<Rcpp::NumericVector>(inputData["Pumping"])),
                                              GWL(Rcpp::as<Rcpp::NumericVector>(inputData["GWL"])),
                                              H2O_SB1(Rcpp::as<Rcpp::NumericVector>(inputData["H2O_SB1"])),
                                              H2O_SB2(Rcpp::as<Rcpp::NumericVector>(inputData["H2O_SB2"]))

{
}
void Model::initializeVector()
{
    // Initialize the vectors with NA values
    SW = Rcpp::NumericVector(simLength, NA_REAL);
    S = Rcpp::NumericVector(simLength, NA_REAL);
    Ia = Rcpp::NumericVector(simLength, NA_REAL);
    CN = Rcpp::NumericVector(simLength, NA_REAL);
    Qsurf = Rcpp::NumericVector(simLength, NA_REAL);
    MC_SB1 = Rcpp::NumericVector(simLength, NA_REAL);
    DoS_SB1 = Rcpp::NumericVector(simLength, NA_REAL);
    finf_SB1 = Rcpp::NumericVector(simLength, NA_REAL);
    finfla_SB1 = Rcpp::NumericVector(simLength, NA_REAL);
    H2O1_SB1 = Rcpp::NumericVector(simLength, NA_REAL);
    IntercH2O = Rcpp::NumericVector(simLength, NA_REAL);
    CanopyH2O = Rcpp::NumericVector(simLength, NA_REAL);
    E0_Int = Rcpp::NumericVector(simLength, NA_REAL);
    YE = Rcpp::NumericVector(simLength, NA_REAL);
    E = Rcpp::NumericVector(simLength, NA_REAL);
    H2O2_SB1 = Rcpp::NumericVector(simLength, NA_REAL);
    QSE = Rcpp::NumericVector(simLength, NA_REAL);
    SWexcess = Rcpp::NumericVector(simLength, NA_REAL);
    Kunsat = Rcpp::NumericVector(simLength, NA_REAL);
    Kunsat2 = Rcpp::NumericVector(simLength, NA_REAL);
    TT = Rcpp::NumericVector(simLength, NA_REAL);
    TT2 = Rcpp::NumericVector(simLength, NA_REAL);
    Wdperc = Rcpp::NumericVector(simLength, NA_REAL);
    Wperc = Rcpp::NumericVector(simLength, NA_REAL);
    Wperc2 = Rcpp::NumericVector(simLength, NA_REAL);
    H2O3_SB1 = Rcpp::NumericVector(simLength, NA_REAL);
    H2O3_SB2 = Rcpp::NumericVector(simLength, NA_REAL);
    MC_SB2 = Rcpp::NumericVector(simLength, NA_REAL);
    DoS_SB2 = Rcpp::NumericVector(simLength, NA_REAL);
    H2O1_SB2 = Rcpp::NumericVector(simLength, NA_REAL);
    YE2 = Rcpp::NumericVector(simLength, NA_REAL);
    Edd = Rcpp::NumericVector(simLength, NA_REAL);
    H2O2_SB2 = Rcpp::NumericVector(simLength, NA_REAL);
    Sw = Rcpp::NumericVector(simLength, NA_REAL);
    Sw2 = Rcpp::NumericVector(simLength, NA_REAL);
    SWexcess2 = Rcpp::NumericVector(simLength, NA_REAL);
    Wrechg = Rcpp::NumericVector(simLength, NA_REAL);
    H2O1_AQ = Rcpp::NumericVector(simLength, NA_REAL);
    WrechgAve = Rcpp::NumericVector(simLength, NA_REAL);
    H2O2_AQ = Rcpp::NumericVector(simLength, NA_REAL);
    H2O3_AQ = Rcpp::NumericVector(simLength, NA_REAL);
    SGD1 = Rcpp::NumericVector(simLength, NA_REAL);
    SGD2 = Rcpp::NumericVector(simLength, NA_REAL);
    xn = Rcpp::NumericVector(simLength, NA_REAL);          // the distance from the coast to the peak of the groundwater mound
    xn1 = Rcpp::NumericVector(simLength, NA_REAL);         // the distance from the coast to the peak of the groundwater mound when x < xT
    xn2 = Rcpp::NumericVector(simLength, NA_REAL);         // the distance from the coast to the peak of the groundwater mound when x > xT
    hn = Rcpp::NumericVector(simLength, NA_REAL);          // the elevation of the groundwater mount (0m AHD)
    hn1 = Rcpp::NumericVector(simLength, NA_REAL);         // the elevation of the groundwater mount (0m AHD) when x < xT
    hn2 = Rcpp::NumericVector(simLength, NA_REAL);         // the elevation of the groundwater mount (0m AHD) when x > xT
    M1 = Rcpp::NumericVector(simLength, NA_REAL);          //
    M2 = Rcpp::NumericVector(simLength, NA_REAL);          //
    xT1 = Rcpp::NumericVector(simLength, NA_REAL);         // positioning of the saltwater toe when x < xT
    xT2 = Rcpp::NumericVector(simLength, NA_REAL);         // positioning of the saltwater toe when x > xT
    SGD = Rcpp::NumericVector(simLength, NA_REAL);         // SGD calculated based on the calculate position of the saltwater toe
    SGDdrop = Rcpp::NumericVector(simLength, NA_REAL);     // drop in water level of the bucket due to SGD (m)
    PumpingDrop = Rcpp::NumericVector(simLength, NA_REAL); // drop in water level of the bucket due to pumping extraction (m)
    dxT = Rcpp::NumericVector(simLength, NA_REAL);         // saltwater toe moving speed m/day
    xT3 = Rcpp::NumericVector(simLength, NA_REAL);         // the position of the saltwater toe at the end of the day
    SWvol = Rcpp::NumericVector(simLength, NA_REAL);       // volume of soil water in the aquifer bucket
    FWLdrop = Rcpp::NumericVector(simLength, NA_REAL);     // change in water level in the bucket
}

void Model::calc_recharge()
{
    initializeVector();
    for (int i = 0; i < simLength; ++i)
    {
        update_head(i);
        // curve number
        calc_cn_runoff(i);
        // bucket1
        calc_bucket1_h2o(i);
        // bucket2
        calc_bucket2_h2o(i);
        // aquifer
        calc_aquifer_recharge(i);
    }
}

void Model::calc_sgd(int windowSize)
{
    WrechgAve = moving_average(Wrechg, windowSize);
    for (int i = 0; i < simLength; ++i)
    {
        update_gwl(i);
        calc_sgds(i);
        update_aq_head(i);
    }
}

void Model::update_head(int i)
{
    if (i > 0)
    {
        H2O_SB1[i] = H2O3_SB1[i - 1];
        H2O_SB2[i] = H2O3_SB2[i - 1];
    }
}

void Model::calc_cn_runoff(int i)
{
    SW[i] = H2O_SB1[i] + H2O_SB2[i] - (parameters.get_bucket1().get_WPmm() + parameters.get_bucket2().get_WPmm());
    S[i] = parameters.get_curveNumber().get_Smax() * (1 - SW[i] / (SW[i] + exp(parameters.get_curveNumber().get_w1() - parameters.get_curveNumber().get_w2() * SW[i])));
    Ia[i] = S[i] * 0.2;
    CN[i] = 25400 / (S[i] + 254);
    Qsurf[i] = (R[i] < Ia[i]) ? 0 : pow((R[i] - Ia[i]), 2) / (R[i] + 0.8 * S[i]);
}

void Model::calc_bucket1_h2o(int i)
{
    MC_SB1[i] = H2O_SB1[i] / parameters.get_bucket1().get_z();
    DoS_SB1[i] = MC_SB1[i] / parameters.get_bucket1().get_phi_soil();
    finf_SB1[i] = (R[i] < Ia[i]) ? 0 : R[i] - Qsurf[i] - Ia[i];
    finfla_SB1[i] = R[i] > Ia[i] ? parameters.get_curveNumber().get_PIa() * Ia[i] : (R[i] > (Ia[i] * (1 - parameters.get_curveNumber().get_PIa())) ? R[i] - (1 - parameters.get_curveNumber().get_PIa()) * Ia[i] : 0);
    H2O1_SB1[i] = H2O_SB1[i] + finf_SB1[i] + finfla_SB1[i];
    IntercH2O[i] = R[i] > ((1 - parameters.get_curveNumber().get_PIa()) * Ia[i]) ? (1 - parameters.get_curveNumber().get_PIa()) * Ia[i] : R[i];
    CanopyH2O[i] = i == 0 ? 0.0 : std::max(CanopyH2O[i - 1] + IntercH2O[i] - E0[i], 0.0);
    E0_Int[i] = std::max(0.0, E0[i] - IntercH2O[i] - CanopyH2O[i]);
    YE[i] = E0_Int[i] * parameters.get_bucket1().get_Y() / 100;
    E[i] = std::min(H2O1_SB1[i] < parameters.get_bucket1().get_FCmm() ? YE[i] * exp(2.5 * (H2O1_SB1[i] - parameters.get_bucket1().get_FCmm()) / (parameters.get_bucket1().get_FCmm() - parameters.get_bucket1().get_WPmm())) : YE[i], 0.8 * (H2O1_SB1[i] - parameters.get_bucket1().get_WPmm()));
    H2O2_SB1[i] = std::min(std::max(0.0, H2O1_SB1[i] - E[i]), parameters.get_bucket1().get_SAT());
    QSE[i] = ((H2O1_SB1[i] - E[i]) > parameters.get_bucket1().get_SAT()) ? H2O1_SB1[i] - E[i] - parameters.get_bucket1().get_SAT() : 0;
    SWexcess[i] = std::max(0.0, H2O2_SB1[i] - parameters.get_bucket1().get_FCmm());
    Sw[i] = std::min(1.0, (H2O2_SB1[i] / parameters.get_bucket1().get_SAT() - parameters.get_bucket1().get_Swres()) / (1 - parameters.get_bucket1().get_Swres()));
    Kunsat[i] = std::max(0.0,
                         pow(Sw[i], 0.5) * pow((1 - pow(1 - pow(Sw[i], parameters.get_bucket1().get_n() / (parameters.get_bucket1().get_n() - 1)), (parameters.get_bucket1().get_n() - 1) / parameters.get_bucket1().get_n())), 2) * parameters.get_bucket1().get_Ksat());
    TT[i] = std::max((parameters.get_bucket1().get_SAT() - parameters.get_bucket1().get_FCmm()) / Kunsat[i], 72.0);
    Wdperc[i] = SWexcess[i] * (1 - exp(-24 / TT[i]));
    Wperc[i] = H2O_SB2[i] >= parameters.get_bucket2().get_FCmm() + parameters.get_bucket2().get_Z() / 100 * (parameters.get_bucket2().get_SAT() - parameters.get_bucket2().get_FCmm()) ? 0.0 : (Wdperc[i] <= (parameters.get_bucket2().get_SAT() - H2O_SB2[i]) ? Wdperc[i] : (parameters.get_bucket2().get_SAT() - H2O_SB2[i]));
    H2O3_SB1[i] = H2O2_SB1[i] - Wperc[i];
}

void Model::calc_bucket2_h2o(int i)
{
    MC_SB2[i] = H2O_SB2[i] / parameters.get_bucket2().get_z();
    DoS_SB2[i] = MC_SB2[i] / parameters.get_bucket2().get_phi_soil();
    H2O1_SB2[i] = H2O_SB2[i] + Wperc[i];
    YE2[i] = E0_Int[i] * (100 - parameters.get_bucket1().get_Y()) / 100;
    Edd[i] = std::min((H2O1_SB2[i] < parameters.get_bucket2().get_FCmm() ? YE2[i] * exp(2.5 * (H2O1_SB2[i] - parameters.get_bucket2().get_FCmm()) / (parameters.get_bucket2().get_FCmm() - parameters.get_bucket2().get_WPmm())) : YE2[i]),
                      0.8 * (H2O1_SB2[i] - parameters.get_bucket2().get_WPmm()));
    H2O2_SB2[i] = H2O1_SB2[i] - Edd[i];
    SWexcess2[i] = std::max(0.0, H2O1_SB2[i] - parameters.get_bucket2().get_FCmm());
    Sw2[i] = std::min(1.0, (H2O2_SB2[i] / parameters.get_bucket2().get_SAT() - parameters.get_bucket2().get_Swres()) / (1.0 - parameters.get_bucket2().get_Swres()));
    Kunsat2[i] = std::max(0.0,
                          pow(Sw2[i], 0.5) *
                              pow(1 -
                                      pow(1 -
                                              pow(Sw2[i], parameters.get_bucket2().get_n() / (parameters.get_bucket2().get_n() - 1)),
                                          (parameters.get_bucket2().get_n() - 1) / parameters.get_bucket2().get_n()),
                                  2) *
                              parameters.get_bucket2().get_Ksat());
    TT2[i] = (parameters.get_bucket2().get_SAT() - parameters.get_bucket2().get_FCmm()) / Kunsat2[i];
    Wperc2[i] = SWexcess2[i] * (1 - exp(-24 / TT2[i]));
    H2O3_SB2[i] = H2O2_SB2[i] - Wperc2[i];
}

void Model::calc_aquifer_recharge(int i)
{
    Wrechg[i] = i == 0 ? (1 - exp(-1 / parameters.get_aquifer().get_delta())) * Wperc2[i] : (1 - exp(-1 / parameters.get_aquifer().get_delta())) * Wperc2[i] + exp(-1 / parameters.get_aquifer().get_delta()) * Wrechg[i - 1];
}

Rcpp::NumericVector Model::moving_average(Rcpp::NumericVector x, int windowSize)
{
    int n = x.size();
    Rcpp::NumericVector avg(n); // Initialize the result vector

    // Initial pass for moving average calculation
    for (int i = 0; i < n; ++i)
    {
        if (i < windowSize - 1)
        {
            avg[i] = NA_REAL; // NA for elements before enough data points are available
        }
        else
        {
            double sum = 0;
            for (int j = i - windowSize + 1; j <= i; ++j)
            {
                sum += x[j];
            }
            avg[i] = sum / windowSize;
        }
    }

    // Fill initial NAs with the first available average value
    if (n >= windowSize)
    {
        double firstAvg = avg[windowSize - 1];
        for (int i = 0; i < windowSize - 1; ++i)
        {
            avg[i] = firstAvg;
        }
    }

    return avg;
}

void Model::update_gwl(int i)
{
    if (i > 0)
    {
        GWL[i] = H2O2_AQ[i - 1];
    }

    H2O1_AQ[i] = GWL[i] + Wrechg[i] / 1000 / parameters.get_aquifer().get_Sy();
}

void Model::calc_sgds(int i)
{
    // Calculate SGD based on the relative distance between saltwater toe and the observation well location.
    // The equations belwo can handle the case when the average recharge is zero.
    SGD1[i] = calc_sgd1(i);
    SGD2[i] = calc_sgd2(i);
    if (WrechgAve[i] > 0.0)
    {
        xn1[i] = calc_xn(SGD1[i], i);
        xn2[i] = calc_xn(SGD2[i], i);
        hn1[i] = calc_hn(SGD1[i], i);
        hn2[i] = calc_hn(SGD2[i], i);
        M1[i] = calc_M(xn1[i], i);
        M2[i] = calc_M(xn2[i], i);
        if (M1[i] < 1)
        {
            xT1[i] = calc_xT(SGD1[i], i);
        }
        if (M2[i] < 1)
        {
            xT2[i] = calc_xT(SGD2[i], i);
        }
    }
    else // if WrechgAve[i] == 0
    {
        xn1[i] = 1000000000;
        xn2[i] = 1000000000;
        hn1[i] = 1000;
        hn2[i] = 1000;
        M1[i] = 0;
        M2[i] = 0;
        xT1[i] = calc_xT0(SGD1[i], i);
        xT2[i] = calc_xT0(SGD2[i], i);
    }

    if (xT2[i] < 0)
    {
        SGD[i] = SGD1[i];
        xn[i] = xn1[i];
        hn[i] = hn1[i];
    }
    else
    {
        if (M1[i] >= 1)
        {
            SGD[i] = SGD1[i];
            xn[i] = xn1[i];
            hn[i] = hn1[i];
        }
        else
        {
            if (parameters.get_constParameter().get_W() <= xT2[i])
            {
                SGD[i] = SGD1[i];
                xn[i] = xn1[i];
                hn[i] = hn1[i];
            }
            else
            {
                if (parameters.get_constParameter().get_x() > xT2[i])
                {
                    SGD[i] = SGD2[i];
                    xn[i] = xn2[i];
                    hn[i] = hn2[i];
                }
                else
                {
                    SGD[i] = SGD1[i];
                    xn[i] = xn1[i];
                    hn[i] = hn1[i];
                }
            }
        }
    }
}

void Model::update_aq_head(int i)
{

    SGDdrop[i] = SGD[i] / parameters.get_constParameter().get_Area() * 1000 / parameters.get_aquifer().get_Sy();
    PumpingDrop[i] = Pumping[i] / parameters.get_constParameter().get_Area() * 1000 / parameters.get_aquifer().get_Sy();
    H2O2_AQ[i] = H2O1_AQ[i] - SGDdrop[i] / 1000 - PumpingDrop[i] / 1000;
    if (i == 0)
    {
        dxT[i] = 0;
        xT3[i] = parameters.get_aquifer().get_xT();
    }
    else
    {
        if (M1[i] > 1)
        {
            dxT[i] = parameters.get_aquifer().get_dxT();
        }
        else
        {
            if (xT1[i] >= xT3[i - 1])
            {
                dxT[i] = std::min(parameters.get_aquifer().get_dxT(), xT1[i] - xT3[i - 1]);
            }
            else
            {
                dxT[i] = std::max(-parameters.get_aquifer().get_dxT(), xT1[i] - xT3[i - 1]);
            }
        }
        xT3[i] = xT3[i - 1] + dxT[i];
    }
    SWvol[i] = parameters.get_aquifer().get_Sy() * parameters.get_aquifer().get_z0() * xT3[i] / 3;
    FWLdrop[i] = i == 0 ? 0 : (SWvol[i - 1] - SWvol[i]) * parameters.get_constParameter().get_W() / parameters.get_constParameter().get_Area() * 1000;
    H2O3_AQ[i] = i == 0 ? H2O2_AQ[i] : H2O2_AQ[i] - FWLdrop[i] / 1000;
}

double Model::calc_sgd1(int i)
{
    return ((parameters.get_aquifer().get_K() / 2 / parameters.get_constParameter().get_x() * parameters.get_aquifer().get_rho_s() / (parameters.get_aquifer().get_rho_s() - parameters.get_aquifer().get_rho_f()) * std::pow(GWL[i], 2)) + (WrechgAve[i] / 1000) * parameters.get_constParameter().get_x() / 2) * parameters.get_constParameter().get_W();
}

double Model::calc_sgd2(int i)
{
    return ((parameters.get_aquifer().get_K() * ((std::pow((GWL[i] + parameters.get_aquifer().get_z0()), 2) - parameters.get_aquifer().get_rho_s() / parameters.get_aquifer().get_rho_f() * std::pow(parameters.get_aquifer().get_z0(), 2))) + (WrechgAve[i] / 1000) * std::pow(parameters.get_constParameter().get_x(), 2)) / 2 / parameters.get_constParameter().get_x() * parameters.get_constParameter().get_W());
}

double Model::calc_xn(double sgd, int i)
{
    return sgd / (WrechgAve[i] / 1000) / parameters.get_constParameter().get_W();
}

double Model::calc_hn(double sgd, int i)
{
    return sqrt(std::pow(sgd / parameters.get_constParameter().get_W(), 2) / (WrechgAve[i] / 1000) / parameters.get_aquifer().get_K() + parameters.get_aquifer().get_rho_s() / parameters.get_aquifer().get_rho_f() * parameters.get_aquifer().get_z0() * parameters.get_aquifer().get_z0()) - parameters.get_aquifer().get_z0();
}

double Model::calc_M(double xn, int i)
{
    return parameters.get_aquifer().get_K() * parameters.get_aquifer().get_a() * (1 + parameters.get_aquifer().get_a()) * parameters.get_aquifer().get_z0() * parameters.get_aquifer().get_z0() / (WrechgAve[i] / 1000) / (xn * xn);
}

double Model::calc_xT(double sgd, int i)
{
    return (sgd / parameters.get_constParameter().get_W()) / (WrechgAve[i] / 1000) - sqrt(pow((sgd / parameters.get_constParameter().get_W()) / (WrechgAve[i] / 1000), 2) - parameters.get_aquifer().get_K() * parameters.get_aquifer().get_a() * (1 + parameters.get_aquifer().get_a()) * (parameters.get_aquifer().get_z0() * parameters.get_aquifer().get_z0()) / (WrechgAve[i] / 1000));
}

double Model::calc_xT0(double sgd, int i)
{
    return parameters.get_aquifer().get_K() * parameters.get_aquifer().get_a() * (1 + parameters.get_aquifer().get_a()) * parameters.get_aquifer().get_z0() * parameters.get_aquifer().get_z0() / 2 / (sgd / parameters.get_constParameter().get_W());
}

void Model::update_parameters(const Rcpp::List &newCalibratableParams, const Rcpp::List &newConstParams)
{
    parameters.update(newCalibratableParams, newConstParams);
}

// Getters
Rcpp::NumericVector Model::get_SW() { return SW; }
Rcpp::NumericVector Model::get_S() { return S; }
Rcpp::NumericVector Model::get_Ia() { return Ia; }
Rcpp::NumericVector Model::get_CN() { return CN; }
Rcpp::NumericVector Model::get_Qsurf() { return Qsurf; }
Rcpp::NumericVector Model::get_MC_SB1() { return MC_SB1; }
Rcpp::NumericVector Model::get_DoS_SB1() { return DoS_SB1; }
Rcpp::NumericVector Model::get_finf_SB1() { return finf_SB1; }
Rcpp::NumericVector Model::get_finfla_SB1() { return finfla_SB1; }
Rcpp::NumericVector Model::get_H2O1_SB1() { return H2O1_SB1; }
Rcpp::NumericVector Model::get_IntercH2O() { return IntercH2O; }
Rcpp::NumericVector Model::get_CanopyH2O() { return CanopyH2O; }
Rcpp::NumericVector Model::get_E0_Int() { return E0_Int; }
Rcpp::NumericVector Model::get_YE() { return YE; }
Rcpp::NumericVector Model::get_E() { return E; }
Rcpp::NumericVector Model::get_H2O2_SB1() { return H2O2_SB1; }
Rcpp::NumericVector Model::get_QSE() { return QSE; }
Rcpp::NumericVector Model::get_SWexcess() { return SWexcess; }
Rcpp::NumericVector Model::get_Kunsat() { return Kunsat; }
Rcpp::NumericVector Model::get_Kunsat2() { return Kunsat2; }
Rcpp::NumericVector Model::get_TT() { return TT; }
Rcpp::NumericVector Model::get_TT2() { return TT2; }
Rcpp::NumericVector Model::get_Wdperc() { return Wdperc; }
Rcpp::NumericVector Model::get_Wperc() { return Wperc; }
Rcpp::NumericVector Model::get_Wperc2() { return Wperc2; }
Rcpp::NumericVector Model::get_H2O3_SB1() { return H2O3_SB1; }
Rcpp::NumericVector Model::get_H2O3_SB2() { return H2O3_SB2; }
Rcpp::NumericVector Model::get_MC_SB2() { return MC_SB2; }
Rcpp::NumericVector Model::get_DoS_SB2() { return DoS_SB2; }
Rcpp::NumericVector Model::get_H2O1_SB2() { return H2O1_SB2; }
Rcpp::NumericVector Model::get_YE2() { return YE2; }
Rcpp::NumericVector Model::get_Edd() { return Edd; }
Rcpp::NumericVector Model::get_H2O2_SB2() { return H2O2_SB2; }
Rcpp::NumericVector Model::get_Sw() { return Sw; }
Rcpp::NumericVector Model::get_Sw2() { return Sw2; }
Rcpp::NumericVector Model::get_SWexcess2() { return SWexcess2; }
Rcpp::NumericVector Model::get_Wrechg() { return Wrechg; }
Rcpp::NumericVector Model::get_H2O1_AQ() { return H2O1_AQ; }
Rcpp::NumericVector Model::get_WrechgAve() { return WrechgAve; }
Rcpp::NumericVector Model::get_H2O2_AQ() { return H2O2_AQ; }
Rcpp::NumericVector Model::get_H2O3_AQ() { return H2O3_AQ; }
Rcpp::NumericVector Model::get_SGD1() { return SGD1; }
Rcpp::NumericVector Model::get_SGD2() { return SGD2; }
Rcpp::NumericVector Model::get_xn() { return xn; }
Rcpp::NumericVector Model::get_xn1() { return xn1; }
Rcpp::NumericVector Model::get_xn2() { return xn2; }
Rcpp::NumericVector Model::get_hn() { return hn; }
Rcpp::NumericVector Model::get_hn1() { return hn1; }
Rcpp::NumericVector Model::get_hn2() { return hn2; }
Rcpp::NumericVector Model::get_M1() { return M1; }
Rcpp::NumericVector Model::get_M2() { return M2; }
Rcpp::NumericVector Model::get_xT1() { return xT1; }
Rcpp::NumericVector Model::get_xT2() { return xT2; }
Rcpp::NumericVector Model::get_SGD() { return SGD; }
Rcpp::NumericVector Model::get_SGDdrop() { return SGDdrop; }
Rcpp::NumericVector Model::get_PumpingDrop() { return PumpingDrop; }
Rcpp::NumericVector Model::get_dxT() { return dxT; }
Rcpp::NumericVector Model::get_xT3() { return xT3; }
Rcpp::NumericVector Model::get_SWvol() { return SWvol; }
Rcpp::NumericVector Model::get_FWLdrop() { return FWLdrop; }
Rcpp::DataFrame Model::get_recharge_output()
{
    return Rcpp::DataFrame::create(
        Rcpp::Named("H2O3_SB1") = H2O3_SB1,
        Rcpp::Named("H2O3_SB2") = H2O3_SB2,
        Rcpp::Named("Wperc2") = Wperc2,
        // Wrechg
        Rcpp::Named("Wrechg") = Wrechg,
        // SW
        Rcpp::Named("SW") = SW,
        // S
        Rcpp::Named("S") = S,
        // Ia
        Rcpp::Named("Ia") = Ia,
        // CN
        Rcpp::Named("CN") = CN,
        // Qsurf
        Rcpp::Named("Qsurf") = Qsurf,
        // finf_SB1
        Rcpp::Named("finf_SB1") = finf_SB1,
        // finfla_SB1
        Rcpp::Named("finfla_SB1") = finfla_SB1);
}

Rcpp::DataFrame Model::get_sgd_output()
{
    return Rcpp::DataFrame::create(
        Rcpp::Named("t") = inputData["t"],
        Rcpp::Named("wl") = H2O3_AQ,
        Rcpp::Named("SGD") = SGD,
        Rcpp::Named("xn") = xn,
        Rcpp::Named("hn") = hn,
        Rcpp::Named("Wrechg") = Wrechg,
        Rcpp::Named("WrechgAve") = WrechgAve);
}

Rcpp::List Model::get_all_params_list()
{
    return parameters.get_all_params_list();
}
Rcpp::DataFrame Model::get_inputData()
{
    return inputData;
}

// RCPP_EXPOSED_CLASS(Model);