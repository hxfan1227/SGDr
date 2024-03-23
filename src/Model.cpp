#include "Model.h"

Model::Model(const Rcpp::NumericVector &Time,
             const Rcpp::NumericVector &Precipitation,
             const Rcpp::NumericVector &Evaporation,
             Rcpp::NumericVector H2OinSB1,
             Rcpp::NumericVector H2OinSB2,
             const Rcpp::List &calibratableParams,
             const Rcpp::List &constParams) : parameters(calibratableParams, constParams),
                                              simLength(Time.size()),
                                              R(Precipitation),
                                              E0(Evaporation),
                                              H2O_SB1(H2OinSB1),
                                              H2O_SB2(H2OinSB2)
{

    initializeVector();
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
}

Rcpp::DataFrame Model::calc_recharge()
{
    for (int i = 0; i < simLength; ++i)
    {
        if (i > 0)
        {
            H2O_SB1[i] = H2O3_SB1[i - 1];
            H2O_SB2[i] = H2O3_SB2[i - 1];
        }
        // curve number
        SW[i] = H2O_SB1[i] + H2O_SB2[i] - (parameters.get_bucket1().get_WPmm() + parameters.get_bucket2().get_WPmm());
        S[i] = parameters.get_curveNumber().get_Smax() * (1 - SW[i] / (SW[i] + exp(parameters.get_curveNumber().get_w1() - parameters.get_curveNumber().get_w2() * SW[i])));
        Ia[i] = S[i] * 0.2;
        CN[i] = 25400 / (S[i] + 254);
        Qsurf[i] = (R[i] < Ia[i]) ? 0 : pow((R[i] - Ia[i]), 2) / (R[i] + 0.8 * S[i]);
        // bucket1
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
        // bucket2
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
        // aquifer
        Wrechg[i] = i == 0 ? (1 - exp(-1 / parameters.get_aquifer().get_delta())) * Wperc2[i] : (1 - exp(-1 / parameters.get_aquifer().get_delta())) * Wperc2[i] + exp(-1 / parameters.get_aquifer().get_delta()) * Wrechg[i - 1];
    }
    // return Wrechg;
    Rcpp::DataFrame results = Rcpp::DataFrame::create(
        Rcpp::Named("H2O3_SB1") = H2O3_SB1,
        Rcpp::Named("H2O3_SB2") = H2O3_SB2,
        Rcpp::Named("Wperc2") = Wperc2,
        Rcpp::Named("Wrechg") = Wrechg,
        Rcpp::Named("SW") = SW,
        Rcpp::Named("S") = S,
        Rcpp::Named("Ia") = Ia,
        Rcpp::Named("CN") = CN,
        Rcpp::Named("Qsurf") = Qsurf,
        Rcpp::Named("finf_SB1") = finf_SB1,
        Rcpp::Named("finfla_SB1") = finfla_SB1);
    return results;
}

RCPP_MODULE(ModelModule)
{
    Rcpp::class_<Model>("Model")
        .constructor<Rcpp::NumericVector, Rcpp::NumericVector, Rcpp::NumericVector, Rcpp::NumericVector, Rcpp::NumericVector, Rcpp::List, Rcpp::List>()
        .method("calc_recharge", &Model::calc_recharge);
}