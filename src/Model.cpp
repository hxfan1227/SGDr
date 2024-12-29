#include "Model.h"

Model::Model(Rcpp::DataFrame inputData,
             const Rcpp::List &calibratableParams,
             const Rcpp::List &constParams, int warmUp) : inputData(inputData),
                                                          parameters(calibratableParams, constParams),
                                                          simLength(inputData.nrows()),
                                                          warmUp(warmUp),
                                                          R(Rcpp::as<Rcpp::NumericVector>(inputData["R"])),
                                                          E0(Rcpp::as<Rcpp::NumericVector>(inputData["E0"])),
                                                          q1(Rcpp::as<Rcpp::NumericVector>(inputData["Pumping"])),
                                                          hf(Rcpp::as<Rcpp::NumericVector>(inputData["GWL"])),
                                                          phi1(Rcpp::as<Rcpp::NumericVector>(inputData["H2O_SB1"])),
                                                          phi2(Rcpp::as<Rcpp::NumericVector>(inputData["H2O_SB2"]))

{
}

void Model::calc(int nw, int period)
{
    for (int i = 0; i < period; ++i)
    {
        bucket1.update(i);
        bucket2.update(i);
        cn.calc_runoff(i, R, bucket1, bucket2);
        bucket1.calc_phi(i, R, E0, cn, bucket2);
        bucket2.calc_phi(i, R, E0, cn, bucket1);
        aquifer.calc_recharge(i, bucket2);
    }
    aquifer.average_recharge(nw);
    for (int i = 0; i < period; ++i)
    {
        aquifer.update(i);
        aquifer.calc_sfgd(i, constpar, q1);
    }
}

void Model::run(int nw, int warmup)
{
    bucket1.initialize(simLength, phi1[0]);
    bucket2.initialize(simLength, phi2[0]);
    aquifer.initialize(simLength, hf[0]);
    cn.initialize(simLength);

    if (warmup >= 1)
    {
        calc(nw, warmup);
        bucket1.reset(warmup);
        bucket2.reset(warmup);
        aquifer.reset(warmup);
        calc(nw, simLength);
    }
    else
    {
        calc(nw, simLength);
    }
}

Rcpp::DataFrame Model::get_recharge_output()
{
    return Rcpp::DataFrame::create(
        Rcpp::Named("H2O3_SB1") = bucket1._phi(),
        Rcpp::Named("H2O3_SB2") = bucket2._phi(),
        Rcpp::Named("Wperc2") = bucket2.P(),
        // Wrechg
        Rcpp::Named("Wrechg") = aquifer.Wnet(),
        // SW
        Rcpp::Named("SW") = cn.Phi(),
        // S
        Rcpp::Named("S") = cn.S(),
        // Ia
        Rcpp::Named("Ia") = cn.Ia(),
        // CN
        Rcpp::Named("CN") = cn.CN(),
        // Qsurf
        Rcpp::Named("Qsurf") = cn.Q(),
        // finf_SB1
        Rcpp::Named("finf_SB1") = bucket1.F(),
        // finfla_SB1
        Rcpp::Named("finfla_SB1") = bucket1.FI());
}

Rcpp::NumericVector Model::get_Qsurf()
{
    return cn.Q();
}


Rcpp::DataFrame Model::get_sgd_output()
{
    Rcpp::NumericVector t_temp = inputData["t"];
    t_temp = t_temp;
    return Rcpp::DataFrame::create(
        Rcpp::Named("t") = t_temp,
        Rcpp::Named("wl") = aquifer._hf(),
        Rcpp::Named("SGD") = aquifer.q0(),
        Rcpp::Named("xn") = aquifer.xn(),
        Rcpp::Named("hn") = aquifer.hn(),
        Rcpp::Named("Wrechg") = aquifer.Wnet(),
        Rcpp::Named("WrechgAve") = aquifer._Wnet());
}

Rcpp::List Model::get_all_params_list()
{
    return parameters.get_all_params_list();
}
Rcpp::DataFrame Model::get_inputData()
{
    Rcpp::NumericVector t_temp = inputData["t"];
    Rcpp::NumericVector R_temp = inputData["R"];
    Rcpp::NumericVector E0_temp = inputData["E0"];
    Rcpp::NumericVector Pump_temp = inputData["Pumping"];
    return Rcpp::DataFrame::create(
        Rcpp::Named("t") = t_temp,
        Rcpp::Named("R") = R_temp,
        Rcpp::Named("E0") = E0_temp,
        Rcpp::Named("Pumping") = Pump_temp);
}

// RCPP_EXPOSED_CLASS(Model);