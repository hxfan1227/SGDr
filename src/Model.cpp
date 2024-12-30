#include "Model.h"

Model::Model(Rcpp::DataFrame inputData,
             const Rcpp::List &calibratableParams,
             const Rcpp::List &constParams, int warmUp) : inputData(inputData),
                                                          parameters(calibratableParams, constParams),
                                                          simLength(inputData.nrows()),
                                                          warmUp(warmUp),
                                                          R(Rcpp::as<Rcpp::NumericVector>(inputData["R"])),
                                                          E0(Rcpp::as<Rcpp::NumericVector>(inputData["E0"])),
                                                          q1(Rcpp::as<Rcpp::NumericVector>(inputData["q1"])),
                                                          hf(Rcpp::as<Rcpp::NumericVector>(inputData["hf"])),
                                                          phi1(Rcpp::as<Rcpp::NumericVector>(inputData["phi1"])),
                                                          phi2(Rcpp::as<Rcpp::NumericVector>(inputData["phi2"]))

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
        Rcpp::Named("phi1") = bucket1._phi(),
        Rcpp::Named("phi2") = bucket2._phi(),
        Rcpp::Named("P") = bucket2.P(),
        // Wrechg
        Rcpp::Named("Wnet") = aquifer.Wnet(),
        // SW
        Rcpp::Named("Phi") = cn.Phi(),
        // S
        Rcpp::Named("S") = cn.S(),
        // Ia
        Rcpp::Named("Ia") = cn.Ia(),
        // CN
        Rcpp::Named("CN") = cn.CN(),
        // Qsurf
        Rcpp::Named("Q") = cn.Q(),
        // finf_SB1
        Rcpp::Named("F") = bucket1.F(),
        // finfla_SB1
        Rcpp::Named("FI") = bucket1.FI());
}

Rcpp::NumericVector Model::Q()
{
    return cn.Q();
}

Rcpp::DataFrame Model::q0()
{
    Rcpp::NumericVector t_temp = inputData["t"];
    t_temp = t_temp;
    return Rcpp::DataFrame::create(
        Rcpp::Named("t") = t_temp,
        Rcpp::Named("hf") = aquifer._hf(),
        Rcpp::Named("q0") = aquifer.q0(),
        Rcpp::Named("xn") = aquifer.xn(),
        Rcpp::Named("hn") = aquifer.hn(),
        Rcpp::Named("Wnet") = aquifer.Wnet(),
        Rcpp::Named("_Wnet") = aquifer._Wnet());
}

Rcpp::List Model::pars()
{
    return parameters.all_pars();
}
Rcpp::DataFrame Model::input()
{
    return Rcpp::DataFrame::create(
        Rcpp::Named("t") = inputData["t"],
        Rcpp::Named("R") = inputData["R"],
        Rcpp::Named("E0") = inputData["E0"],
        Rcpp::Named("q1") = inputData["q1"]);
}

// RCPP_EXPOSED_CLASS(Model);