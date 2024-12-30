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
        // [timeseries] amount of water held in the soil profile at the begining of the day [mm H2O]
        Rcpp::Named("phi1") = bucket1.phi(),
        // [timeseries] amount of water held in the soil profile at the beginning of the day [mm H2O]
        Rcpp::Named("phi2") = bucket2.phi(),
        // [timeseries] actual percolation from soil bucket i to the underlying bucket limited to field capacity and available space of the underlying bucket [mm H2O]
        Rcpp::Named("P") = bucket2.P(),
        // [timeseries] recharge entering the aquifer bucket [mm H2O]
        Rcpp::Named("Wnet") = aquifer.Wnet(),
        // [timeseries] amount of water held in the soil (summing the soil water amounts in soil buckets 1 and 2) excluding the water held in the soil profile at wilting point [mm H2O]
        Rcpp::Named("Phi") = cn.Phi(),
        // [timeseries] retention parameter [mm H2O]
        Rcpp::Named("S") = cn.S(),
        // [timeseries] initial abstraction that includes surface storage, interception and infiltration prior to runoff [mm H2O]
        Rcpp::Named("Ia") = cn.Ia(),
        // [timeseries] curve number [-]
        Rcpp::Named("CN") = cn.CN(),
        // [timeseries] runoff [mm H2O]
        Rcpp::Named("Q") = cn.Q(),
        // [timeseries] infiltration without considering the pre-runoff infiltration [mm H2O]
        Rcpp::Named("F") = bucket1.F(),
        // [timeseries] pre-runoff infiltration [mm H2O]
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
        Rcpp::Named("hf") = aquifer.hf(),
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