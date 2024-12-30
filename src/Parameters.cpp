#include "Parameters.h"

Parameters::Parameters(Rcpp::List calibratableParams, Rcpp::List constParams) : bucket1_(Rcpp::as<Rcpp::List>(calibratableParams["bucket1"]), constParams),
                                                                                bucket2_(Rcpp::as<Rcpp::List>(calibratableParams["bucket2"]), constParams),
                                                                                cn_(Rcpp::as<Rcpp::List>(calibratableParams["cn"]), bucket1_, bucket2_),
                                                                                aquifer_(Rcpp::as<Rcpp::List>(calibratableParams["aquifer"])),
                                                                                constpar_(constParams)
{
}

Bucket &Parameters::bucket1() { return bucket1_; }

Bucket &Parameters::bucket2() { return bucket2_; }

CurveNumber &Parameters::cn() { return cn_; }

Aquifer &Parameters::aquifer() { return aquifer_; }

ConstParameter &Parameters::constpar() { return constpar_; }

Rcpp::List Parameters::all_pars()
{
    return Rcpp::List::create(
        Rcpp::Named("bucket1") = bucket1_.all_pars(),
        Rcpp::Named("bucket2") = bucket2_.all_pars(),
        Rcpp::Named("cn") = cn_.all_pars(),
        Rcpp::Named("aquifer") = aquifer_.all_pars());
}

Rcpp::List Parameters::const_pars()
{
    return constpar_.const_pars();
}

Rcpp::List Parameters::cali_pars()
{
    return Rcpp::List::create(
        Rcpp::Named("bucket1") = bucket1_.cali_pars(),
        Rcpp::Named("bucket2") = bucket2_.cali_pars(),
        Rcpp::Named("cn") = cn_.cali_pars(),
        Rcpp::Named("aquifer") = aquifer_.cali_pars());
}

void Parameters::update(const Rcpp::List &newCalibratableParams, const Rcpp::List &newConstParams)
{
    bucket1_ = Bucket(Rcpp::as<Rcpp::List>(newCalibratableParams["bucket1"]), newConstParams);
    bucket2_ = Bucket(Rcpp::as<Rcpp::List>(newCalibratableParams["bucket2"]), newConstParams);
    cn_ = CurveNumber(Rcpp::as<Rcpp::List>(newCalibratableParams["cn"]), bucket1_, bucket2_);
    aquifer_ = Aquifer(Rcpp::as<Rcpp::List>(newCalibratableParams["aquifer"]));
    constpar_ = ConstParameter(newConstParams);
}

// RCPP_EXPOSED_CLASS(Parameters);
