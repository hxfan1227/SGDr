#ifndef Parameters_H
#define Parameters_H
#include <Rcpp.h>
#include "bucket.h"
#include "cn.h"
#include "aquifer.h"
#include "constpar.h"

class Parameters
{
private:
    Bucket bucket1_;
    Bucket bucket2_;
    CurveNumber cn_;
    Aquifer aquifer_;
    ConstParameter constpar_;

public:
    Parameters(Rcpp::List cali_pars, Rcpp::List const_pars);
    Bucket &bucket1();
    Bucket &bucket2();
    CurveNumber &cn();
    Aquifer &aquifer();
    ConstParameter &constpar();
    Rcpp::List cali_pars();
    Rcpp::List const_pars();
    Rcpp::List all_pars();
    void update(const Rcpp::List &new_cali_pars, const Rcpp::List &new_const_pars);
};

#endif // Parameters_H