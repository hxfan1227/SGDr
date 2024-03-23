#ifndef MODEL_H
#define MODEL_H
#include <Rcpp.h>
#include "Paramters.h"

class Model
{
private:
    Parameters parameters;
    const int simLength;
    const Rcpp::NumericVector R;          // precipitation (mm)
    const Rcpp::NumericVector E0;         // potential evapotranspiration (mm)
    Rcpp::NumericVector H2O_SB1;    // head in soil bucket 1 (mm)
    Rcpp::NumericVector H2O_SB2;    // head in soil bucket 2 (mm)
    Rcpp::NumericVector SW;         // water in entire soil profile excluding wilting point water (mm)
    Rcpp::NumericVector S;          // retention parameter for a given day (mm)
    Rcpp::NumericVector Ia;         // interception (includes surface storage in rills) and infiltration prior to runoff (mm)
    Rcpp::NumericVector CN;         // curve number for a given day
    Rcpp::NumericVector Qsurf;      // surface runoff for a given day (mm)
    Rcpp::NumericVector MC_SB1;     // moisture content in SB1 (maximum = porosity, -)
    Rcpp::NumericVector DoS_SB1;    // degree of saturation in SB1 = MC/Porosity
    Rcpp::NumericVector finf_SB1;   // infiltration rate at time t (mm/hr)
    Rcpp::NumericVector finfla_SB1; //
    Rcpp::NumericVector H2O1_SB1;   // water in bucket after day's infiltration (mm)
    Rcpp::NumericVector IntercH2O;  //
    Rcpp::NumericVector CanopyH2O;
    Rcpp::NumericVector E0_Int;   // Penman monteith methodology
    Rcpp::NumericVector YE;       //
    Rcpp::NumericVector E;        // amount of water removed from layer ly by evaporation (mm H2O)
    Rcpp::NumericVector H2O2_SB1; // water in bucket after day's infiltration and ET (mm)
    Rcpp::NumericVector QSE;
    Rcpp::NumericVector SWexcess;
    Rcpp::NumericVector Kunsat;
    Rcpp::NumericVector Kunsat2;
    Rcpp::NumericVector TT;
    Rcpp::NumericVector TT2;
    Rcpp::NumericVector Wdperc;
    Rcpp::NumericVector Wperc;
    Rcpp::NumericVector Wperc2;
    Rcpp::NumericVector H2O3_SB1; // water in bucket after infiltration, ET and percolation (mm)
    Rcpp::NumericVector H2O3_SB2; // water in bucket after percolation in, ET and percolation out (mm)
    Rcpp::NumericVector MC_SB2;   // moisture content in SB2 (maximum = porosity, -)
    Rcpp::NumericVector DoS_SB2;  // degree of saturation in SB2 = MC/Porosity
    Rcpp::NumericVector H2O1_SB2; // water in SB2 after percolation.
    Rcpp::NumericVector YE2;      // (100-Y)%E'' in SB2
    Rcpp::NumericVector Edd;      // E'' in SB2
    Rcpp::NumericVector H2O2_SB2; //  Water in SB2 after percolation.
    Rcpp::NumericVector Sw;
    Rcpp::NumericVector Sw2;
    Rcpp::NumericVector SWexcess2;
    Rcpp::NumericVector Wrechg;
    Rcpp::NumericVector H2O1_AQ;
    void initializeVector();

public:
    Model(const Rcpp::NumericVector& Time,
          const Rcpp::NumericVector& Precipitation,
          const Rcpp::NumericVector& Evaporation,
          Rcpp::NumericVector H2OinSB1,
          Rcpp::NumericVector H2OinSB2,
          const Rcpp::List& calibratableParams,
          const Rcpp::List& constParams);
          Rcpp::DataFrame calc_recharge();
    // Rcpp::NumericVector calc_sgd();
};

#endif