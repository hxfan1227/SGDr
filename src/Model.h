#ifndef MODEL_H
#define MODEL_H
#include "Parameters.h"
class Model
{
private:
    Rcpp::DataFrame inputData;
    Parameters parameters;
    int simLength;
    Rcpp::NumericVector R;       // precipitation (mm)
    Rcpp::NumericVector E0;      // potential evapotranspiration (mm)
    Rcpp::NumericVector Pumping; // pumping rate (m3/day)
    Rcpp::NumericVector GWL;           // groundwater level (m)
    Rcpp::NumericVector H2O_SB1;       // head in soil bucket 1 (mm)
    Rcpp::NumericVector H2O_SB2;       // head in soil bucket 2 (mm)
    Rcpp::NumericVector SW;            // water in entire soil profile excluding wilting point water (mm)
    Rcpp::NumericVector S;             // retention parameter for a given day (mm)
    Rcpp::NumericVector Ia;            // interception (includes surface storage in rills) and infiltration prior to runoff (mm)
    Rcpp::NumericVector CN;            // curve number for a given day
    Rcpp::NumericVector Qsurf;         // surface runoff for a given day (mm)
    Rcpp::NumericVector MC_SB1;        // moisture content in SB1 (maximum = porosity, -)
    Rcpp::NumericVector DoS_SB1;       // degree of saturation in SB1 = MC/Porosity
    Rcpp::NumericVector finf_SB1;      // infiltration rate at time t (mm/hr)
    Rcpp::NumericVector finfla_SB1;    //
    Rcpp::NumericVector H2O1_SB1;      // water in bucket after day's infiltration (mm)
    Rcpp::NumericVector IntercH2O;     //
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
    Rcpp::NumericVector WrechgAve;
    Rcpp::NumericVector H2O2_AQ;
    Rcpp::NumericVector H2O3_AQ;
    Rcpp::NumericVector SGD1;
    Rcpp::NumericVector SGD2;
    Rcpp::NumericVector xn;
    Rcpp::NumericVector xn1;
    Rcpp::NumericVector xn2;
    Rcpp::NumericVector hn;
    Rcpp::NumericVector hn1;
    Rcpp::NumericVector hn2;
    Rcpp::NumericVector M1;
    Rcpp::NumericVector M2;
    Rcpp::NumericVector xT1;
    Rcpp::NumericVector xT2;
    Rcpp::NumericVector SGD;
    Rcpp::NumericVector SGDdrop;
    Rcpp::NumericVector PumpingDrop;
    Rcpp::NumericVector dxT;
    Rcpp::NumericVector xT3;
    Rcpp::NumericVector SWvol;
    Rcpp::NumericVector FWLdrop;
    void initializeVector();
    Rcpp::NumericVector moving_average(Rcpp::NumericVector x, int windowSize);

public:
    Model(Rcpp::DataFrame  inputData, 
          const Rcpp::List &calibratableParams,
          const Rcpp::List &constParams);
    void calc_recharge();
    void update_head(int i);
    void calc_cn_runoff(int i);
    void calc_bucket1_h2o(int i);
    void calc_bucket2_h2o(int i);
    void calc_aquifer_recharge(int i);
    void calc_sgd(int windowSize);
    void update_gwl(int i);
    void calc_sgds(int i);
    void update_aq_head(int i);
    void update_parameters(const Rcpp::List &newCalibratableParams, const Rcpp::List &newConstParams);
    Rcpp::List get_all_params_list();
    // getters
    Rcpp::NumericVector get_SW();
    Rcpp::NumericVector get_S();
    Rcpp::NumericVector get_Ia();
    Rcpp::NumericVector get_CN();
    Rcpp::NumericVector get_Qsurf();
    Rcpp::NumericVector get_MC_SB1();
    Rcpp::NumericVector get_DoS_SB1();
    Rcpp::NumericVector get_finf_SB1();
    Rcpp::NumericVector get_finfla_SB1();
    Rcpp::NumericVector get_H2O1_SB1();
    Rcpp::NumericVector get_IntercH2O();
    Rcpp::NumericVector get_CanopyH2O();
    Rcpp::NumericVector get_E0_Int();
    Rcpp::NumericVector get_YE();
    Rcpp::NumericVector get_E();
    Rcpp::NumericVector get_H2O2_SB1();
    Rcpp::NumericVector get_QSE();
    Rcpp::NumericVector get_SWexcess();
    Rcpp::NumericVector get_Kunsat();
    Rcpp::NumericVector get_Kunsat2();
    Rcpp::NumericVector get_TT();
    Rcpp::NumericVector get_TT2();
    Rcpp::NumericVector get_Wdperc();
    Rcpp::NumericVector get_Wperc();
    Rcpp::NumericVector get_Wperc2();
    Rcpp::NumericVector get_H2O3_SB1();
    Rcpp::NumericVector get_H2O3_SB2();
    Rcpp::NumericVector get_MC_SB2();
    Rcpp::NumericVector get_DoS_SB2();
    Rcpp::NumericVector get_H2O1_SB2();
    Rcpp::NumericVector get_YE2();
    Rcpp::NumericVector get_Edd();
    Rcpp::NumericVector get_H2O2_SB2();
    Rcpp::NumericVector get_Sw();
    Rcpp::NumericVector get_Sw2();
    Rcpp::NumericVector get_SWexcess2();
    Rcpp::NumericVector get_Wrechg();
    Rcpp::NumericVector get_H2O1_AQ();
    Rcpp::NumericVector get_WrechgAve();
    Rcpp::NumericVector get_H2O2_AQ();
    Rcpp::NumericVector get_H2O3_AQ();
    Rcpp::NumericVector get_SGD1();
    Rcpp::NumericVector get_SGD2();
    Rcpp::NumericVector get_xn();
    Rcpp::NumericVector get_xn1();
    Rcpp::NumericVector get_xn2();
    Rcpp::NumericVector get_hn();
    Rcpp::NumericVector get_hn1();
    Rcpp::NumericVector get_hn2();
    Rcpp::NumericVector get_M1();
    Rcpp::NumericVector get_M2();
    Rcpp::NumericVector get_xT1();
    Rcpp::NumericVector get_xT2();
    Rcpp::NumericVector get_SGD();
    Rcpp::NumericVector get_SGDdrop();
    Rcpp::NumericVector get_PumpingDrop();
    Rcpp::NumericVector get_dxT();
    Rcpp::NumericVector get_xT3();
    Rcpp::NumericVector get_SWvol();
    Rcpp::NumericVector get_FWLdrop();
    Rcpp::DataFrame get_recharge_output();
    Rcpp::DataFrame get_sgd_output();
};

#endif // MODEL_H