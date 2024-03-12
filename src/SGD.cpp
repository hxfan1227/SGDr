#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;

//' 
//' Return a numeric vector of moving average with leading NAs filled with first average value
//' 
//' @param x A numeric vector with no NAs
//' @param windowSize A integer of window size
//' @return A numeric vector of moving average with leading NAs filled with first average value
//' @export
// [[Rcpp::export]]

NumericVector moving_average(NumericVector x, int windowSize)
{
  int n = x.size();
  NumericVector avg(n); // Initialize the result vector
  
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



//' 
//' Calculate the surface water recharge using 1D model based on SWAT procedure
//' 
//' @param t A integer vector of simulation day
//' @param R A numeric vector of daily precipitation (mm)
//' @param E0 A numeric vector of daily potential evaporation calculated with Penman monteith method (mm)
//' @param H2O_SB1 A numeric vector of depth of water in SB1 (initial condition)
//' @param H2O_SB2 A numeric vector of depth of water in SB2 (initial condition)
//' @param params A list for parameters of the model
//' @param calibration If true return only the recharge value (for computational efficiency)
//' @return A data.frame with daily surface water recharge
//' @export
// [[Rcpp::export]]

DataFrame cal_recharge(NumericVector t, // Simulation day 
                       NumericVector R, // daily precipitation (mm)
                       NumericVector E0,// daily potential evaporation calculated with Penman monteith method (mm)
                       NumericVector H2O_SB1, // depth of water in SB1 (mm)
                       NumericVector H2O_SB2, // depth of water in SB2 (mm)
                       List params, // a list for parameters of the model
                       bool calibration // If true return only the recharge value (for computational efficiency)
                       ) 
{
  List bucket1 = as<List>(params["bucket1"]); // a list of parameters in SB1
  List bucket2 = as<List>(params["bucket2"]); // a list of parameters in SB2
  List cn = as<List>(params["cn"]); // a list of parameters in curve number approach
  List aq = as<List>(params["aq"]); // a list of parameters in aquifer
  double bucket1_WP1 = as<double>(bucket1["WP1"]); // wilting point in SB1 (mm)
  double bucket1_z1 = as<double>(bucket1["z1"]); // soil depth in SB1 (mm)
  double bucket2_WP2 = as<double>(bucket2["WP2"]); // wilting point in SB2 (mm)
  double bucket2_z2 = as<double>(bucket2["z2"]); //soil depth in SB2 (mm)
  double cn_Smax = as<double>(cn["Smax"]); // maximum retention parameter (mm H20)
  double cn_w1 = as<double>(cn["w1"]); // first shape coefficient (-)
  double cn_w2 = as<double>(cn["w2"]); // second shape coefficient (-)
  double bucket1_phi_soil = as<double>(bucket1["phi_soil"]); // soil porosity in SB1 (mg/m3)
  double bucket2_phi_soil = as<double>(bucket2["phi_soil"]); // soil porosity in SB2 (mg/m3)
  double cn_PIa = as<double>(cn["PIa"]); // proportion of Ia that is pre-runoff infiltration (-)
  double bucket1_Y = as<double>(bucket1["Y"]); // percent of evaporation from SB1 (%)
  double bucket1_FCmm = as<double>(bucket1["FCmm"]); // water at field capacity in SB1 (mm)
  double bucket1_SAT = as<double>(bucket1["SAT"]); // water in soil when saturated in SB1 (mm H2O)
  double bucket1_Swres = as<double>(bucket1["Swres"]); // vol H20 per vol voids in SB1 (-)
  double bucket2_Swres = as<double>(bucket2["Swres"]); // vol H20 per vol voids in SB2 (-)
  double bucket1_n = as<double>(bucket1["n"]); // van Genuchten n in SB1 (-)
  double bucket2_n = as<double>(bucket2["n"]); // van Genuchten n in SB2 (-)
  double bucket1_Ksat = as<double>(bucket1["Ksat"]); // saturated hydraulic conductivity in SB1 (mm/hr)
  double bucket2_Ksat = as<double>(bucket2["Ksat"]); // saturated hydraulic conductivity in SB2 (mm/hr)
  double bucket2_FCmm = as<double>(bucket2["FCmm"]);// water at field capacity in SB2 (mm)
  double bucket2_SAT = as<double>(bucket2["SAT"]); // water in soil when saturated in SB2 (mm H2O)
  double bucket1_Z = as<double>(bucket1["Z"]); // percent of FC to SAT (in SB2) when perc stops (%)
  double aq_delta = as<double>(aq["delta"]); // delay time (d)
  
  NumericVector SW(t.size()); // water in entire soil profile excluding wilting point water (mm)
  NumericVector S(t.size()); // retention parameter for a given day (mm)
  NumericVector Ia(t.size()); // interception (includes surface storage in rills) and infiltration prior to runoff (mm)
  NumericVector CN(t.size()); // curve number for a given day
  NumericVector Qsurf(t.size()); // surface runoff for a given day (mm)
  NumericVector MC_SB1(t.size()); // moisture content in SB1 (maximum = porosity, -)
  NumericVector DoS_SB1(t.size()); // degree of saturation in SB1 = MC/Porosity 
  NumericVector finf_SB1(t.size()); // infiltration rate at time t (mm/hr)
  NumericVector finfla_SB1(t.size());  // 
  NumericVector H2O1_SB1(t.size()); // water in bucket after day's infiltration (mm)
  NumericVector IntercH2O(t.size()); // 
  NumericVector CanopyH2O(t.size());
  NumericVector E0_Int(t.size()); // Penman monteith methodology 
  NumericVector YE(t.size()); //
  NumericVector E(t.size()); // amount of water removed from layer ly by evaporation (mm H2O)
  NumericVector H2O2_SB1(t.size()); // water in bucket after day's infiltration and ET (mm)
  NumericVector QSE(t.size());
  NumericVector SWexcess(t.size());
  NumericVector Kunsat(t.size());
  NumericVector Kunsat2(t.size());
  NumericVector TT(t.size());
  NumericVector TT2(t.size());
  NumericVector Wdperc(t.size());
  NumericVector Wperc(t.size());
  NumericVector Wperc2(t.size());
  NumericVector H2O3_SB1(t.size()); // water in bucket after infiltration, ET and percolation (mm)
  NumericVector H2O3_SB2(t.size()); // water in bucket after percolation in, ET and percolation out (mm)
  NumericVector MC_SB2(t.size()); // moisture content in SB2 (maximum = porosity, -)
  NumericVector DoS_SB2(t.size()); // degree of saturation in SB2 = MC/Porosity 
  NumericVector H2O1_SB2(t.size()); // water in SB2 after percolation.
  NumericVector YE2(t.size()); // (100-Y)%E'' in SB2
  NumericVector Edd(t.size()); // E'' in SB2
  NumericVector H2O2_SB2(t.size()); //  Water in SB2 after percolation.
  NumericVector Sw(t.size());
  NumericVector Sw2(t.size());
  NumericVector SWexcess2(t.size());
  NumericVector Wrechg(t.size());
  NumericVector H2O1_AQ(t.size());
  
  for (int i = 0; i < t.size(); ++i)
  {
    if (i > 0)
    {
      H2O_SB1[i] = H2O3_SB1[i - 1];
      H2O_SB2[i] = H2O3_SB2[i - 1];
    }
    // curve number
    SW[i] = H2O_SB1[i] + H2O_SB2[i] - (bucket1_WP1 * bucket1_z1 + bucket2_WP2 * bucket2_z2);
    S[i] = cn_Smax * (1 - SW[i] / (SW[i] + exp(cn_w1 - cn_w2 * SW[i])));
    Ia[i] = S[i] * 0.2;
    CN[i] = 25400 / (S[i] + 254);
    Qsurf[i] = (R[i] < Ia[i]) ? 0 : pow((R[i] - Ia[i]), 2) / (R[i] + 0.8 * S[i]);
    // bucket1
    MC_SB1[i] = H2O_SB1[i] / bucket1_z1;
    DoS_SB1[i] = MC_SB1[i] / bucket1_phi_soil;
    finf_SB1[i] = (R[i] < Ia[i]) ? 0 : R[i] - Qsurf[i] - Ia[i];
    finfla_SB1[i] = R[i] > Ia[i] ? cn_PIa * Ia[i] : (R[i] > (Ia[i] * (1 - cn_PIa)) ? R[i] - (1 - cn_PIa) * Ia[i] : 0);
    H2O1_SB1[i] = H2O_SB1[i] + finf_SB1[i] + finfla_SB1[i];
    IntercH2O[i] = R[i] > ((1 - cn_PIa) * Ia[i]) ? (1 - cn_PIa) * Ia[i] : R[i];
    CanopyH2O[i] = i == 0 ? 0.0 : std::max(CanopyH2O[i - 1] + IntercH2O[i] - E0[i], 0.0);
    E0_Int[i] = std::max(0.0, E0[i] - IntercH2O[i] - CanopyH2O[i]);
    YE[i] = E0_Int[i] * bucket1_Y / 100;
    E[i] = std::min(H2O1_SB1[i] < bucket1_FCmm ? YE[i] * exp(2.5 * (H2O1_SB1[i] - bucket1_FCmm) / (bucket1_FCmm - bucket1_WP1 * bucket1_z1)) : YE[i], 0.8 * (H2O1_SB1[i] - bucket1_WP1 * bucket1_z1));
    H2O2_SB1[i] = std::min(std::max(0.0, H2O1_SB1[i] - E[i]), bucket1_SAT);
    QSE[i] = ((H2O1_SB1[i] - E[i]) > bucket1_SAT) ? H2O1_SB1[i] - E[i] - bucket1_SAT : 0;
    SWexcess[i] = std::max(0.0, H2O2_SB1[i] - bucket1_FCmm);
    Sw[i] = std::min(1.0, (H2O2_SB1[i] / bucket1_SAT - bucket1_Swres) / (1 - bucket1_Swres));
    Kunsat[i] = std::max(0.00000001,
                         pow(Sw[i], 0.5) * pow((1 - pow(1 - pow(Sw[i], bucket1_n / (bucket1_n - 1)), (bucket1_n - 1) / bucket1_n)), 2) * bucket1_Ksat);
    TT[i] = std::max((bucket1_SAT - bucket1_FCmm) / Kunsat[i], 72.0);
    Wdperc[i] = SWexcess[i] * (1 - exp(-24 / TT[i]));
    Wperc[i] = H2O_SB2[i] >= bucket2_FCmm + bucket1_Z / 100 * (bucket2_SAT - bucket2_FCmm) ? 0.0 : (Wdperc[i] <= (bucket2_SAT - H2O_SB2[i]) ? Wdperc[i] : (bucket2_SAT - H2O_SB2[i]));
    H2O3_SB1[i] = H2O2_SB1[i] - Wperc[i];
    // bucket2
    MC_SB2[i] = H2O_SB2[i] / bucket2_z2;
    DoS_SB2[i] = MC_SB2[i] / bucket2_phi_soil;
    H2O1_SB2[i] = H2O_SB2[i] + Wperc[i];
    YE2[i] = E0_Int[i] * (100 - bucket1_Y) / 100;
    Edd[i] = std::min((H2O1_SB2[i] < bucket2_FCmm ? YE2[i] * exp(2.5 * (H2O1_SB2[i] - bucket2_FCmm) / (bucket2_FCmm - bucket2_WP2 * bucket2_z2)) : YE2[i]),
                      0.8 * (H2O1_SB2[i] - bucket2_WP2 * bucket2_z2));
    H2O2_SB2[i] = H2O1_SB2[i] - Edd[i];
    SWexcess2[i] = std::max(0.0, H2O1_SB2[i] - bucket2_FCmm);
    Sw2[i] = std::min(1.0, (H2O2_SB2[i] / bucket2_SAT - bucket2_Swres) / (1.0 - bucket2_Swres));
    Kunsat2[i] = std::max(0.00000001,
                          pow(Sw2[i], 0.5) *
                            pow(1 -
                            pow(1 -
                            pow(Sw2[i], bucket2_n / (bucket2_n - 1)),
                            (bucket2_n - 1) / bucket2_n),
                            2) *
                              bucket2_Ksat);
    TT2[i] = (bucket2_SAT - bucket2_FCmm) / Kunsat2[i];
    Wperc2[i] = SWexcess2[i] * (1 - exp(-24 / TT2[i]));
    H2O3_SB2[i] = H2O2_SB2[i] - Wperc2[i];
    // aquifer
    Wrechg[i] = i == 0 ? (1 - exp(-1 / aq_delta)) * Wperc2[i] : (1 - exp(-1 / aq_delta)) * Wperc2[i] + exp(-1 / aq_delta) * Wrechg[i - 1];
  }
  
  // List results = List::create(Named("SW") = SW,
  //                             Named("S") = S,
  //                             Named("Ia") = Ia,
  //                             Named("CN") = CN,
  //                             Named("Qsurf") = Qsurf,
  //                             Named("MC_SB1") = MC_SB1,
  //                             Named("DoS_SB1") = DoS_SB1,
  //                             Named("finf_SB1") = finf_SB1,
  //                             Named("finfla_SB1") = finfla_SB1,
  //                             Named("H2O3_SB1") = H2O3_SB1);
  
  DataFrame results = DataFrame::create(
    Named("H2O3_SB1") = H2O3_SB1,
    Named("H2O3_SB2") = H2O3_SB2,
    Named("Wperc2") = Wperc2,
    Named("Wrechg") = Wrechg
  // Named("SW") = SW,
  // Named("S") = S,
  // Named("Ia") = Ia,
  // Named("CN") = CN,
  // Named("Qsurf") = Qsurf,
  // Named("finf_SB1") = finf_SB1,
  // Named("finfla_SB1") = finfla_SB1
  );
  if (calibration)
  {
    return DataFrame::create(Named("Wrechg") = Wrechg);
  }
  
  return results;
}

//' 
//' Calculate the Submarine Groundwater Discharge using an analytical solution for sharp-interface steady state
//' 
//' @param GWL A numeric vector of daily groundwater level (mm)
//' @param Wrechg A numeric vector of daily surface water recharge (mm) calculated from `cal_recharge()`.
//' @param WrechgAve A numeric vector of moving average of `Wrechg`
//' @param Pumping A numeric vector of daily pumping rate
//' @param H2O_SB2 A numeric vector of depth of water in SB2 (initial condition)
//' @param params A list for parameters of the model
//' @param calibration If true return only the water level value (for computational efficiency)
//' @return A data.frame with daily SGD
//' @export
// [[Rcpp::export]]
DataFrame cal_sgd(NumericVector GWL,
                  NumericVector Wrechg,
                  NumericVector WrechgAve,
                  NumericVector Pumping,
                  List params,
                  bool calibration)
{
  NumericVector H2O1_AQ(GWL.size(), NA_REAL);
  NumericVector H2O2_AQ(GWL.size(), NA_REAL);
  NumericVector SGD1(GWL.size(), NA_REAL);
  NumericVector SGD2(GWL.size(), NA_REAL);
  NumericVector xn1(GWL.size(), NA_REAL);
  NumericVector xn2(GWL.size(), NA_REAL);
  NumericVector hn1(GWL.size(), NA_REAL);
  NumericVector hn2(GWL.size(), NA_REAL);
  NumericVector M1(GWL.size(), NA_REAL);
  NumericVector M2(GWL.size(), NA_REAL);
  NumericVector xT1(GWL.size(), NA_REAL);
  NumericVector xT2(GWL.size(), NA_REAL);
  NumericVector SGD(GWL.size(), NA_REAL);
  NumericVector SGDdrop(GWL.size(), NA_REAL);
  NumericVector PumpingDrop(GWL.size(), NA_REAL);
  
  List aq = as<List>(params["aq"]);
  
  double aq_Sy = as<double>(aq["Sy"]);
  double aq_K = as<double>(aq["K"]);
  double aq_x = as<double>(aq["x"]);
  double aq_rho_s = as<double>(aq["rho_s"]);
  double aq_rho_f = as<double>(aq["rho_f"]);
  double aq_W = as<double>(aq["W"]);
  double aq_z0 = as<double>(aq["z0"]);
  double aq_a = as<double>(aq["a"]);
  double aq_Area = as<double>(aq["Area"]);
  
  for (int i = 0; i < GWL.size(); ++i)
  {
    if (i > 0)
    {
      GWL[i] = H2O2_AQ[i - 1];
    }
    
    H2O1_AQ[i] = GWL[i] + Wrechg[i] / 1000 / aq_Sy;
    if (WrechgAve[i] > 1E-30)
    {
      SGD1[i] = ((aq_K / 2 / aq_x * aq_rho_s / (aq_rho_s - aq_rho_f) * std::pow(GWL[i], 2)) + (WrechgAve[i] / 1000) * aq_x / 2) * aq_W;
      SGD2[i] = (aq_K * ((std::pow((GWL[i] + aq_z0), 2) - aq_rho_s / aq_rho_f * std::pow(aq_z0, 2))) + (WrechgAve[i] / 1000) * std::pow(aq_x, 2)) / 2 / aq_x * aq_W;
      xn1[i] = SGD1[i] / (WrechgAve[i] / 1000) / aq_W;
      xn2[i] = SGD2[i] / (WrechgAve[i] / 1000) / aq_W;
      hn1[i] = sqrt(std::pow(SGD1[i] / aq_W, 2) / (WrechgAve[i] / 1000) / aq_K + aq_rho_s / aq_rho_f * aq_z0 * aq_z0) - aq_z0;
      hn2[i] = sqrt(std::pow(SGD2[i] / aq_W, 2) / (WrechgAve[i] / 1000) / aq_K + aq_rho_s / aq_rho_f * aq_z0 * aq_z0) - aq_z0;
      M1[i] = aq_K * aq_a * (1 + aq_a) * aq_z0 * aq_z0 / (WrechgAve[i] / 1000) / (xn1[i] * xn1[i]);
      M2[i] = aq_K * aq_a * (1 + aq_a) * aq_z0 * aq_z0 / (WrechgAve[i] / 1000) / (xn2[i] * xn2[i]);
      if (M1[i] < 1)
      {
        xT1[i] = (SGD1[i] / aq_W) / (WrechgAve[i] / 1000) - sqrt(pow((SGD1[i] / aq_W) / (WrechgAve[i] / 1000), 2) - aq_K * aq_a * (1 + aq_a) * (aq_z0 * aq_z0) / (WrechgAve[i] / 1000));
      }
      if (M2[i] < 1)
      {
        xT2[i] = (SGD2[i] / aq_W) / (WrechgAve[i] / 1000) - sqrt(pow((SGD2[i] / aq_W) / (WrechgAve[i] / 1000), 2) - aq_K * aq_a * (1 + aq_a) * (aq_z0 * aq_z0) / (WrechgAve[i] / 1000));
      }
    }
    else
    {
      SGD1[i] = (aq_K / 2 / aq_x * aq_rho_s / (aq_rho_s - aq_rho_f) * std::pow(GWL[i], 2)) * aq_W;
      SGD2[i] = (aq_K * ((std::pow((GWL[i] + aq_z0), 2) - aq_rho_s / aq_rho_f * std::pow(aq_z0, 2))) / 2 / aq_x * aq_W);
      xn1[i] = 1000000000;
      xn2[i] = 1000000000;
      hn1[i] = 1000;
      hn2[i] = 1000;
      M1[i] = 0;
      M2[i] = 0;
      if (M1[i] < 1)
      {
        xT1[i] = aq_K * aq_a * (1 + aq_a) * aq_z0 * aq_z0 / 2 / (SGD1[i] / aq_W);
      }
      if (M2[i] < 1)
      {
        xT2[i] = aq_K * aq_a * (1 + aq_a) * aq_z0 * aq_z0 / 2 / (SGD2[i] / aq_W);
      }
    }
    
    if (xT2[i] < 0)
    {
      SGD[i] = SGD1[i];
    }
    else
    {
      if (M1[i] >= 1)
      {
        SGD[i] = SGD1[i];
      }
      else
      {
        if (aq_W <= xT2[i])
        {
          SGD[i] = SGD1[i];
        }
        else
        {
          if (aq_x > xT2[i])
          {
            SGD[i] = SGD2[i];
          }
          else
          {
            SGD[i] = SGD1[i];
          }
        }
      }
    }
    
    SGDdrop[i] = SGD[i] / aq_Area * 1000 / aq_Sy;
    PumpingDrop[i] = Pumping[i] / aq_Area * 1000 / aq_Sy;
    H2O2_AQ[i] = H2O1_AQ[i] - SGDdrop[i] / 1000 - PumpingDrop[i] / 1000;
    
    // SGD1[i] = (WrechgAve[i] > 1E-30) ?
    // ((aq_K / 2 / aq_x * aq_rho_s / (aq_rho_s - aq_rho_f) * std::pow(GWL[i], 2)) + (WrechgAve[i] / 1000) * aq_x / 2) * aq_W :
    //   (aq_K / 2 / aq_x * aq_rho_s / (aq_rho_s - aq_rho_f) * std::pow(GWL[i], 2)) * aq_W;
    // SGD2[i] = (WrechgAve[i] > 1E-30) ?
    // ((aq_K * ((std::pow((GWL[i] + aq_z0), 2) - aq_rho_s / aq_rho_f * std::pow(aq_z0, 2)) + (WrechgAve[i] / 1000) * std::pow(aq_x, 2))) / 2 / aq_x * aq_W) :
    //   (aq_K * ((std::pow((GWL[i] + aq_z0), 2) - aq_rho_s / aq_rho_f * std::pow(aq_z0, 2))) / 2 / aq_x * aq_W);
    // xn1[i] = (WrechgAve[i] > 1E-30) ? SGD1[i] / (WrechgAve[i] / 1000) / aq_W : 1000000000;
    // xn2[i] = (WrechgAve[i] > 1E-30) ? SGD2[i] / (WrechgAve[i] / 1000) / aq_W : 1000000000;
    // hn1[i] = (WrechgAve[i] > 1E-30) ?
    // sqrt(std::pow(SGD1[i] / aq_W, 2) / (WrechgAve[i] / 1000) / aq_K + aq_rho_s / aq_rho_f * aq_z0 * aq_z0) - aq_z0 : 1000;
    // hn2[i] = (WrechgAve[i] > 1E-30) ?
    // sqrt(std::pow(SGD2[i] / aq_W, 2) / (WrechgAve[i] / 1000) / aq_K + aq_rho_s / aq_rho_f * aq_z0 * aq_z0) - aq_z0 : 1000;
  }
  DataFrame result = DataFrame::create(Named("H2O1_AQ") = H2O1_AQ,
                                       Named("SGD1") = SGD1,
                                       Named("SGD2") = SGD2,
                                       Named("xn1") = xn1,
                                       Named("xn2") = xn2,
                                       Named("hn1") = hn1,
                                       Named("hn2") = hn2,
                                       Named("M1") = M1,
                                       Named("M2") = M2,
                                       Named("xT1") = xT1,
                                       Named("xT2") = xT2,
                                       Named("SGD") = SGD,
                                       Named("SGDdrop") = SGDdrop,
                                       Named("PumpingDrop") = PumpingDrop,
                                       Named("H2O2_AQ") = H2O2_AQ
                                         
  );
  if (calibration)
  {
    return DataFrame::create(Named("pred") = H2O2_AQ);
  }
  return result;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

// /*** R
// test_data <- read_csv('./Data/test_data.csv')
// test_data %>%
//   mutate(cal_recharge(t, R, E0, H2O_SB1, H2O_SB2, full_params, calibration = T)) %>%
//   mutate(Wrechg180 = moving_average(Wrechg, 180)) %>%
//   mutate(cal_sgd(GWL, Wrechg, Wrechg180, Pumping, full_params, calibration = T))
// */
