#include "Parameters.h"

Parameters::Parameters(Rcpp::List calibratableParams, Rcpp::List constParams) : bucket1(Rcpp::as<Rcpp::List>(calibratableParams["bucket1"]), constParams),
                                                                                bucket2(Rcpp::as<Rcpp::List>(calibratableParams["bucket2"]), constParams),
                                                                                curveNumber(Rcpp::as<Rcpp::List>(calibratableParams["curveNumber"]), bucket1, bucket2),
                                                                                aquifer(Rcpp::as<Rcpp::List>(calibratableParams["aquifer"])),
                                                                                constParameter(constParams)
{
}

Bucket Parameters::get_bucket1() { return bucket1; }

Bucket Parameters::get_bucket2() { return bucket2; }

CurveNumber Parameters::get_curveNumber() { return curveNumber; }

Aquifer Parameters::get_aquifer() { return aquifer; }

ConstParameter Parameters::get_constParameter() { return constParameter; }

Rcpp::List Parameters::get_all_params_list()
{
    return Rcpp::List::create(
        Rcpp::Named("bucket1") = bucket1.get_all_params_list(),
        Rcpp::Named("bucket2") = bucket2.get_all_params_list(),
        Rcpp::Named("curveNumber") = curveNumber.get_all_params_list(),
        Rcpp::Named("aquifer") = aquifer.get_all_params_list());
}

Rcpp::List Parameters::get_const_params_list()
{
    return constParameter.get_all_params_list();
}

Rcpp::List Parameters::get_calibratable_params_list()
{
    return Rcpp::List::create(
        Rcpp::Named("bucket1") = bucket1.get_calibratable_params_list(),
        Rcpp::Named("bucket2") = bucket2.get_calibratable_params_list(),
        Rcpp::Named("curveNumber") = curveNumber.get_calibratable_params_list(),
        Rcpp::Named("aquifer") = aquifer.get_calibratable_params_list());
}

void Parameters::update(const Rcpp::List& newCalibratableParams, const Rcpp::List& newConstParams)
{
    bucket1 = Bucket(Rcpp::as<Rcpp::List>(newCalibratableParams["bucket1"]), newConstParams);
    bucket2 = Bucket(Rcpp::as<Rcpp::List>(newCalibratableParams["bucket2"]), newConstParams);
    curveNumber = CurveNumber(Rcpp::as<Rcpp::List>(newCalibratableParams["curveNumber"]), bucket1, bucket2);
    aquifer = Aquifer(Rcpp::as<Rcpp::List>(newCalibratableParams["aquifer"]));
    constParameter = ConstParameter(newConstParams);
}

// RCPP_EXPOSED_CLASS(Parameters);

RCPP_MODULE(ParametersModule)
{
    Rcpp::class_<Parameters>("Parameters")
        .constructor<Rcpp::List, Rcpp::List>()
        .method("get_all_params_list", &Parameters::get_all_params_list)
        .method("get_const_params_list", &Parameters::get_const_params_list)
        .method("get_calibratable_params_list", &Parameters::get_calibratable_params_list)
        .method("update", &Parameters::update);
}

double Bucket::interp1(Rcpp::NumericVector x, Rcpp::NumericVector y, double xi)
{
    Rcpp::Environment pracmaEnv = Rcpp::Environment::namespace_env("pracma");
    Rcpp::Function interp1_r = pracmaEnv["interp1"];
    return Rcpp::as<double>(interp1_r(x, y, xi));
}

Bucket::Bucket(Rcpp::List bucketParams, Rcpp::List Params)
{
    // ructor implementation remains the same
    if (!bucketParams.containsElementNamed("layer"))
    {
        Rcpp::stop("Bucket parameter layer not found");
    }
    layer = Rcpp::as<int>(bucketParams["layer"]);
    if (!bucketParams.containsElementNamed("z"))
    {
        Rcpp::stop("Bucket parameter z not found");
    }
    z = Rcpp::as<double>(bucketParams["z"]);
    if (!bucketParams.containsElementNamed("rho_b"))
    {
        Rcpp::stop("Bucket parameter rho_b not found");
    }
    rho_b = Rcpp::as<double>(bucketParams["rho_b"]);
    if (!bucketParams.containsElementNamed("mc"))
    {
        Rcpp::stop("Bucket parameter mc not found");
    }
    mc = Rcpp::as<double>(bucketParams["mc"]);
    if (!bucketParams.containsElementNamed("rho_s"))
    {
        Rcpp::stop("Bucket parameter rho_s not found");
    }
    rho_s = Rcpp::as<double>(bucketParams["rho_s"]);
    if (!bucketParams.containsElementNamed("Ksat"))
    {
        Rcpp::stop("Bucket parameter Ksat not found");
    }
    Ksat = Rcpp::as<double>(bucketParams["Ksat"]);
    if (!bucketParams.containsElementNamed("n"))
    {
        Rcpp::stop("Bucket parameter n not found");
    }
    n = Rcpp::as<double>(bucketParams["n"]);
    if (!bucketParams.containsElementNamed("Y"))
    {
        if (layer == 1)
        {
            Rcpp::stop("Bucket parameter Y not found");
        }
        else
        {
            Y = 0;
        }
    }
    else
    {
        Y = Rcpp::as<double>(bucketParams["Y"]);
    }
    if (!bucketParams.containsElementNamed("Z"))
    {
        if (layer > 1)
        {
            Rcpp::stop("Bucket parameter Z not found");
        }
        else
        {
            Z = 0;
        }
    }
    else
    {
        Z = Rcpp::as<double>(bucketParams["Z"]);
    }
    calculate_bucket_params(Params);
}

int Bucket::get_layer() { return layer; }
double Bucket::get_z() { return z; }
double Bucket::get_rho_b() { return rho_b; }
double Bucket::get_mc() { return mc; }
double Bucket::get_rho_s() { return rho_s; }
double Bucket::get_SAT() { return SAT; }
double Bucket::get_n() { return n; }
double Bucket::get_WPmm() { return WPmm; }
double Bucket::get_WP() { return WP; }
double Bucket::get_FC() { return FC; }
double Bucket::get_FCmm() { return FCmm; }
double Bucket::get_Swres() { return Swres; }
double Bucket::get_Ksat() { return Ksat; }
double Bucket::get_phi_soil() { return phi_soil; }
double Bucket::get_Y() { return Y; }
double Bucket::get_Z() { return Z; }
double Bucket::get_Ksatmd() { return Ksatmd; }
double Bucket::get_TTperc() { return TTperc; }
double Bucket::get_AWClyr() { return AWClyr; }
Rcpp::List Bucket::get_all_params_list()
{
    return Rcpp::List::create(Rcpp::Named("layer") = layer,
                              Rcpp::Named("z") = z,
                              Rcpp::Named("rho_b") = rho_b,
                              Rcpp::Named("mc") = mc,
                              Rcpp::Named("rho_s") = rho_s,
                              Rcpp::Named("SAT") = SAT,
                              Rcpp::Named("n") = n,
                              Rcpp::Named("WPmm") = WPmm,
                              Rcpp::Named("WP") = WP,
                              Rcpp::Named("FC") = FC,
                              Rcpp::Named("FCmm") = FCmm,
                              Rcpp::Named("Swres") = Swres,
                              Rcpp::Named("Ksat") = Ksat,
                              Rcpp::Named("phi_soil") = phi_soil,
                              Rcpp::Named("Y") = Y,
                              Rcpp::Named("Z") = Z,
                              Rcpp::Named("Ksatmd") = Ksatmd,
                              Rcpp::Named("TTperc") = TTperc,
                              Rcpp::Named("AWClyr") = AWClyr);
}

void Bucket::calculate_bucket_params(Rcpp::List Params)
{
    if (!Params.containsElementNamed("waterContent"))
    {
        Rcpp::stop(" parameter waterContent not found");
    }
    Rcpp::List waterContent = Rcpp::as<Rcpp::List>(Params["waterContent"]);
    phi_soil = 1 - rho_b / rho_s;
    SAT = phi_soil * z;
    WP = 0.4 * ((mc * rho_b) / 100);
    WPmm = WP * z;
    Swres = WPmm / SAT;
    FC = interp1(waterContent["ClayContent"], waterContent["FC"], mc);
    FCmm = FC * z;
    AWClyr = FC - WP;
    Ksatmd = Ksat / 1000 * 24;
    TTperc = (SAT - FCmm) / Ksat;
}

Rcpp::List Bucket::get_calibratable_params_list()
{
    return Rcpp::List::create(Rcpp::Named("layer") = layer,
                              Rcpp::Named("z") = z,
                              Rcpp::Named("rho_b") = rho_b,
                              Rcpp::Named("mc") = mc,
                              Rcpp::Named("rho_s") = rho_s,
                              Rcpp::Named("Ksat") = Ksat,
                              Rcpp::Named("n") = n,
                              Rcpp::Named("Y") = Y,
                              Rcpp::Named("Z") = Z);
}

Aquifer::Aquifer(const Rcpp::List aquiferParams)
{
    if (!aquiferParams.containsElementNamed("delta"))
    {
        Rcpp::stop("Aquifer parameter delta not found");
    }
    delta = Rcpp::as<double>(aquiferParams["delta"]);
    if (!aquiferParams.containsElementNamed("rho_s"))
    {
        Rcpp::stop("Aquifer parameter rho_s not found");
    }
    rho_s = Rcpp::as<double>(aquiferParams["rho_s"]);
    if (!aquiferParams.containsElementNamed("rho_f"))
    {
        Rcpp::stop("Aquifer parameter rho_f not found");
    }
    rho_f = Rcpp::as<double>(aquiferParams["rho_f"]);
    if (!aquiferParams.containsElementNamed("K"))
    {
        Rcpp::stop("Aquifer parameter K not found");
    }
    K = Rcpp::as<double>(aquiferParams["K"]);
    if (!aquiferParams.containsElementNamed("z0"))
    {
        Rcpp::stop("Aquifer parameter z0 not found");
    }
    z0 = Rcpp::as<double>(aquiferParams["z0"]);
    if (!aquiferParams.containsElementNamed("Sy"))
    {
        Rcpp::stop("Aquifer parameter Sy not found");
    }
    Sy = Rcpp::as<double>(aquiferParams["Sy"]);
    if (!aquiferParams.containsElementNamed("xT"))
    {
        Rcpp::stop("Aquifer parameter xT not found");
    }
    xT = Rcpp::as<double>(aquiferParams["xT"]);
    if (!aquiferParams.containsElementNamed("dxT"))
    {
        Rcpp::stop("Aquifer parameter dxT not found");
    }
    dxT = Rcpp::as<double>(aquiferParams["dxT"]);
    a = (rho_s - rho_f) / rho_f;
    he = z0 * rho_s / rho_f - z0;
}

double Aquifer::get_delta()
{
    return {delta};
}

double Aquifer::get_rho_s()
{
    return {rho_s};
}

double Aquifer::get_rho_f()
{
    return {rho_f};
}

double Aquifer::get_K()
{
    return {K};
}

double Aquifer::get_z0()
{
    return {z0};
}

double Aquifer::get_Sy()
{
    return {Sy};
}

double Aquifer::get_a()
{
    return {a};
}

double Aquifer::get_he()
{
    return {he};
}

double Aquifer::get_xT()
{
    return {xT};
}

double Aquifer::get_dxT()
{
    return {dxT};
}

Rcpp::List Aquifer::get_all_params_list()
{
    return Rcpp::List::create(Rcpp::Named("delta") = delta,
                              Rcpp::Named("rho_s") = rho_s,
                              Rcpp::Named("rho_f") = rho_f,
                              Rcpp::Named("K") = K,
                              Rcpp::Named("z0") = z0,
                              Rcpp::Named("Sy") = Sy,
                              Rcpp::Named("a") = a,
                              Rcpp::Named("he") = he,
                              Rcpp::Named("xT") = xT,
                              Rcpp::Named("dxT") = dxT);
}

Rcpp::List Aquifer::get_calibratable_params_list()
{
    return Rcpp::List::create(Rcpp::Named("delta") = delta,
                              Rcpp::Named("K") = K,
                              Rcpp::Named("z0") = z0,
                              Rcpp::Named("Sy") = Sy,
                              Rcpp::Named("xT") = xT,
                              Rcpp::Named("dxT") = dxT);
}

ConstParameter::ConstParameter(Rcpp::List constParams)
{
    if (!constParams.containsElementNamed("x"))
    {
        Rcpp::stop("Const parameter x not found");
    }
    x = Rcpp::as<double>(constParams["x"]);
    if (!constParams.containsElementNamed("W"))
    {
        Rcpp::stop("Const parameter W not found");
    }
    W = Rcpp::as<double>(constParams["W"]);
    if (!constParams.containsElementNamed("Area"))
    {
        Rcpp::stop("Const parameter Area not found");
    }
    Area = Rcpp::as<double>(constParams["Area"]);
    if (!constParams.containsElementNamed("L"))
    {
        Rcpp::stop("Const parameter L not found");
    }
    L = Rcpp::as<double>(constParams["L"]);
    if (!constParams.containsElementNamed("waterContent"))
    {
        Rcpp::stop("Const parameter waterContent not found");
    }
    waterContent = Rcpp::as<Rcpp::List>(constParams["waterContent"]);
}

double ConstParameter::get_x()
{
    return {x};
}

double ConstParameter::get_W()
{
    return {W};
}

double ConstParameter::get_Area()
{
    return {Area};
}

double ConstParameter::get_L()
{
    return {L};
}

Rcpp::List ConstParameter::get_all_params_list()
{
    return Rcpp::List::create(Rcpp::Named("x") = x,
                              Rcpp::Named("W") = W,
                              Rcpp::Named("Area") = Area,
                              Rcpp::Named("L") = L,
                              Rcpp::Named("waterContent") = waterContent);
}

CurveNumber::CurveNumber(const Rcpp::List cnParams, Bucket &SB1, Bucket &SB2)
{
    if (!cnParams.containsElementNamed("CN2"))
    {
        Rcpp::stop("Curve number parameter CN2 not found");
    }
    CN2 = Rcpp::as<double>(cnParams["CN2"]);
    if (!cnParams.containsElementNamed("PIa"))
    {
        Rcpp::stop("Curve number parameter PIa not found");
    }
    PIa = Rcpp::as<double>(cnParams["PIa"]);
    CN1 = CN2 - 20 * (100 - CN2) / (100 - CN2 + exp(2.533 - 0.0636 * (100 - CN2)));
    CN3 = CN2 * exp(0.00673 * (100 - CN2));
    Smax = 25.4 * ((1000 / CN1) - 10);
    S3 = 25.4 * ((1000 / CN3) - 10);
    w2 = (log(((SB1.get_FCmm() + SB2.get_FCmm()) / (1 - S3 / Smax)) - (SB1.get_FCmm() + SB2.get_FCmm())) - log((SB1.get_SAT() + SB2.get_SAT()) / (1 - 2.54 / Smax) - (SB1.get_SAT() + SB2.get_SAT()))) / ((SB1.get_SAT() + SB2.get_SAT()) - (SB1.get_FCmm() + SB2.get_FCmm()));
    w1 = log((SB1.get_FCmm() + SB2.get_FCmm()) / (1 - S3 / Smax) - (SB1.get_FCmm() + SB2.get_FCmm())) + w2 * (SB1.get_FCmm() + SB2.get_FCmm());
}

double CurveNumber::get_CN2()
{
    return CN2;
}

double CurveNumber::get_PIa()
{
    return PIa;
}

double CurveNumber::get_CN1()
{
    return CN1;
}

double CurveNumber::get_CN3()
{
    return CN3;
}

double CurveNumber::get_Smax()
{
    return Smax;
}

double CurveNumber::get_S3()
{
    return S3;
}

double CurveNumber::get_w1()
{
    return w1;
}

double CurveNumber::get_w2()
{
    return w2;
}

Rcpp::List CurveNumber::get_all_params_list()
{

    return Rcpp::List::create(Rcpp::Named("CN2") = CN2,
                              Rcpp::Named("PIa") = PIa,
                              Rcpp::Named("CN1") = CN1,
                              Rcpp::Named("CN3") = CN3,
                              Rcpp::Named("Smax") = Smax,
                              Rcpp::Named("S3") = S3,
                              Rcpp::Named("w1") = w1,
                              Rcpp::Named("w2") = w2);
}

Rcpp::List CurveNumber::get_calibratable_params_list()
{
    return Rcpp::List::create(Rcpp::Named("CN2") = CN2,
                              Rcpp::Named("PIa") = PIa);
}
