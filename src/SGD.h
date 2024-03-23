#ifndef SGD_H_
#define SGD_H_
#include <Rcpp.h>

Rcpp::NumericVector moving_average(Rcpp::NumericVector x, int windowSize);


class Bucket {
private:
  int layer; // soil layer number in bucket (-). In this model we only have 2 layers (1 and 2).
  double WP; // wilting point in bucket (mm)
  double z; // soil depth in bucket (mm)
  double FCmm; // water at field capacity in bucket (mm)
  double SAT; // water in soil when saturated in bucket (mm H2O)
  double Swres; // vol H20 per vol voids in bucket (-)
  double Ksat; // saturated hydraulic conductivity in bucket (mm/hr)
  double n; // van Genuchten n in bucket (-)
  double phi_soil; // soil porosity in bucket (mg/m3)
  double Y; // percent of evaporation from bucket (%)
  double Z; // percent of FC to SAT (in SB2) when perc stops (%)

public:
  // Constructor that initializes the class with an Rcpp::List
  Bucket(const Rcpp::List bucketParams) {
    if (!bucketParams.containsElementNamed("layer")) {
      Rcpp::stop("Bucket parameter layer not found");
    }
    if (!bucketParams.containsElementNamed("WP")) {
      Rcpp::stop("Bucket parameter WP not found");
    }
    WP = Rcpp::as<double>(bucketParams["WP"]);
    if (!bucketParams.containsElementNamed("z")) {
      Rcpp::stop("Bucket parameter z not found");
    }
    z = Rcpp::as<double>(bucketParams["z"]);
    if (!bucketParams.containsElementNamed("FCmm")) {
      Rcpp::stop("Bucket parameter FCmm not found");
    }
    FCmm = Rcpp::as<double>(bucketParams["FCmm"]);
    if (!bucketParams.containsElementNamed("SAT")) {
      Rcpp::stop("Bucket parameter SAT not found");
    }
    SAT = Rcpp::as<double>(bucketParams["SAT"]);
    if (!bucketParams.containsElementNamed("Swres")) {
      Rcpp::stop("Bucket parameter Swres not found");
    }
    Swres = Rcpp::as<double>(bucketParams["Swres"]);
    if (!bucketParams.containsElementNamed("Ksat")) {
      Rcpp::stop("Bucket parameter Ksat not found");
    }
    Ksat = Rcpp::as<double>(bucketParams["Ksat"]);
    if (!bucketParams.containsElementNamed("n")) {
      Rcpp::stop("Bucket parameter n not found");
    }
    n = Rcpp::as<double>(bucketParams["n"]);
    if (!bucketParams.containsElementNamed("phi_soil")) {
      Rcpp::stop("Bucket parameter phi_soil not found");
    }
    phi_soil = Rcpp::as<double>(bucketParams["phi_soil"]);
    if (!bucketParams.containsElementNamed("Y")) {
      if (layer == 1) {
        Rcpp::stop("Bucket parameter Y not found");
      } else {
        Y = 0;
      }      
    } else {
        Y = Rcpp::as<double>(bucketParams["Y"]);  
    } 
    if (!bucketParams.containsElementNamed("Z")) {
      if (layer == 1) {
        Rcpp::stop("Bucket parameter Z not found");
      } else {
        Z = 0;
      }      
    } else {
        Z = Rcpp::as<double>(bucketParams["Z"]);  
    }      
  }

  // Getter methods for the class
  double get_WP() { return WP; }
  double get_z() { return z; }
  double get_FCmm() { return FCmm; }
  double get_SAT() { return SAT; }
  double get_Swres() { return Swres; }
  double get_Ksat() { return Ksat; }
  double get_n() { return n; }
  double get_phi_soil() { return phi_soil; }
  double get_Y() { return Y; }
  double get_Z() { return Z; }
};
// Expose the class to R with the specific methods instantiated for expected types
RCPP_MODULE(BucketModule) {
  Rcpp::class_<Bucket>("Bucket")
    .constructor<Rcpp::List>()
    .method("get_WP", &Bucket::get_WP)
    .method("get_z", &Bucket::get_z)
    .method("get_FCmm", &Bucket::get_FCmm)
    .method("get_SAT", &Bucket::get_SAT)
    .method("get_Swres", &Bucket::get_Swres)
    .method("get_Ksat", &Bucket::get_Ksat)
    .method("get_n", &Bucket::get_n)
    .method("get_phi_soil", &Bucket::get_phi_soil)
    .method("get_Y", &Bucket::get_Y)
    .method("get_Z", &Bucket::get_Z);
}
#endif // SGD_H_
