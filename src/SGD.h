#ifndef SGD_H_
#define SGD_H_
#include <Rcpp.h>
Rcpp::NumericVector moving_average(Rcpp::NumericVector x, int windowSize);
class Bucket {
private:
  Rcpp::List bucketParamters;

public:
  // Constructor that initializes the class with an Rcpp::List
  Bucket(Rcpp::List bucket) : bucketParamters(bucket) {}

  // Template method to get a parameter of various types
  template <typename T>
  T getParameter(const std::string& key) {
    if (bucketParamters.containsElementNamed(key.c_str())) {
      return Rcpp::as<T>(bucketParamters[key]);
    } else {
      throw std::runtime_error("Key not found in bucket");
    }
  }
};
// Expose the class to R with the specific methods instantiated for expected types
RCPP_MODULE(BucketModule) {
  Rcpp::class_<Bucket>("Bucket")
    .constructor<Rcpp::List>()
    .method("getParameterDouble", &Bucket::getParameter<double>)
    .method("getParameterInt", &Bucket::getParameter<int>)
    .method("getParameterString", &Bucket::getParameter<std::string>);
}
#endif // SGD_H_
