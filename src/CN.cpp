#include "SGD.h"
#include <Rcpp.h>
//' @export
//' @param sb1 A list of parameters
// [[Rcpp::export]]
double calc_cn_runoff (Rcpp::List sb1){
    Bucket sb1_bucket(sb1);
    return sb1_bucket.get_FCmm();
} 