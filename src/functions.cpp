#include <Rcpp.h>
#include "Model.h"

//' Estimate the SFGD (q0 \[m2/d\]) .
//' @param inputData A data frame containing the input data for the model.
//' @param calibratableParams A list of calibratable parameters.
//' @param constParams A list of constant parameters.
//' @param nw A integer indicating the period you want to average the recharge.
//' @param warmup A integer indicating the number of days to warm up the model.
//' @description Estimate the SGD volume based on Strack (1976) analytical solution.
//' @details Estimate the SGD volume based on Strack (1976) analytical solution.
//' * The \code{inputData} should be a data.frame containing the required columns for the model (R, E0, q1, hf, phi1, phi2).
//'   * \code{R} is a numeric vector of the daily precipitation (mm).
//'   * \code{E0} is a numeric vector of the daily potential evapotranspiration (mm).
//'   * \code{q1} is a numeric vector of the daily pumping rate (m3/d).
//'   * \code{hf} is a numeric vector of the initial groundwater level (mm), Only the first value is used as the initial value for the SB1.
//'   * \code{phi1} is a numeric vector of the initial water level in the first soil bucket (mm). Only the first value is used as the initial value for the SB1.
//'   * \code{phi2} is a numeric vector of the inital water level in the second soil bucket (mm). Only the first value is used as the initial value for the SB2.
//' * \code{calibratableParams} should be a list of calibratable parameters. It's recommended to use \code{\link{json_to_paramter_list}} to create this list.
//' * \code{nw} A integer indicating the period you want to average the recharge.
//' * \code{warmup} A integer indicating the number of days to warm up the model.
//' @return A SGD_ESTIMATION_DF class.
//' @export
// [[Rcpp::export]]
Rcpp::List estimate_sgd(const Rcpp::DataFrame &inputData, const Rcpp::List &calibratableParams, const Rcpp::List &constParams, int nw = 120, int warmup = 1500)
{
    Model model(inputData, calibratableParams, constParams, warmup);
    model.run(nw, warmup);
    Rcpp::List output = Rcpp::List::create(Rcpp::Named("results") = model.q0(),
                                           Rcpp::Named("parameters") = model.pars(),
                                           Rcpp::Named("input") = model.input(),
                                           Rcpp::Named("runoff") = model.Q());
    output.attr("class") = Rcpp::CharacterVector::create("SGD_ESTIMATION_DF");
    return output;
}

//' Prepare the warm-up data
//' @rdname prepare_warm_up
//' @param data A data.frame containing the input data for the model.
//' @param length An integer indicating the length of the warm-up period (days).
//' @return A data.frame containing the input data for the model with the warm-up period.
//' @export
// [[Rcpp::export]]

Rcpp::DataFrame prepare_warm_up(const Rcpp::DataFrame &data, int length)
{
    if (length > data.nrows())
    {
        Rcpp::stop("n is greater than the number of rows in the DataFrame");
    }
    int total_rows = length + data.nrows();
    Rcpp::CharacterVector col_names = data.names();
    Rcpp::List new_cols(data.size());
    for (int i = 0; i < data.size(); i++)
    {
        // 取出当前列
        switch (TYPEOF(data[i]))
        {
        case INTSXP:
        {
            Rcpp::IntegerVector col = data[i];
            Rcpp::IntegerVector new_col(total_rows);
            std::copy(col.begin(), col.begin() + length, new_col.begin());
            std::copy(col.begin(), col.end(), new_col.begin() + length);
            new_cols[i] = new_col;
            break;
        }
        case REALSXP:
        {
            Rcpp::NumericVector col = data[i];
            Rcpp::NumericVector new_col(total_rows);
            std::copy(col.begin(), col.begin() + length, new_col.begin());
            std::copy(col.begin(), col.end(), new_col.begin() + length);
            new_cols[i] = new_col;
            break;
        }
        case STRSXP:
        {
            Rcpp::CharacterVector col = data[i];
            Rcpp::CharacterVector new_col(total_rows);
            std::copy(col.begin(), col.begin() + length, new_col.begin());
            std::copy(col.begin(), col.end(), new_col.begin() + length);
            new_cols[i] = new_col;
            break;
        }
        default:
            Rcpp::stop("Unsupported column type");
        }
    }
    Rcpp::DataFrame combined = Rcpp::DataFrame::create(new_cols, Rcpp::Named("stringsAsFactors") = false);
    combined.attr("names") = col_names;
    combined.attr("class") = "data.frame";
    combined.attr("row.names") = Rcpp::seq(1, total_rows);
    return combined;
}
