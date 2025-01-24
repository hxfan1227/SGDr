# SGDr 0.1.1

## Major changes

* GUI for SFGD calculator is availble now.

* `SGDr_example()` makes it easy to access example files bundled with `SGDr`.

## Minor changes

* Fix typos in the documents

# SGDr 0.1.0

## Major changes

* First public release of `SGDr`

* Restructured in Rcpp, with major modifications made to each class.

* 3 more classes have been added (i.e., `CurveNumber` for Curve Number method, `Bucket` for the soil layer and `Aquifer` for the aquifer layer)

## Minor changes

* slight modification made to the `plot` method.

# SGDr (development version)

* The unknown pumping rate is now available for calibration. (not a general purpose function, just for case study)

* The window size is now available for calibration.

* The warm up period is now available for the model, implemented in `C++`.


* `plot`, `print` and `summary` methods are available for the `SGD_ESTIMATION_DF` class.

* The package is now using `S3` method for `SGD_ESTIMATION_DF` class.

* The package has been completely restructured using Object-Oriented Programming (OOP). 
It now comprises two primary classes, namely `Model` and `Parameters`. 
Despite these structural modifications, the results and efficiency remain consistent with previous implementations. 
However, the changes have significantly enhanced the system's extensibility.
