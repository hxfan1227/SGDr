# SGDr (development version)

* The unknown pumping rate is now available for calibration. (Still not a general purpose function, just for case study)

* The window size is now available for calibration.

* The warm up period is now available for the model, implemented in `C++`.

* The package is now using S3 method for `SGD_ESTIMATION_DF` class.

* `plot`, `print` and `summary` methods are available for the `SGD_ESTIMATION_DF` class.

* The package has been totally re-constructed using OOP. 
Now we have 2 classes which are `Model` and `Parameters`, respectively. 
We have achieved same results and efficiency when applying these changes.
However the extensibility is much larger. 

* The `Model` class and `Parameters` class is still under development.
