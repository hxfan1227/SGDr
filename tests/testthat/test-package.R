test_that("all input data have been coerced to correct type", {
  expect_s3_class(json_to_parameter_list(test_path('testdata', 'preferred.json')), 'list')
  expect_s3_class(json_to_parameter_list(test_path('testdata', 'consts.json')), 'list')
})

test_that("calibratable parameters are correctly loaded", {
  x <- json_to_parameter_list(test_path('testdata', 'preferred.json'))
  expect_contains(names(x), c('bucket1', 'bucket2', 'curveNumber', 'aquifer'))
  expect_contains(names(x$bucket1), c('layer', 'z', 'rho_b', 'mc', 'rho_s', 'Ksat', 'n', 'Y'))
  expect_contains(names(x$bucket2), c('layer', 'z', 'rho_b', 'mc', 'rho_s', 'Ksat', 'n', 'Z'))
  expect_contains(names(x$curveNumber), c('CN2', 'PIa'))
  expect_contains(names(x$aquifer), c('delta', 'rho_s', 'rho_f', 'K', 'z0', 'Sy', 'xT', 'dxT'))
})

test_that("constant parameters are correctly loaded", {
  x <- json_to_parameter_list(test_path('testdata', 'consts.json'))
  expect_contains(names(x), c('x', 'W', 'Area', 'L', 'waterContent'))
  expect_type(x$waterContent, 'list')
})

test_that("estimate_sgd() checks parameters", {
  expect_error(estimate_sgd(1:2, 1:2, 1:2, 120, 1500))
})

test_that("estimate_sgd() returns correct type", {
  inputData <- read.csv(test_path('testdata', 'test_data.csv'))
  calibratableParamsList <- json_to_parameter_list(test_path('testdata', 'preferred.json'))
  constParamsList <- json_to_parameter_list(test_path('testdata', 'consts.json'))
  expect_s3_class(estimate_sgd(inputData, calibratableParamsList, constParamsList), 'SGD_ESTIMATION_DF')
})

test_that("print() warns if no base_date is available", {
  inputData <- read.csv(test_path('testdata', 'test_data.csv'))
  calibratableParamsList <- json_to_parameter_list(test_path('testdata', 'preferred.json'))
  constParamsList <- json_to_parameter_list(test_path('testdata', 'consts.json'))
  x <- estimate_sgd(inputData, calibratableParamsList, constParamsList)
  expect_warning(print(x), 'No base date available. Please provide the base date to calculate the simulation period. \n')
  expect_warning(summary(x), 'No base date available. Please provide the base date to calculate the simulation period. \n')
})

test_that("plot() returns correct object type", {
  inputData <- read.csv(test_path('testdata', 'test_data.csv'))
  obsData <- read.csv(test_path('testdata', 'obs_data.csv')) %>% 
    mutate(date = ymd(obs_date)) %>%
      select(date, wl) %>%
      distinct(date, .keep_all = T) %>%
      filter(wl <= 2)
  calibratableParamsList <- json_to_parameter_list(test_path('testdata', 'preferred.json'))
  constParamsList <- json_to_parameter_list(test_path('testdata', 'consts.json'))
  x <- estimate_sgd(inputData, calibratableParamsList, constParamsList)
  y <- plot(x, y = obsData, type = 'comp', base_date = 19670101)
  expect_s3_class(y, 'SGD_ESTIMATION_DF')
  expect_equal(x, y)
})

test_that("plot() throws error when no base_date is available", {
  inputData <- read.csv(test_path('testdata', 'test_data.csv'))
  obsData <- read.csv(test_path('testdata', 'obs_data.csv')) %>% 
    mutate(date = ymd(obs_date)) %>%
      select(date, wl) %>%
      distinct(date, .keep_all = T) %>%
      filter(wl <= 2)
  calibratableParamsList <- json_to_parameter_list(test_path('testdata', 'preferred.json'))
  constParamsList <- json_to_parameter_list(test_path('testdata', 'consts.json'))
  x <- estimate_sgd(inputData, calibratableParamsList, constParamsList)
  expect_error(plot(x, y = obsData, type = 'comp'))
  expect_error(plot(x, y = obsData, type = 'input'))
  expect_error(plot(x, y = obsData, type = 'pred'))
})

test_that("estimate_sgd_from_pars() returns same reults as estimate_sgd()", {
  inputData <- read.csv(test_path('testdata', 'test_data.csv'))
  calibratableParamsList <- json_to_parameter_list(test_path('testdata', 'preferred.json'))
  calibratableParams <- parameter_list_to_vector(calibratableParamsList)
  constParamsList <- json_to_parameter_list(test_path('testdata', 'consts.json'))
  parNames <- names(calibratableParams)[which(!(names(calibratableParams) %in% c('bucket1.layer', 'bucket2.layer')))]
  x <- estimate_sgd(inputData, calibratableParamsList, constParamsList, warmUp = 0)
  y <- estimate_sgd_from_pars(pars = calibratableParams[parNames],
                              input = inputData,
                              const_par_list = constParamsList,
                              parset = parNames,
                              default_pars = calibratableParams,
                              skeleton = calibratableParamsList, 
                              parnames = parNames, warm_up = 0)
  expect_equal(x, y)
})



