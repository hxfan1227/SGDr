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
  inputData <- read.csv(test_path('testdata', 'test_data.csv'))
  expect_snapshot_error(estimate_sgd(inputData, constParams = c(1, 2, 3), calibratableParams = c(1, 2, 3), warmUp = 1e5))
  expect_snapshot_error(estimate_sgd(inputData, constParams = c(1, 2, 3), calibratableParams = c(1, 2, 3), warmUp = 120))
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
  expect_snapshot_warning(print(x))
  expect_snapshot_warning(summary(x))
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
  expect_snapshot_error(plot(x, y = obsData, type = 'comp'))
  expect_snapshot_error(plot(x, y = obsData, type = 'input'))
  expect_snapshot_error(plot(x, y = obsData, type = 'pred'))
})

test_that("estimate_sgd_from_pars() returns same reults as estimate_sgd()", {
  inputData <- read.csv(test_path('testdata', 'test_data.csv'))
  calibratableParamsList <- json_to_parameter_list(test_path('testdata', 'preferred.json'))
  calibratableParams <- parameter_list_to_vector(calibratableParamsList)
  constParamsList <- json_to_parameter_list(test_path('testdata', 'consts.json'))
  parNames <- names(calibratableParams)[which(!(names(calibratableParams) %in% c('bucket1.layer', 'bucket2.layer')))]
  x <- estimate_sgd(inputData, calibratableParamsList, constParamsList, warmUp = 1500, windowSize = 180)
  y <- estimate_sgd_from_pars(pars = c(calibratableParams[parNames], nw = 180),
                              parset = parNames,
                              parnames = c(parNames, 'nw'),
                              inputDf = inputData,
                              calibratableParams = calibratableParamsList,
                              constParams = constParamsList,
                              skeleton = calibratableParamsList, 
                              warmUp = 1500)
  expect_equal(x, y)
})

test_that('estimate_sgd_from_pars() checks parameter names', {
  inputData <- read.csv(test_path('testdata', 'test_data.csv'))
  calibratableParamsList <- json_to_parameter_list(test_path('testdata', 'preferred.json'))
  calibratableParams <- parameter_list_to_vector(calibratableParamsList)
  constParamsList <- json_to_parameter_list(test_path('testdata', 'consts.json'))
  parNames <- names(calibratableParams)[which(!(names(calibratableParams) %in% c('bucket1.layer', 'bucket2.layer')))]
  expect_snapshot_error(estimate_sgd_from_pars(pars = c(calibratableParams[parNames], nw = 120),
                                      parset = parNames,
                                      parnames = parNames,
                                      inputDf = inputData,
                                      calibratableParams = calibratableParamsList,
                                      constParams = constParamsList,
                                      skeleton = calibratableParamsList, 
                                      warmUp = 1500))
  expect_snapshot_error(estimate_sgd_from_pars(pars = c(calibratableParams[parNames], nw = 120),
                                               parset = c(parNames, 'nw'),
                                               parnames = c(parNames, 'nw'),
                                               inputDf = inputData,
                                               calibratableParams = calibratableParamsList,
                                               constParams = constParamsList,
                                               skeleton = calibratableParamsList, 
                                               warmUp = 1500))
  expect_snapshot_error(estimate_sgd_from_pars(pars = c(calibratableParams[parNames], nw = 120),
                                               parset = c(parNames, 'nw'),
                                               parnames = c(parNames, 'wn'),
                                               inputDf = inputData,
                                               calibratableParams = calibratableParamsList,
                                               constParams = constParamsList,
                                               skeleton = calibratableParamsList, 
                                               warmUp = 1500))
})

test_that('change_unknown_pumping() works', {
  inputData <- read.csv(test_path('testdata', 'test_data.csv')) %>% 
    setDT()
  calibratableParamsList <- json_to_parameter_list(test_path('testdata', 'preferred.json'))
  calibratableParams <- parameter_list_to_vector(calibratableParamsList)
  constParamsList <- json_to_parameter_list(test_path('testdata', 'consts.json'))
  parNames <- names(calibratableParams)[which(!(names(calibratableParams) %in% c('bucket1.layer', 'bucket2.layer')))]
  yearly_pumping_df <- read.csv(test_path('testdata', 'yearly_pumping_data.csv')) %>% 
    setDT()
  monthly_pumping_df <- read.csv(test_path('testdata', 'monthly_pumping_data.csv')) %>% 
    setDT()
  pumping_data <- data.table(date = seq(ymd(19670101), by = '1 day', length.out = length(inputData$t)))
  pumping_data[, ':='(year = lubridate::year(date), month = lubridate::month(date))]
  pumping_data <- pumping_data[monthly_pumping_df,  percent := i.percent, on = .(month)]
  pumping_data <- pumping_data[yearly_pumping_df, on = .(year)]
  pumping_data[, daily_pumping := pumping * percent * 1000 / lubridate::days_in_month(date)]
  pumping_data[, t := as.integer(date - ymd(19670101)) + 1]
  inputData[, Pumping := NULL]
  preferred_pumping <- 15
  new_input <- change_unknown_pumping(x = preferred_pumping, input_df = inputData, pumping_df = pumping_data, yearly_df = yearly_pumping_df) 
  tempfile <- tempfile(tmpdir = test_path('testdata'), fileext = '.csv')
  write.csv(new_input, tempfile, row.names = F)
  input_from_file <- read.csv(tempfile)
  x <- estimate_sgd(input_from_file, calibratableParamsList, constParamsList, warmUp = 1500, windowSize = 180)
  y <- estimate_sgd_from_pars(pars = c(calibratableParams[parNames], 
                                       nw = 180, 
                                       pumping = preferred_pumping),
                              parset = parNames,
                              parnames = c(parNames, 'nw', 'pumping'),
                              yearlyPumping = yearly_pumping_df,
                              pumpingDf = pumping_data,
                              skeleton = calibratableParamsList, 
                              inputDf = input_from_file,
                              calibratableParams = calibratableParamsList,
                              constParams = constParamsList,
                              warmUp = 1500)
  unlink(tempfile)
  expect_equal(x, y)
})



