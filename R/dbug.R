# library(DEoptim)
# library(tidyverse)
# library(SGDr)
# library(scales)
# 
# DEoptim2 <- function (fn, lower, upper, control = DEoptim.control(), ..., 
#                       fnMap = NULL) 
# {
#   if (length(lower) != length(upper)) 
#     stop("'lower' and 'upper' are not of same length")
#   if (!is.vector(lower)) 
#     lower <- as.vector(lower)
#   if (!is.vector(upper)) 
#     upper <- as.vector(upper)
#   if (any(lower > upper)) 
#     stop("'lower' > 'upper'")
#   if (any(lower == "Inf")) 
#     warning("you set a component of 'lower' to 'Inf'. May imply 'NaN' results", 
#             immediate. = TRUE)
#   if (any(lower == "-Inf")) 
#     warning("you set a component of 'lower' to '-Inf'. May imply 'NaN' results", 
#             immediate. = TRUE)
#   if (any(upper == "Inf")) 
#     warning("you set a component of 'upper' to 'Inf'. May imply 'NaN' results", 
#             immediate. = TRUE)
#   if (any(upper == "-Inf")) 
#     warning("you set a component of 'upper' to '-Inf'. May imply 'NaN' results", 
#             immediate. = TRUE)
#   if (!is.null(names(lower))) 
#     nam <- names(lower)
#   else if (!is.null(names(upper)) && is.null(names(lower))) 
#     nam <- names(upper)
#   else nam <- paste("par", 1:length(lower), sep = "")
#   ctrl <- do.call(DEoptim.control, as.list(control))
#   ctrl$npar <- length(lower)
#   if (is.na(ctrl$NP)) 
#     ctrl$NP <- 10 * length(lower)
#   if (ctrl$NP < 4) {
#     warning("'NP' < 4; set to default value 10*length(lower)\n", 
#             immediate. = TRUE)
#     ctrl$NP <- 10 * length(lower)
#   }
#   if (ctrl$NP < 10 * length(lower)) 
#     warning("For many problems it is best to set 'NP' (in 'control') to be at least ten times the length of the parameter vector. \n", 
#             immediate. = TRUE)
#   if (!is.null(ctrl$initialpop)) {
#     ctrl$specinitialpop <- TRUE
#     if (!identical(as.numeric(dim(ctrl$initialpop)), as.numeric(c(ctrl$NP, 
#                                                                   ctrl$npar)))) 
#       stop("Initial population is not a matrix with dim. NP x length(upper).")
#   }
#   else {
#     ctrl$specinitialpop <- FALSE
#     ctrl$initialpop <- 0
#   }
#   ctrl$trace <- as.numeric(ctrl$trace)
#   ctrl$specinitialpop <- as.numeric(ctrl$specinitialpop)
#   ctrl$initialpop <- as.numeric(ctrl$initialpop)
#   if (!is.null(ctrl$cluster)) {
#     if (!inherits(ctrl$cluster, "cluster")) 
#       stop("cluster is not a 'cluster' class object")
#     parallel::clusterExport(ctrl$cluster, ctrl$parVar)
#     if (!is.null(ctrl$parallelArgs)) 
#       fnPop <- function(`*params`, ...) {
#         parallel::parApply(cl = ctrl$cluster, X = `*params`, 
#                            MARGIN = 1, FUN = fn, ctrl$parallelArgs, ...)
#       }
#     else fnPop <- function(`*params`, ...) {
#       parallel::parApply(cl = ctrl$cluster, X = `*params`, 
#                          MARGIN = 1, FUN = fn, ...)
#     }
#   }
#   else if (ctrl$parallelType == "foreach") {
#     use.foreach <- "foreach" %in% installed.packages()
#     if (!use.foreach) 
#       stop("foreach package not available but parallelType set to 'foreach'")
#     if (!foreach::getDoParRegistered()) {
#       foreach::registerDoSEQ()
#     }
#     args <- ctrl$foreachArgs
#     fnPop <- function(`*params`, ...) {
#       my_chunksize <- ceiling(NROW(`*params`)/foreach::getDoParWorkers())
#       my_iter <- iterators::iter(`*params`, by = "row", 
#                                  chunksize = my_chunksize)
#       args$i <- my_iter
#       if (is.null(args$.combine)) 
#         args$.combine <- c
#       if (!is.null(args$.export)) 
#         args$.export = c(args$.export, "fn")
#       else args$.export = "fn"
#       if (is.null(args$.errorhandling)) 
#         args$.errorhandling = c("stop", "remove", "pass")
#       if (is.null(args$.verbose)) 
#         args$.verbose = FALSE
#       if (is.null(args$.inorder)) 
#         args$.inorder = TRUE
#       if (is.null(args$.multicombine)) 
#         args$.multicombine = FALSE
#       foreach::"%dopar%"(do.call(foreach::foreach, args), 
#                          apply(X = i, MARGIN = 1, FUN = fn, ...))
#     }
#   }
#   else if (ctrl$parallelType == "parallel") {
#     if (!requireNamespace("parallelly", quietly = TRUE)) {
#       stop("the parallelly package is required for parallelType = 'parallel'\n", 
#            "please install it via install.packages('parallelly')")
#     }
#     cl <- parallel::makeCluster(parallelly::availableCores())
#     packFn <- function(packages) {
#       for (i in packages) library(i, character.only = TRUE)
#     }
#     parallel::clusterCall(cl, packFn, ctrl$packages)
#     if (is.null(ctrl$parVar)) 
#       ctrl$parVar <- ls()
#     parallel::clusterExport(cl = cl, varlist = ctrl$parVar, 
#                             envir = environment())
#     fnPop <- function(`*params`, ...) {
#       parallel::parApply(cl = cl, X = `*params`, MARGIN = 1, 
#                          FUN = fn, ...)
#     }
#   }
#   else {
#     fnPop <- function(`*params`, ...) {
#       apply(X = `*params`, MARGIN = 1, FUN = fn, ...)
#     }
#   }
#   fnMapC <- NULL
#   if (!is.null(fnMap)) {
#     fnMapC <- function(`*params`, ...) {
#       mappedPop <- t(apply(X = `*params`, MARGIN = 1, FUN = fnMap))
#       if (all(dim(mappedPop) != dim(`*params`))) 
#         stop("mapping function did not return an object with ", 
#              "dim NP x length(upper).")
#       dups <- duplicated(mappedPop)
#       np <- NCOL(mappedPop)
#       tries <- 0
#       while (tries < 5 && any(dups)) {
#         nd <- sum(dups)
#         newPop <- matrix(runif(nd * np), ncol = np)
#         newPop <- rep(lower, each = nd) + newPop * rep(upper - 
#                                                          lower, each = nd)
#         mappedPop[dups, ] <- t(apply(newPop, MARGIN = 1, 
#                                      FUN = fnMap))
#         dups <- duplicated(mappedPop)
#         tries <- tries + 1
#       }
#       if (tries == 5) 
#         warning("Could not remove ", sum(dups), " duplicates from the mapped ", 
#                 "population in 5 tries. Evaluating population with duplicates.", 
#                 call. = FALSE, immediate. = TRUE)
#       storage.mode(mappedPop) <- "double"
#       mappedPop
#     }
#   }
#   outC <- .Call("DEoptimC", lower, upper, fnPop, ctrl, new.env(), 
#                 fnMapC, PACKAGE = "DEoptim")
#   # if (ctrl$parallelType == "parallel") 
#   #   parallel::stopCluster(cl)
#   if (length(outC$storepop) > 0) {
#     nstorepop <- floor((outC$iter - ctrl$storepopfrom)/ctrl$storepopfreq)
#     storepop <- list()
#     cnt <- 1
#     for (i in 1:nstorepop) {
#       idx <- cnt:((cnt - 1) + (ctrl$NP * ctrl$npar))
#       storepop[[i]] <- matrix(outC$storepop[idx], nrow = ctrl$NP, 
#                               ncol = ctrl$npar, byrow = TRUE)
#       cnt <- cnt + (ctrl$NP * ctrl$npar)
#       dimnames(storepop[[i]]) <- list(1:ctrl$NP, nam)
#     }
#   }
#   else {
#     storepop = NULL
#   }
#   names(outC$bestmem) <- nam
#   iter <- max(1, as.numeric(outC$iter))
#   names(lower) <- names(upper) <- nam
#   bestmemit <- matrix(outC$bestmemit[1:(iter * ctrl$npar)], 
#                       nrow = iter, ncol = ctrl$npar, byrow = TRUE)
#   dimnames(bestmemit) <- list(1:iter, nam)
#   storepop <- as.list(storepop)
#   outR <- list(optim = list(bestmem = outC$bestmem, bestval = outC$bestval, 
#                             nfeval = outC$nfeval, iter = outC$iter), member = list(lower = lower, 
#                                                                                    upper = upper, bestmemit = bestmemit, bestvalit = outC$bestvalit, 
#                                                                                    pop = t(outC$pop), storepop = storepop))
#   attr(outR, "class") <- "DEoptim"
#   return(outR)
# }
# 
# lowerB <- json_to_paramter_list('./tests/testthat/lower_limits.json') %>% unlist()
# upperB <- json_to_paramter_list('./tests/testthat/upper_limits.json') %>% unlist()
# ws <- json_to_paramter_list('./tests/testthat/evidence_weights.json') %>% unlist()
# initial_params <- json_to_paramter_list('./tests/testthat/tanya.json') %>% unlist()
# preferred_params <- json_to_paramter_list('./tests/testthat/preferred.json') %>% unlist()
# const_params <- json_to_paramter_list('./tests/testthat/consts.json')
# params_skeleton <- json_to_paramter_list('./tests/testthat/tanya.json')
# 
# weighted.rmse <- function(actual, predicted, weight){
#   sqrt(sum((predicted-actual)^2*weight)/sum(weight))
# }
# 
# test_data <- read.csv('C:/RProjects/DailyWork/CoffinBay/Scripts/SGD/Data/test_data.csv')
# obs_dt <- read.csv('C:/RProjects/DailyWork/CoffinBay/Scripts/SGD/Data/obs_data.csv') %>% 
#   select(obs_date, wl) %>% 
#   rename('date' = 'obs_date') %>% 
#   pivot_longer(-date, values_to = 'obs') %>% 
#   mutate(name = factor(name, levels = c('wl', 'hn', 'xn', 'Wrechg', 'WrechgAve'), 
#                        labels = c('WL~(m~AHD)', 
#                                   'h[n]~(m~AHD)',
#                                   'x[n]~(m)', 
#                                   'w[recharge]~(mm)', 
#                                   'w[recharge[ave]]~(mm)')),
#          date = dmy(date))
# 
# ids <- as.integer(obs_dt$date - ymd('1967-01-01') + 1)
# obs_wl <- obs_dt$obs
# 
# export_model <- function() {
#   sgd_model <- create_sgd_model(test_data, relist(initial_params), const_params)
#   sgd_model$calc_recharge()
#   sgd_model$calc_sgd(180)
#   sgd_model
# }
# 
# 
# obj_fun <- function(params, 
#                     skeleton = params_skeleton, 
#                     weights = ws, 
#                     obs = obs_wl, 
#                     xn = 3000, hn = 23, 
#                     default_values = preferred_params, 
#                     lower_bound = lowerB, 
#                     upper_bound = upperB, ...) {
#   model <- create_sgd_model(test_data, relist(params, skeleton = skeleton), const_params)
#   model$calc_recharge()
#   model$calc_sgd(180)
#   if (model$get_all_params_list()$aquifer$Sy > model$get_all_params_list()$bucket2$phi_soil) {
#     return(Inf)
#   }
#   h_n <- model$get_hn()
#   h_n[h_n <= hn] <- hn
#   x_n <- model$get_xn()
#   x_n[x_n <= x_n] <- xn
#   
#   hydroGOF::nrmse(model$get_H2O3_AQ()[ids], obs, norm = 'maxmin') * 0.5 + 
#     hydroGOF::rmse(h_n[ids], rep(hn, length(obs))) / hn * 0.5 + 
#     hydroGOF::rmse(x_n[ids], rep(xn, length(obs))) / xn * 5 + 
#     weighted.rmse(default_values/(upper_bound - lower_bound + 1E-30), params/(upper_bound - lower_bound + 1E-30), weights) * 0.5
# }
# 
# obj_fun(initial_params)
# 
# model_diagnosis <- function(model, obs_data) {
#   x <- data.frame(
#     t = model$get_inputData()$t,
#     wl = model$get_H2O3_AQ(), 
#     hn = model$get_hn(), 
#     xn = model$get_xn(),
#     Wrechg = model$get_Wrechg(),
#     WrechgAve = model$get_WrechgAve()
#   )
#   plot_df <- x %>% pivot_longer(cols = -t, values_to = 'pred') %>% 
#     mutate(date = ymd(19670101) - 1 + t) %>% 
#     mutate(name = factor(name, levels = c('wl', 'hn', 'xn', 'Wrechg', 'WrechgAve'), 
#                          labels = c('WL~(m~AHD)', 
#                                     'h[n]~(m~AHD)',
#                                     'x[n]~(m)', 
#                                     'w[recharge]~(mm)', 
#                                     'w[recharge[ave]]~(mm)'))) %>% 
#     left_join(obs_data, by = c('date', 'name'))
#   ggplot(plot_df, aes(x = date, y = pred)) +
#     geom_line(linewidth = 0.8) +
#     geom_point(data = drop_na(plot_df), aes(y = obs, color = 'Observed Water Level'), size = 1.2) +
#     geom_point(data = drop_na(plot_df), aes(y = pred, color = 'Estimated Water Level'), size = 1.2) +
#     facet_wrap(~name, ncol = 1, strip.position = 'left',
#                labeller = label_parsed, scales = 'free_y') +
#     scale_y_continuous(breaks = pretty_breaks(5),
#                        labels = scales::label_number(accuracy = 0.01, scale_cut = cut_long_scale())) +
#     scale_x_date(breaks = pretty_breaks(n = 10)) +
#     scale_color_manual(name = NULL, 
#                        values = c('Observed Water Level' = 'red', 
#                                   'Estimated Water Level' = 'black')) +
#     theme_bw(base_size = 12, base_family = 'serif') +
#     theme(strip.placement = 'outside', 
#           strip.background = element_blank(), 
#           axis.text = element_text(color = 'black'), 
#           legend.position = 'bottom') +
#     labs(x = 'Date', y = '') 
#   
# }
# 
# sgd_mod <- create_sgd_model(test_data, tanya_params, const_params)
# sgd_mod$calc_recharge()
# sgd_mod$calc_sgd(365)
# model_diagnosis(sgd_mod, obs_dt)
# 
# cl <- parallel::makeCluster(16) # this is a little time consuming (to be fixed)
# 
# packFn <- function(packages) {
#   for (i in packages) library(i, character.only = TRUE)
# }
# parallel::clusterCall(cl, packFn, c('Rcpp', 'tidyverse', 'pracma', 'magrittr', 'DEoptim', 'SGDr'))
# parallel::clusterExport(cl = cl, varlist = ls(), 
#                         envir = environment())
# 
# 
# optDE <- DEoptim2(fn = obj_fun,
#                   weights = ws,
#                   obs = obs_wl, 
#                   xn = 3000, hn = 23, 
#                   default_values = preferred_params,
#                   lower = lowerB, upper = upperB,
#                   skeleton = params_skeleton,
#                   control = DEoptim.control(NP = NA,
#                                             trace = 50, 
#                                             itermax = 500, 
#                                             parallelType = 'parallel', 
#                                             cluster = cl
#                   ))
# opt_param <- optDE$optim$bestmem %>% relist(skeleton = params_skeleton)
# opt_mod <- create_sgd_model(test_data, opt_param, const_params)
# opt_mod$calc_recharge()
# opt_mod$calc_sgd(180)
# model_diagnosis(opt_mod, obs_dt)
