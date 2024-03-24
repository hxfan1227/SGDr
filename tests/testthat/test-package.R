# library(SGDr)
# 
# test_data <- read.csv('/mnt/c/RProjects/DailyWork/CoffinBay/Scripts/SGD/Data/test_data.csv')
# 
# 
# const_params <- list(
#   x = 3253, 
#   W = 5864,
#   Area = 36885913,
#   L = 6760,
#   waterContent = list(Texture = c("Sand", "Loam", "Clay"), 
#                       ClayContent = c(3, 22, 47), 
#                       Saturation = c(0.4, 0.5, 0.6), 
#                       FC = c(0.06, 0.29, 0.41), 
#                       PWP = c(0.02, 0.05, 0.2)
#   )
# )
# 
# jsonlite::write_json(const_params, './tests/testthat/consts.json', pretty = T)
# 
# tanya_params <- json_to_paramter_list('./tests/testthat/tanya.json')
# const_params <- json_to_paramter_list('./tests/testthat/consts.json')
# preferred_params <- json_to_paramter_list('./tests/testthat/preferred.json')
# 
# sgd_mod <- create_sgd_model(test_data, tanya_params, const_params)
# sgd_mod$calc_recharge()
# sgd_mod$calc_sgd(30)
# plot(sgd_mod$get_H2O3_AQ(), ylim = c(0, 2))
# sgd_mod$update_parameters(preferred_params, const_params)
# sgd_mod$calc_recharge()
# sgd_mod$calc_sgd(180)
# lines(sgd_mod$get_H2O3_AQ(), col = 'red')
# 
# plot(sgd_mod$get_xn(), col = 'blue')
# 
# 
