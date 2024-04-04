# # Load the preferred parameters 
# preferred_params_list <- json_to_paramter_list('/mnt/c/RProjects/DailyWork/CoffinBay/Scripts/SGD/Data/json/preferred.json')
# # Load Tanya's calibrated parameters
# tanya_params_list <- json_to_paramter_list('/mnt/c/RProjects/DailyWork/CoffinBay/Scripts/SGD/Data/json/tanya.json')
# # Load the constant parameters
# const_params_list <- json_to_paramter_list('/mnt/c/RProjects/DailyWork/CoffinBay/Scripts/SGD/Data/json/consts.json')
# # Load the lower bounds of the parameters
# lower_bounds_list <- json_to_paramter_list('/mnt/c/RProjects/DailyWork/CoffinBay/Scripts/SGD/Data/json/lower_limits.json')
# # Load the upper bounds of the parameters
# upper_bounds_list <- json_to_paramter_list('/mnt/c/RProjects/DailyWork/CoffinBay/Scripts/SGD/Data/json/upper_limits.json')
# # Load the maximum bounds of the parameters
# max_bounds_list <- json_to_paramter_list('/mnt/c/RProjects/DailyWork/CoffinBay/Scripts/SGD/Data/json/maximum.json')
# # Load the minimum bounds of the parameters
# min_bounds_list <- json_to_paramter_list('/mnt/c/RProjects/DailyWork/CoffinBay/Scripts/SGD/Data/json/minimum.json')
# # Load the skeleton of the parameters
# params_skeleton <- json_to_paramter_list('/mnt/c/RProjects/DailyWork/CoffinBay/Scripts/SGD/Data/json/preferred.json')
# # Load the observed data
# obs_data <- read_csv('/mnt/c/RProjects/DailyWork/CoffinBay/Scripts/SGD/Data/obs_data.csv') %>% 
#   mutate(date = dmy(obs_date)) %>% 
#   select(date, wl) %>% 
#   distinct(date, .keep_all = T) %>% 
#   filter(wl <= 2) # to eliminate the outliers
# # load the input data
# test_data <- read_csv('/mnt/c/RProjects/DailyWork/CoffinBay/Scripts/SGD/Data/test_data.csv')
# ids <- as.integer(obs_data$date - ymd('1967-01-01') + 1)

# obs_dt <- read_csv('/mnt/c/RProjects/DailyWork/CoffinBay/Scripts/SGD/Data/obs_data.csv') %>% 
#   select(obs_date, wl) %>%
#   filter(wl <= 2) %>%
#   rename('date' = 'obs_date') %>%
#   pivot_longer(-date, values_to = 'obs') %>%
#   mutate(name = factor(name, levels = c('wl', 'SGD', ' hn', 'xn', 'Wrechg', 'WrechgAve'),
#                        labels = c('WL~(m~AHD)',
#                                   'SGD~(m^3/d)',
#                                   'h[n]~(m~AHD)',
#                                   'x[n]~(m)',
#                                   'w[recharge]~(mm)',
#                                   'w[recharge[ave]]~(mm)')),
#          date = dmy(date))

# x1 <- estimate_sgd(test_data, tanya_params_list, const_params_list, 30)
# x2 <- estimate_sgd_fixed_hxn(test_data, tanya_params_list, const_params_list, 30, xn = 7000, hn = 28)
# plot(x1$wl, type = 'l', lty = 2)
# lines(x2$wl, col = 'red')
# legend('topright', c("Without fixed hn and xn", "With fixed hn and xn"),
#        col = c('black', 'red'), lty = 2)
# model_diagnose(x2, obs_dt, '1967-01-01')
