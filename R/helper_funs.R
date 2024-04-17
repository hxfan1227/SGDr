#' @import jsonlite scales dplyr tidyr ggplot2 patchwork data.table
#' @importFrom lubridate ymd date days_in_month
#' @importFrom glue glue

NULL 
#'  Convert a JSON file to a list of parameters
#' @name json_to_paramter_list
#' @aliases json_to_paramter_list
#' @rdname json_to_paramter_list
#' @description This function converts a JSON file to a list of parameters.
#' @param json_file A JSON file containing the parameters
#' @return A list of parameters (invisible)
#' @export

json_to_parameter_list <- function(json_file) {
    invisible(as.relistable(jsonlite::fromJSON(json_file)))
}


#' Convert a vector of parameters to a list of parameters (useful for calibrating parameters)
#' @name paramter_vec_to_list
#' @aliases paramter_vec_to_list
#' @rdname paramter_vec_to_list
#' @description This function converts a vector of parameters to a list of parameters.
#' @param parameter_vec A vector of parameters
#' @param prameter_skeleton A list of parameters to use as a skeleton for the new list
#' @return A list of parameters (invisible)
#' @export

parameter_vec_to_list <- function(parameter_vec, prameter_skeleton) {
    invisible(relist(parameter_vec, skeleton = prameter_skeleton))
}

#' Convert a list of parameters to a vector of parameters
#' @name parameter_list_to_vector
#' @aliases parameter_list_to_vector
#' @rdname parameter_list_to_vector
#' @description This function converts a list of parameters to a vector of parameters.
#' @param parameter_list A list of parameters
#' @return A named vector of parameters
#' @export

parameter_list_to_vector <- function(parameter_list) {
  params <- unlist(parameter_list) %>% as.vector()
  names(params) <- names(unlist(parameter_list))
  params
}


#' Print the SGD estimation results.
#' @rdname print.SGD_ESTIMATION_DF
#' @param x An object of class \code{SGD_ESTIMATION_DF}. Created by \code{\link{estimate_sgd}}.
#' @param ... other arguments not used by this method
#' @export
print.SGD_ESTIMATION_DF <- function(x, ...) {
    summary(x, ...)
}


#' Summary the SGD estimation results.
#' @rdname summary
#' @param object An object of class \code{sgdEstimation}. Created by \code{\link{estimate_sgd}}.
#' @param base_date A valid string that can be transformed into data format using \code{\link[lubridate]{ymd}}.
#' @param ... other arguments not used by this method
#' @export
summary.SGD_ESTIMATION_DF <- function(object, base_date, ...) {
  cat("SGD Estimation Summary\n")
  if (missing(base_date)) {
    warning("No base date available. Please provide the base date to calculate the simulation period. \n")
    cat("Simulation length: ", NROW(object$results), 'days\n')
  } else {
    cat(glue("Simlation period: {date(ymd(base_date))} - {date(ymd(base_date) + days(NROW(object$results) - 1))}"), "\n")
  }
  cat("Estimated SGD Volume: ", mean(object$results$SGD, na.rm = T), "m3/d\n")
  cat("Estimated Recharge: ", mean(object$results$Wrechg, na.rm = T), "mm/d\n")
  cat("Estimated Runoff: ", mean(object$runoff, na.rm = T), "mm/d\n")
  cat("Median hn:", median(object$results$hn, na.rm = T), "m\n")
  cat("Median xn:", median(object$results$xn, na.rm = T), "m\n")
}

NULL 
#' Plot the SGD estimation results
#' @rdname plot.SGD_ESTIMATION_DF
#' @param x An object of class \code{sgdEstimation}. Created by \code{\link{estimate_sgd}}.
#' @param y A data.frame of daily observed groundwater level data. See Details.
#' @param vars A character vector of the variables to plot. Can be 'all' or a subset of 'wl', 'SGD', 'hn', 'xn', 'Wrechg', 'WrechgAve'.
#' @param base_date A valid string that can be transformed into data format using \code{\link[lubridate]{ymd}}.
#' @param type A character string indicating the type of plot. Can be one of 'pred', 'comp', 'input'.
#' @param obs_x A character string indicating the column name of the date in the observed data.
#' @param obs_y A character string indicating the column name of the groundwater level in the observed data.
#' @param xbreaks An integer indicating the number of breaks on the x-axis.
#' @param ybreaks An integer indicating the number of breaks on the y-axis.
#' @param names A character vector of the variable names to plot.
#' @param labels A character vector of the variable labels to plot.
#' @param maxRange A numeric indicating the maximum range of the first axis (groundwater level).
#' @param coeff A numeric indicating the shrink coefficient of the precipitation.
#' @param colors A named character vector of the colors to use
#' @param color_labels A named character vector of the labels for the colors.
#' @param ... Additional arguments passed to the plot function. (currently not implemented)
#' @return Invisibly returns the original estimated results.
#' @export
plot.SGD_ESTIMATION_DF <- function(x, y, 
                                   vars = c('wl', 'SGD', 'Wrechg'), 
                                   base_date, 
                                   type = c('pred', 'comp', 'input'),
                                   obs_x = 'date',
                                   obs_y = 'wl',
                                   xbreaks = 10,
                                   ybreaks = 5,
                                   maxRange = 600,
                                   coeff = 0.3,
                                   names = c('wl', 'SGD', 'hn', 'xn', 'Wrechg', 'WrechgAve'),
                                   colors = c('obs' = 'red', 'precp' = 'blue' , 'pumping' = 'green', 'est' = 'blue'),
                                   color_labels = c('obs' = 'Observed GWL', 'precp' = 'Precipitation', 
                                                    'pumping' = 'Pumping', 'est' = 'Estimated GWL'),
                                   labels = c('WL~(m~AHD)',
                                              'SGD~(m^3/d)',
                                              'h[n]~(m~AHD)',
                                              'x[n]~(m)',
                                              'w[recharge]~(mm)',
                                              'w[recharge[ave]]~(mm)'),
                                   ...) {
  args <- list(...)
  results <- x$results
  if (all(vars == 'all')) {
    vars <- colnames(results)[!colnames(results) %in% c('t')]
  }
  stopifnot(all(vars %in% names), length(names) == length(labels))
  type = match.arg(type)
  if (type == 'pred') {
    plot_df <- data.table::setDT(results)
    plot_df[, date := ymd(base_date) + t - 1]
    plot_df[, t := NULL]
    plot_df <- data.table::melt(plot_df, id.vars = 'date', measure.vars = vars)
    plot_df[, variable := factor(variable, levels = names, labels = labels)]
    p <- ggplot(plot_df, aes(x = date, y = value)) +
      geom_line(linewidth = 0.8) +
      facet_wrap(~ variable, ncol = 1, strip.position = 'left',
                 labeller = label_parsed, scales = 'free_y') +
      scale_y_continuous(breaks = scales::pretty_breaks(ybreaks),
                         labels = scales::label_number(accuracy = 0.01, scale_cut = scales::cut_si(''))) +
      scale_x_date(name = NULL, breaks = scales::pretty_breaks(n = xbreaks)) +
      theme_bw(base_size = 12, base_family = 'serif') +
      theme(strip.placement = 'outside',
            strip.background = element_blank(),
            axis.text = element_text(color = 'black'),
            legend.position = 'bottom') +
      labs(x = 'Date', y = '')
  }
  
  if (type == 'comp') {
    if(missing(y)) {
      stop('You must provide the observed data')
    }
    obs_df <- data.table::setDT(y)
    plot_df <- data.table::setDT(results)
    plot_df[, date := ymd(base_date) + t - 1]
    plot_df[, t := NULL]
    p1 <- ggplot() +
      geom_line(data = plot_df, aes(x = date, y = wl, color = 'est'), linewidth = 0.8, show.legend = F) +
      geom_point(data = plot_df[date %in% obs_df[[obs_x]]], aes(x = date, y = wl, color = 'est'), size = 1.2, show.legend = F) +
      geom_point(data = obs_df, aes(x = .data[[obs_x]], y = .data[[obs_y]], color = 'obs'), size = 1.2, show.legend = F) +
      scale_y_continuous(name = parse(text = "GWL~(m~AHD)"),
                         expand = c(0, 0),
                         limits = args$gwl_range,
                         breaks = scales::pretty_breaks(ybreaks)) +
      scale_x_date(name = NULL, breaks = scales::pretty_breaks(n = xbreaks), expand = c(0, 0)) +
      scale_color_manual(name = NULL, values = colors, labels = color_labels) +
      guides(color = guide_legend(nrow = 1), fill =guide_legend(nrow = 1)) +
      theme_bw(base_size = 12, base_family = 'serif') +
      theme(strip.placement = 'outside',
            strip.background = element_blank(),
            axis.text = element_text(color = 'black'),
            legend.position = 'top')
    # scatter plot
    plot_df <- plot_df %>% 
      rename('est' = 'wl')
    obs_df <- obs_df %>%
      rename('obs' = all_of(obs_y))
    scatter_df <- merge(plot_df[, .(date, est)], obs_df, by = c('date' = obs_x))
    p2 <- ggplot(scatter_df, aes(x = obs, y = est)) +
      geom_point(size = 1.2) +
      geom_abline(intercept = 0, slope = 1, linetype = 2) +
      scale_x_continuous(name = parse(text = "Observed~GWL~(m~AHD)"),
                         expand = c(0, 0),
                         breaks = scales::pretty_breaks(ybreaks)) +
      scale_y_continuous(name = parse(text = "Estimated~GWL~(m~AHD)"),
                         expand = c(0, 0),
                         breaks = scales::pretty_breaks(ybreaks)) +
      theme_bw(base_size = 12, base_family = 'serif')
    p3 <- scatter_df %>% 
      pivot_longer(cols = c('obs', 'est')) %>%
      ggplot(aes(x = value, color = name)) + 
      stat_ecdf(geom = 'line', linewidth = 0.8) +
      scale_y_continuous(name = 'Cumulative probability (-)')+
      scale_x_continuous(name = parse(text = "GWL~(m~AHD)"),
                         expand = c(0, 0),
                         breaks = scales::pretty_breaks(ybreaks)) +
      scale_color_manual(name = NULL, values = colors, labels = color_labels) +
      guides(color = guide_legend(nrow = 1), fill =guide_legend(nrow = 1)) +
      theme_bw(base_size = 12, base_family = 'serif') +
      theme(strip.placement = 'outside',
            strip.background = element_blank(),
            axis.text = element_text(color = 'black'),
            legend.position = 'top')
    p <- p1 / (p2 + p3) +
      patchwork::plot_layout(guides = 'collect') &
      theme(legend.position = 'top')
  }
  
  if (type == 'input') {
    plot_df <- data.table::setDT(x$input)
    plot_df[, date := ymd(base_date) + t - 1]
    p <- ggplot(plot_df) +
      geom_segment(aes(x = date, y = maxRange, yend = maxRange - R/coeff, xend = date, color = 'precp'), linewidth = 0.5) +
      geom_line(aes(x = date, y = Pumping, color = 'pumping'), linewidth = 0.5) +
      scale_y_continuous(name = parse(text = "Pumping~(m^3/d)"),
                         limit = c(0, maxRange),
                         expand = c(0, 0),
                         breaks = scales::pretty_breaks(ybreaks),
                         sec.axis = sec_axis(~ (maxRange - .) * coeff,
                                             name = expression(Precipitation~(mm/d)), 
                                             breaks = scales::pretty_breaks(ybreaks))) +
      scale_x_date(name = NULL, breaks = scales::pretty_breaks(n = xbreaks), expand = c(0, 0)) +
      scale_color_manual(name = NULL, values = colors, labels = color_labels) +
      guides(color = guide_legend(nrow = 1), fill =guide_legend(nrow = 1)) +
      theme_bw(base_size = 12, base_family = 'serif') +
      theme(strip.placement = 'outside',
            strip.background = element_blank(),
            axis.text = element_text(color = 'black'),
            legend.position = 'top',
            axis.title.y.right = element_text(color = colors['precp']),
            axis.line.y.right = element_line(color = colors['precp']), 
            axis.text.y.right = element_text(color = colors['precp']))
  }
  print(p)
  invisible(x)
}

#' A wrapper function to call the model and return result for calibration.
#' @rdname estimate_sgd_from_pars
#' @param pars A numeric vector of parameters to calibrate
#' @param parnames A character vector of names of the parameters
#' @param parset A character vector of names of the to-be-calibrated parameters
#' @param input A data.table of the input data
#' @param warm_up A numeric value of the warm-up period
#' @param default_pars A named numeric vector of the default parameters
#' @param skeleton A named list of the parameter skeleton
#' @param const_par_list A named list of the constant parameters
#' @return A SGD_ESTIMATION_DF class
#' @seealso [estimate_sgd()]
#' @export

estimate_sgd_from_pars <- function(pars, 
                                   parnames = c(params_to_calibrate, 'nw'),
                                   parset = params_to_calibrate, 
                                   skeleton,
                                   inputDf,
                                   yearlyPumping,
                                   pumpingDf,
                                   ...) {
  args <- list(...)
  if(length(pars) != length(parnames)) {
    stop('The length of the parameters and the parameter names should be the same.')
  }
  if (!is.null(names(pars)) & !all(parnames == names(pars))) {
    stop('A named parameter vector is supplied, but the names do not match the parnames.')
  }
  names(pars) <- parnames
  current_pars = parameter_list_to_vector(args$calibratableParams)
  if (any(!parset %in% names(current_pars))) {
    stop('The parset should be a subset of the names of the calibratable parameters.')
  }
  current_pars[parset] = pars[parset]
  current_pars_list <- parameter_vec_to_list(current_pars, prameter_skeleton = skeleton)
  # additional 'calibratable' parameter nw (i.e., window size for averaging recharge)
  default_args <- formals(estimate_sgd)
  nw = default_args$windowSize
  if ('nw' %in% names(pars)) {
    nw = pars['nw']
  } 
  if ('pumping' %in% names(pars)) {
    inputDf <- change_unknown_pumping(pars['pumping'], inputDf, yearlyPumping, pumpingDf)
  }
  return(estimate_sgd(inputData = inputDf,  windowSize = nw, calibratableParams = current_pars_list, 
                      constParams = args$constParams, warmUp = args$warmUp))
}

#' A function to change the unknown pumping rate to a preferred value.
#' @rdname change_unknown_pumping
#' @param x A numeric value of the preferred pumping rate
#' @param input_df A data.frame of the input data. See \link{estimate_sgd} for details.
#' @param yearly_df A data.frame of the yearly pumping data. (2 columns: year, pumping)
#' @param pumping_df A data.frame of the pumping data
#' @return A data.frame of the input data with the preferred pumping rate
#' @export
change_unknown_pumping <- function(x, input_df, yearly_df, pumping_df) {
  input_data <- setDT(copy(input_df))
  yearly_pumping_data <- setDT(copy(yearly_df))
  pumping_data <- setDT(copy(pumping_df))
  yearly_pumping_data[year <= 1985, pumping := x]
  pumping_data <- pumping_data[yearly_pumping_data, pumping := i.pumping, on = .(year)]
  pumping_data[, daily_pumping := pumping * percent * 1000 / days_in_month(date)]
  input_data[pumping_data, Pumping := i.daily_pumping, on = .(t)]
  setDF(input_data)
  return(input_data)
}
