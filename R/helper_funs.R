#' @import jsonlite scales dplyr tidyr data.table ggplot2 patchwork
#' @importFrom lubridate ymd date
#' @importFrom glue glue

NULL 
#' @name json_to_paramter_list
#' @aliases json_to_paramter_list
#' @title Convert a JSON file to a list of parameters
#' @rdname json_to_paramter_list
#' @description This function converts a JSON file to a list of parameters.
#' @param json_file A JSON file containing the parameters
#' @return A list of parameters (invisible)
#' @export

json_to_paramter_list <- function(json_file) {
    invisible(as.relistable(jsonlite::fromJSON(json_file)))
}

NULL 
#' @name paramter_vec_to_list
#' @aliases paramter_vec_to_list
#' @title Convert a vector of parameters to a list of parameters (useful for calibrating parameters)
#' @rdname paramter_vec_to_list
#' @description This function converts a vector of parameters to a list of parameters.
#' @param parameter_vec A vector of parameters
#' @param prameter_skeleton A list of parameters to use as a skeleton for the new list
#' @return A list of parameters (invisible)
#' @export

paramter_vec_to_list <- function(parameter_vec, prameter_skeleton) {
    invisible(relist(parameter_vec, skeleton = prameter_skeleton))
}


NULL
#' @title Print the SGD estimation results.
#' @rdname print
#' @export

print <- function(x, ...) {
    UseMethod("print")
}

#' @rdname print
#' @param x An object of class \code{SGD_ESTIMATION_DF}. Created by \code{\link{estimate_sgd}}.
#' @export
print.SGD_ESTIMATION_DF <- function(x, ...) {
    summary(x, ...)
}

NULL 
#' @title Summary the SGD estimation results.
#' @rdname summary
#' @export
summary <- function(object, ...) {
    UseMethod("summary")
}

#' @rdname summary
#' @param object An object of class \code{sgdEstimation}. Created by \code{\link{estimate_sgd}}.
#' @param base_date A valid string that can be transformed into data format using \code{\link[lubridate]{ymd}}.
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
#' @title Plot the SGD estimation results
#' @rdname plot
#' @export
plot <- function(x, y, ...) {
  UseMethod("plot")
}

#' @rdname plot
#' @param x An object of class \code{sgdEstimation}. Created by \code{\link{estimate_sgd}}.
#' @param y A data.frame of daily observed groundwater level data. See Details.
#' @param base_date A valid string that can be transformed into data format using \code{\link[lubridate]{ymd}}.
#' @param type A character string indicating the type of plot. Can be one of 'pred', 'obs', 'comp', 'input'.
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
#' @export
plot.SGD_ESTIMATION_DF <- function(x, y, 
                                   vars = c('wl', 'SGD', 'Wrechg'), 
                                   base_date, 
                                   type = c('pred', 'obs', 'comp', 'input'),
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
      rename('obs' = obs_y)
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
  p
}