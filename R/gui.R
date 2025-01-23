#' @import shiny
#' @import shinyFiles
#' @import plotly
#' @importFrom glue glue
#' @importFrom shinycssloaders hideSpinner showSpinner withSpinner
#' @importFrom shiny.fluent DatePicker.shinyInput
#' @importFrom purrr walk

withConsoleRedirect <- function(containerId, expr) {
  # Change type="output" to type="message" to catch stderr
  # (messages, warnings, and errors) instead of stdout.
  txt <- capture.output(results <- expr, type = "output")
  if (length(txt) > 0) {
    insertUI(paste0("#", containerId), where = "beforeEnd",
             ui = paste0(txt, "\n", collapse = "")
    )
  }
  results
}

ui <- fluidPage(
  shinyjs::useShinyjs(),
  # theme = bslib::bs_theme(bootswatch ='flatly', version = '3'),
  titlePanel("SFGD calculator V1.0"),
  # chooseSliderSkin("Flat", color = "#112446"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Simulation start date:"),
      shiny.fluent::DatePicker.shinyInput('start_date', value = "1967-01-01T00:00:00.000Z"),
      h4("Input files:"),
      shinyFilesButton("bound_input", "Boundaries", "Please choose the input file for the boundary conditions", F, icon = icon("upload")),
      shinyFilesButton("json_input", "Parameters", "Please choose the json files for the parameters", T, icon = icon("upload")),
      shinyFilesButton("obs_input", "Observations", "Please choose the input file for the well observations", F, icon = icon("upload")),
      br(),
      br(),
      conditionalPanel(
        condition = "output.is_par_uploaded", 
        div(id = 'param_input',
            tabsetPanel(
              tabPanel("Bucket 1",
                       sliderInput("bucket1_z", withMathJax("\\(z~\\mathrm{[mm]}\\)"), min = 0, max = 1000, value = 100, step = 0.1, animate = T),
                       # bsTooltip("bucket1_z", "Depth of soil bucket 1", placement = "top"),
                       sliderInput("bucket1_rho_b", withMathJax("\\(\\rho_b~\\mathrm{[Mg/m^3]}\\)"), min = 0, max = 10, value = 1.56, step = 0.01, animate = T),
                       # bsTooltip("bucket1_rho_b", "Bulk density", placement = "top"),
                       sliderInput("bucket1_m", withMathJax("\\(m~\\mathrm{[\\%]}\\)"), min = 0, max = 100, value = 11.1, step = 0.01, animate = T),
                       sliderInput("bucket1_rho_s", withMathJax("\\(\\rho_s~\\mathrm{[Mg/m^3]}\\)"), min = 0, max = 10, value = 2.72, step = 0.01, animate = T),
                       sliderInput("bucket1_Ks", withMathJax("\\(K_{SAT}~\\mathrm{[mm/hr]}\\)"), min = 0, max = 1000, value = 75.53, step = 0.01, animate = T),
                       sliderInput("bucket1_n", withMathJax("\\(\\mathrm{van~Genuchten}~n~\\mathrm{[-]}\\)"), min = 0, max = 10, value = 1.52, step = 0.01, animate = T),
                       sliderInput("bucket1_pE", withMathJax("\\(p_{E}~\\mathrm{[\\%]}\\)"), min = 0, max = 100, value = 95, step = 0.1, animate = T)
              ),
              tabPanel("Bucket 2",
                       sliderInput("bucket2_z", withMathJax("\\(z~\\mathrm{[mm]}\\)"), min = 0, max = 1000, value = 1000, step = 0.1, animate = T),
                       sliderInput("bucket2_rho_b", withMathJax("\\(\\rho_b~\\mathrm{[Mg/m^3]}\\)"), min = 0, max = 10, value = 1.56, step = 0.01, animate = T),
                       sliderInput("bucket2_m", withMathJax("\\(m~\\mathrm{[\\%]}\\)"), min = 0, max = 100, value = 11.1, step = 0.01, animate = T),
                       sliderInput("bucket2_rho_s", withMathJax("\\(\\rho_s~\\mathrm{[Mg/m^3]}\\)"), min = 0, max = 10, value = 2.72, step = 0.01, animate = T),
                       sliderInput("bucket2_Ks", withMathJax("\\(K_{SAT}~\\mathrm{[mm/hr]}\\)"), min = 0, max = 1000, value = 75.53, step = 0.01, animate = T),
                       sliderInput("bucket2_n", withMathJax("\\(\\mathrm{van~Genuchten}~n~\\mathrm{[-]}\\)"), min = 0, max = 10, value = 1.52, step = 0.01, animate = T),
                       sliderInput("bucket2_pZ", withMathJax("\\(p_{Z}~\\mathrm{[\\%]}\\)"), min = 0, max = 100, value = 90, step = 0.1, animate = T)
              ),
              tabPanel("Aquifer",
                       sliderInput("aquifer_Td", withMathJax("\\(\\delta~\\mathrm{{d}}\\)"), min = 0, max = 31, value = 1, step = 0.01, animate = T),
                       sliderInput("aquifer_rho_s", withMathJax("\\(\\rho_s~\\mathrm{[kg/m^3]}\\)"), min = 1000, max = 1030, value = 1025, step = 0.0, animate = T),
                       sliderInput("aquifer_rho_f", withMathJax("\\(\\rho_f~\\mathrm{[kg/m^3]}\\)"), min = 990, max = 1000, value = 998.949, step = 0.01, animate = T),
                       sliderInput("aquifer_Ka", withMathJax("\\(K_{a}~\\mathrm{[m/d]}\\)"), min = 0, max = 3000, value = 200, step = 1, animate = T),
                       sliderInput("aquifer_z0", withMathJax("\\(z_0~\\mathrm{[m]}\\)"), min = 0, max = 100, value = 30, step = 0.1, animate = T),
                       sliderInput("aquifer_Sy", withMathJax("\\(Sy~\\mathrm{[-]}\\)"), min = 0, max = 1, value = 0.172, step = 0.01, animate = T),
                       sliderInput("aquifer_xT", withMathJax("\\(x_{T}~\\mathrm{[m]}\\)"), min = 0, max = 10000, value = 3000, step = 10, animate = T),
                       sliderInput("aquifer_dxT", withMathJax("\\(\\delta x_{T, max}~\\mathrm{[m/d]}\\)"), min = 0, max = 1, value = 0.2, step = 0.01, animate = T)
              ),
              tabPanel("Others",
                       sliderInput("curveNumber_CN2", withMathJax("\\(CN2~\\mathrm{[-]}\\)"), min = 0, max = 100, value = 75, step = 0.1, animate = T),
                       sliderInput("curveNumber_pF", withMathJax("\\(p_{F}~\\mathrm{[-]}\\)"), min = 0, max = 1, value = 0.8, step = 0.01, animate = T),
                       sliderInput("windowSize", withMathJax("\\(n_w~\\mathrm{[d]}\\)"), min = 30, max = 180, value = 120, step = 1, animate = T),
                       # sliderInput("Pumping", withMathJax("\\(Pumping~\\mathrm{[m^3/d]}\\)"), min = 10, max = 60, value = 50, step = 1, animate = T)
              ),
              # tabPanel("Constants",
              #          textInput("const_x", withMathJax("\\(x~\\mathrm{[m]}\\)"), value = 3253),
              #          textInput("const_W", withMathJax("\\(W~\\mathrm{[m]}\\)"), value = 5864),
              #          textInput("const_A", withMathJax("\\(A~\\mathrm{[m^2]}\\)"), value = 36885913),
              #          textInput("const_L", withMathJax("\\(L~\\mathrm{[m]}\\)"), value = 6760)
              # )
            )
        ),
        div(actionButton("reset", "Reset"))
      ),
      h4('Logs:'),
      pre(id = 'log', '')
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Simulation Plots", 
                 fluidRow(
                   column(width = 9, shinycssloaders::withSpinner(plotlyOutput("wlPlot"))),
                   column(width = 3, tableOutput('wlMetrics'))
                 ),
                 fluidRow(
                   # column(width = 6, plotlyOutput("wrPlot")),
                   column(width = 9, plotlyOutput("q0Plot"))
                 )
        ),
        tabPanel("Results",
                 DT::DTOutput("resultTable"))
      ), width = 9)
  )
)


server <- function(input, output, session) {

  shinycssloaders::hideSpinner("wlPlot")
  shinycssloaders::hideSpinner("q0Plot")
  roots <- getVolumes()() 
  roots <- c(wd='.')
  
  vars <- reactiveValues()
  
  shinyFileChoose(input, "obs_input", roots = roots, filetypes=c('', 'csv'))
  shinyFileChoose(input, "bound_input", roots = roots, filetypes=c('', 'csv'))
  shinyFileChoose(input, "json_input", roots = roots, filetypes=c('', 'json'))
  
  observeEvent(input$obs_input, {
    inFile <- parseFilePaths(roots, input$obs_input)
    if (NROW(inFile) == 0) {
      return(NULL)
    }
    vars$obs_data <- read.csv(inFile$datapath) %>% 
      mutate(date = ymd(date))
    withConsoleRedirect("log", {
      walk(inFile$datapath, ~ cat(glue('Finish reading {.x}'), '\n'))
    }) 
    plotlyProxy("wlPlot", session) %>% 
      plotlyProxyInvoke("addTraces", list(
        x = vars$obs_data$date,
        y = vars$obs_data$hf,
        type = "scatter",
        mode = "markers",
        marker = list(color = "red"),
        name = "Observed GWL",
        yaxis = "y3"   # Ensure it targets the second subplot's y-axis
      ))
  })
  output$is_obs_uploaded <- reactive({
    !is.null(vars$obs_data)
  })
  outputOptions(output, "is_obs_uploaded", suspendWhenHidden = FALSE)
  
  observeEvent(input$bound_input, {
    inFile <- parseFilePaths(roots, input$bound_input)
    if (NROW(inFile) == 0) {
      return(NULL)
    }
    vars$bound_data <- read.csv(inFile$datapath) %>% 
      mutate(date = as.Date(input$start_date) - 1 + t)
    withConsoleRedirect("log", {
      walk(inFile$datapath, ~ cat(glue('Finish reading {.x}'), '\n'))
    }) 
  })
  output$is_bound_uploaded <- reactive({
    !is.null(vars$bound_data)
  })
  outputOptions(output, "is_bound_uploaded", suspendWhenHidden = FALSE)
  
  observeEvent(input$json_input, {
    inFile <- parseFilePaths(roots, input$json_input)
    if (NROW(inFile) == 0) {
      return(NULL)
    }
    vars$cali_par <- inFile %>% 
      filter(name == 'preferred.json') %>% 
      pull('datapath') %>% 
      SGDr::json_to_parameter_list()
    vars$cnst_par <- inFile %>% 
      filter(name == 'consts.json') %>% 
      pull('datapath') %>% 
      SGDr::json_to_parameter_list()
    withConsoleRedirect("log", {
      walk(inFile$datapath, ~ cat(glue('Finish reading {.x}'), '\n'))
    }) 
    showSpinner('wlPlot')
    updateSliderInput(session, "bucket1_z", value = vars$cali_par$bucket1$z)
    updateSliderInput(session, "bucket1_rho_b", value = vars$cali_par$bucket1$rho_b)
    updateSliderInput(session, "bucket1_m", value = vars$cali_par$bucket1$m)
    updateSliderInput(session, "bucket1_rho_s", value = vars$cali_par$bucket1$rho_s)
    updateSliderInput(session, "bucket1_Ks", value = vars$cali_par$bucket1$Ks)
    updateSliderInput(session, "bucket1_n", value = vars$cali_par$bucket1$n)
    updateSliderInput(session, "bucket1_pE", value = vars$cali_par$bucket1$pE)
    updateSliderInput(session, "bucket2_z", value = vars$cali_par$bucket2$z)
    updateSliderInput(session, "bucket2_rho_b", value = vars$cali_par$bucket2$rho_b)
    updateSliderInput(session, "bucket2_m", value = vars$cali_par$bucket2$m)
    updateSliderInput(session, "bucket2_rho_s", value = vars$cali_par$bucket2$rho_s)
    updateSliderInput(session, "bucket2_Ks", value = vars$cali_par$bucket2$Ks)
    updateSliderInput(session, "bucket2_n", value = vars$cali_par$bucket2$n)
    updateSliderInput(session, "bucket2_pZ", value = vars$cali_par$bucket2$pZ)
    updateSliderInput(session, "curveNumber_CN2", value = vars$cali_par$cn$CN2)
    updateSliderInput(session, "curveNumber_pF", value = vars$cali_par$cn$pF)
    
    withConsoleRedirect("log", {
      tictoc::tic(msg = 'SFGD estimation')
      vars$sgd_res <- estimate_sgd(vars$bound_data, vars$cali_par, vars$cnst_par, 
                                   nw = input$windowSize)
      tictoc::toc()
      vars$sgd_res$results <- vars$sgd_res$results %>% 
        mutate(date = as.Date(input$start_date) - 1 + t)
    }) 
    output$wlPlot <- renderPlotly({
      p1 <- plot_ly(data = vars$bound_data) %>% 
        add_bars(x = ~ date, y = ~ R, yaxis = 'y2', name = 'Precipitation', color = I('blue')) %>% 
        add_lines(x = ~ date ,y = ~ q1, yaxis = 'y1', name = 'Pumping', color = I('green')) %>%
        layout(barmode = 'group', 
               xaxis = list(title = 'Date', autorange = T, automargin = T),
               yaxis = list(title = 'Pumping (m<sup>3</sup>/d)', side = "left", range = c(75, 600)),
               yaxis2 = list(title ='Precpitation (mm/d)', automargin = T,
                             overlaying = "y", autorange='reversed', side = "right"))
      p2 <- plot_ly() %>% 
        # add_markers(data = obs_data, x = ~ date, y = ~ hf, yaxis = 'y1', name = 'Observed GWL', color = I('red')) %>%
        add_lines(data = vars$sgd_res$results, x = ~ date, y = ~ hf, yaxis = 'y1', name = 'Preferred GWL', color = I('grey')) %>% 
        add_lines(data = vars$sgd_res$results, x = ~ date, y = ~ hf, yaxis = 'y1', name = 'Current GWL', color = I('orange'))
      
      subplot(p1, p2, nrows = 2, margin = 0.05, shareX = T, titleY = T) %>% 
        layout(xaxis = list(title = '', autorange = T, automargin = T),
               legend = list(y = 0.5, orientation = 'h', x = 0.45, bgcolor = 'rgba(255, 255, 255, 0.5)', 
                             xanchor = 'center', yanchor = 'bottom', yref = 'container'))
    })
    output$q0Plot <- renderPlotly({
      plot_ly() %>% 
        add_lines(data = vars$sgd_res$results, x = ~ date, y = ~ q0, yaxis = 'y1', name = 'Preferred SFGD', color = I('grey')) %>% 
        add_lines(data = vars$sgd_res$results, x = ~ date, y = ~ q0, yaxis = 'y1', name = 'Current SFGD', color = I('orange'))
    })
  })
  
  observeEvent(input$reset, {
    updateSliderInput(session, "bucket1_z", value = vars$cali_par$bucket1$z)
    updateSliderInput(session, "bucket1_rho_b", value = vars$cali_par$bucket1$rho_b)
    updateSliderInput(session, "bucket1_m", value = vars$cali_par$bucket1$m)
    updateSliderInput(session, "bucket1_rho_s", value = vars$cali_par$bucket1$rho_s)
    updateSliderInput(session, "bucket1_Ks", value = vars$cali_par$bucket1$Ks)
    updateSliderInput(session, "bucket1_n", value = vars$cali_par$bucket1$n)
    updateSliderInput(session, "bucket1_pE", value = vars$cali_par$bucket1$pE)
    updateSliderInput(session, "bucket2_z", value = vars$cali_par$bucket2$z)
    updateSliderInput(session, "bucket2_rho_b", value = vars$cali_par$bucket2$rho_b)
    updateSliderInput(session, "bucket2_m", value = vars$cali_par$bucket2$m)
    updateSliderInput(session, "bucket2_rho_s", value = vars$cali_par$bucket2$rho_s)
    updateSliderInput(session, "bucket2_Ks", value = vars$cali_par$bucket2$Ks)
    updateSliderInput(session, "bucket2_n", value = vars$cali_par$bucket2$n)
    updateSliderInput(session, "bucket2_pZ", value = vars$cali_par$bucket2$pZ)
    updateSliderInput(session, "curveNumber_CN2", value = vars$cali_par$cn$CN2)
    updateSliderInput(session, "curveNumber_pF", value = vars$cali_par$cn$pF)
    updateSliderInput(session, "aquifer_Td", value = vars$cali_par$aquifer$Td)
    updateSliderInput(session, "aquifer_rho_s", value = vars$cali_par$aquifer$rho_s)
    updateSliderInput(session, "aquifer_rho_f", value = vars$cali_par$aquifer$rho_f)
    updateSliderInput(session, "aquifer_Ka", value = vars$cali_par$aquifer$Ka)
    updateSliderInput(session, "aquifer_z0", value = vars$cali_par$aquifer$z0)
    updateSliderInput(session, "aquifer_Sy", value = vars$cali_par$aquifer$Sy)
    updateSliderInput(session, "aquifer_xT", value = vars$cali_par$aquifer$XT)
    updateSliderInput(session, "aquifer_dxT", value = vars$cali_par$aquifer$dxT_max)
    updateSliderInput(session, "windowSize", value = 120)
    # updateSliderInput(session, "Pumping", value = 50)
    
  })
  output$is_par_uploaded <- reactive({
    !is.null(vars$cali_par) & !is.null(vars$cnst_par)
  })
  outputOptions(output, "is_par_uploaded", suspendWhenHidden = FALSE)
  
  
  
  
  
  modelResults <- reactive({
    
    if (is.null(vars$bound_data) | is.null(vars$cali_par) | is.null(vars$cnst_par)) {
      return(NULL)
    }
    
    tryCatch({
      cali_params <- c(1, input$bucket1_z, input$bucket1_rho_b, input$bucket1_m, input$bucket1_rho_s, input$bucket1_Ks, input$bucket1_n, input$bucket1_pE,
                       2, input$bucket2_z, input$bucket2_rho_b, input$bucket2_m, input$bucket2_rho_s, input$bucket2_Ks, input$bucket2_n, input$bucket2_pZ,
                       input$curveNumber_CN2, input$curveNumber_pF, 
                       input$aquifer_Td, input$aquifer_rho_s, input$aquifer_rho_f, input$aquifer_Ka, input$aquifer_z0, input$aquifer_Sy, 
                       input$aquifer_xT, input$aquifer_dxT)
      cali_parlist <- parameter_vec_to_list(cali_params, prameter_skeleton = vars$cali_par)
      
      return(estimate_sgd(vars$bound_data, cali_parlist, vars$cnst_par, nw = input$windowSize)$results %>% 
               mutate(date = as.Date(input$start_date) - 1 + t))
    }, error = function(e) {
      shinyjs::alert("Error in model run. Check the parameter values.")
      NULL
    })
  })
  
 
  observe({plotlyProxy("wlPlot", session) %>%
      plotlyProxyInvoke( "restyle", list(
        y = list(modelResults()$hf),
        x = list(modelResults()$date)
      ), list(3)
      )
  })
  # 
  observe({plotlyProxy("q0Plot", session) %>%
      plotlyProxyInvoke( "restyle", list(
        y = list(modelResults()$q0),
        x = list(modelResults()$date)
      ), list(1)
      )
  })
}

NULL
#' GUI of SFGD calculator
#' 
#' @export
launch_gui <- function() {
  shinyApp(ui = ui, server = server)
}

