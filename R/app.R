library(shiny)
library(ggplot2)
library(ggimage)
library(grid)
library(gridExtra)

options(shiny.maxRequestSize = 30*1024^2,
        shiny.launch.browser = .rs.invokeShinyWindowExternal)

ui <- fluidPage(
  fluidRow(column(width = 2,
                  h4("1. Choose exported file to upload"),
                  fileInput("upload", NULL),
                  h4("2. Choose year to display"),
                  sliderInput("year",
                              "",
                              min = 1,
                              max = 1,
                              value = 1,
                              step = 1
                  )
  ),
  column(width = 5,
         plotOutput("plot", height = "800px",)
  )
  )
)

server = function(input, output, session){
  
  session$onSessionEnded(function() {
    stopApp()
  })
  
  ## 1. load up image set ##
  img_file <- reactive({
    req(input$upload)
    path <- input$upload$datapath
    img_file <- get(load(file = path))
    img_file
  })
  
  img_len <- reactive({
    req(img_file())
    img_len <- length(img_file())
    img_len
  })
  
  observe({
    req(img_len())
    updateSliderInput(session, "year", max = img_len())
  })
  
  ## 2. Create a plot ##
  observe({
    output$plot <- renderPlot({
      req(img_file())
      grid.draw(img_file()[[input$year]])
    })
  })
}

shinyApp(ui, server)