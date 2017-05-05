# Sources:
# Zoom: http://shiny.rstudio.com/gallery/plot-interaction-zoom.html
# Hover: http://stackoverflow.com/questions/38917101/how-do-i-show-the-y-value-on-tooltip-while-hover-in-ggplot2
options(shiny.maxRequestSize=120*1024^2) 
require(shiny)
require(data.table)
library(ggplot2)
library(Cairo)   # For nicer ggplot2 output when deployed on Linux
# setwd("methylation seq")
# dt1 <- fread("results.csv")
# dt1
# setkey(dt1, start)
# lvls <- unique(dt1$chr)[order(as.numeric(substr(unique(dt1$chr), 4, 6)))]
# dt1$chr <- factor(dt1$chr,
#                   levels = lvls)

ui <- fluidPage(
  fluidRow(
    style = "position:relative",
    fileInput(inputId = "fileIn", 
              label = "Upload File",
              accept = ".csv"),
    uiOutput("chr"),
    plotOutput("plot1", 
               height = 400,
               dblclick = "plot1_dblclick",
               hover = "plot_hover",
               brush = brushOpts(
                 id = "plot1_brush",
                 resetOnNew = TRUE)),
    uiOutput("hover_info")
  )
)

server <- function(input, output) {
  output$chr <- renderUI({
    validate(need(input$fileIn != "", ""))
    dt1 <- fread(input$fileIn$datapath, sep = ",")
    
    setkey(dt1, start)
    lvls <- unique(dt1$chr)[order(as.numeric(substr(unique(dt1$chr), 4, 6)))]
    dt1$chr <- factor(dt1$chr,
                      levels = lvls)
    
    selectInput(inputId = "chr",
                label = "Select Chromosome",
                choices = levels(dt1$chr))
  })
  
  ranges <- reactiveValues(x = NULL)
  
  output$plot1 <- renderPlot({
    validate(need(input$fileIn != "", ""))
    dt1 <- fread(input$fileIn$datapath, sep = ",")
    
    setkey(dt1, start)
    lvls <- unique(dt1$chr)[order(as.numeric(substr(unique(dt1$chr), 4, 6)))]
    dt1$chr <- factor(dt1$chr,
                      levels = lvls)
    
    tmp <- droplevels(subset(dt1, 
                             chr == input$chr))
    tmp$ystart <- rep(0, nrow(tmp))
    tmp$yend <- rep(1, nrow(tmp))
    tmp$red.meth <- "green"
    tmp$red.meth[tmp$`Control->Exptl:diff` < 0] <- "red"
    tmp$red.meth[tmp$`Control->Exptl:pval` > 0.01] <- "black"

    ggplot(tmp) +
      geom_segment(aes(x = start,
                       xend = start,
                       y = ystart,
                       yend = CpG),
                   col = tmp$red.meth) +
      scale_x_continuous("Starting Location of Regions") +
      scale_y_continuous("Total CpG") +
      ggtitle(paste("Chromosome:", tmp$chr[1])) +
      coord_cartesian(xlim = ranges$x)
  })
  
  output$hover_info <- renderUI({
    validate(need(input$fileIn != "", ""))
    dt1 <- fread(input$fileIn$datapath, sep = ",")
    
    setkey(dt1, start)
    lvls <- unique(dt1$chr)[order(as.numeric(substr(unique(dt1$chr), 4, 6)))]
    dt1$chr <- factor(dt1$chr,
                      levels = lvls)
    
    tmp <- droplevels(subset(dt1, 
                             chr == input$chr))
    tmp$ystart <- rep(0, nrow(tmp))
    tmp$yend <- rep(1, nrow(tmp))
    tmp$red.meth <- "green"
    tmp$red.meth[tmp$`Control->Exptl:diff` < 0] <- "red"
    tmp$red.meth[tmp$`Control->Exptl:pval` > 0.01] <- "black"
    
    hover <- input$plot_hover
    point <- nearPoints(tmp, 
                        hover, 
                        threshold = 5,
                        maxpoints = 1, 
                        addDist = TRUE)
    if (nrow(point) == 0) return(NULL)
    
    # calculate point position INSIDE the image as percent of total dimensions
    # from left (horizontal) and from top (vertical)
    left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
    top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
    
    # calculate distance from left and bottom side of the picture in pixels
    left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
    top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
    
    # create style property fot tooltip
    # background color is set so tooltip is a bit transparent
    # z-index is set so we are sure are tooltip will be on top
    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                    "left:", left_px + 2, "px; top:", top_px + 2, "px;")
    
    # actual tooltip created as wellPanel
    wellPanel(
      style = style,
      p(HTML(paste0("<b> Car: </b>", rownames(point), "<br/>",
                    "<b> mpg: </b>", point$mpg, "<br/>",
                    "<b> hp: </b>", point$hp, "<br/>",
                    "<b> Distance from left: </b>", left_px, "<b>, from top: </b>", top_px)))
    )
  })
  
  # When a double-click happens, check if there's a brush on the plot.
  # If so, zoom to the brush bounds; if not, reset the zoom.
  observeEvent(input$plot1_dblclick, {
    brush <- input$plot1_brush
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      
    } else {
      ranges$x <- NULL
    }
  })
}

shinyApp(ui, server)