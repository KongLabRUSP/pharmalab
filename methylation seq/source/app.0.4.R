# Project: MethylSeq Data Visualization
# Author: Davit Sargsyan
# Created: 03/21/2017
# Sources:
# Zoom: http://shiny.rstudio.com/gallery/plot-interaction-zoom.html
# Hover: http://stackoverflow.com/questions/38917101/how-do-i-show-the-y-value-on-tooltip-while-hover-in-ggplot2

# Header----
options(shiny.maxRequestSize=100*1024^2) 
require(shiny)
require(data.table)
library(ggplot2)
# For nicer ggplot2 output when deployed on Linux
library(Cairo)

# Load genes----
genes <- fread("genes.csv", header = FALSE)

dt.g <- subset(genes, select = c(1, 3:5, 9))
names(dt.g) <- c("chr", "reg", "start", "stop", "gene")

setkey(dt.g, start)
lvls <- unique(dt.g$chr)[order(as.numeric(substr(unique(dt.g$chr), 4, 6)))]
dt.g$chr <- factor(dt.g$chr,
                   levels = lvls)

dt.g$reg <- factor(dt.g$reg)

dt.g$gene <- substr(dt.g$gene, 11, nchar(dt.g$gene) - 2)
dt.g$gene <- factor(dt.g$gene, levels = unique(dt.g$gene))

# UI----
ui <- fluidPage(
  fluidRow(
    style = "position:relative",
    fileInput(inputId = "fileIn", 
              label = "Upload File",
              accept = ".csv"),
    uiOutput("chr"),
    uiOutput("sample"),
    uiOutput("gene"),
    plotOutput("plot1", 
               dblclick = "plot1_dblclick",
               brush = brushOpts(id = "plot1_brush",
                                 resetOnNew = TRUE)),
    plotOutput("plot2", 
               height = 100)
  )
)

# Server----
server <- function(input, output) {
  dd <- reactive({
    validate(need(input$fileIn != "", ""))
    dt1 <- fread(input$fileIn$datapath)
    
    setkey(dt1, start)
    lvls <- unique(dt1$chr)[order(as.numeric(substr(unique(dt1$chr), 4, 6)))]
    dt1$chr <- factor(dt1$chr,
                      levels = lvls)
    
    tmp <- droplevels(subset(dt1,
                             chr == input$chr))
    
    tmp$ystart <- rep(0, nrow(tmp))
    tmp$yend <- rep(1, nrow(tmp))
    
    tmp <- melt.data.table(tmp,
                           id.vars = c("chr", "start", "end", "CpG", "ystart", "yend"),
                           variable.name = "experiment",
                           value.name = "nreadns")
    tmp$experiment <- as.character(tmp$experiment)
    
    tmp$sample <- substr(tmp$experiment, 
                         1, 
                         nchar(tmp$experiment) - 2)
    tmp$sample <- factor(tmp$sample,
                         levels = unique(tmp$sample))
    
    if (!is.null(input$sample)) {
      tmp <- droplevels(tmp[sample %in% input$sample, ])
    }
    
    tmp$what <- substr(tmp$experiment, nchar(tmp$experiment), nchar(tmp$experiment))
    tmp$what <- factor(tmp$what,
                       levels = c("X", "N"))
    tmp.X <- subset(tmp, what == "X", select = -c(7, 10))
    tmp.N <- subset(tmp, what == "N", select = -c(7, 10))
    tmp <- merge(tmp.X, tmp.N,
                 by = c("chr", "start", "end", "CpG", "ystart", "yend", "sample"))
    tmp$pct.meth <- 100*tmp$nreadns.x/tmp$nreadns.y
    tmp
  })
  # Select chromosome----  
  output$chr <- renderUI({
    validate(need(input$fileIn != "", ""))
    dt1 <- fread(input$fileIn$datapath)
    
    setkey(dt1, start)
    lvls <- unique(dt1$chr)[order(as.numeric(substr(unique(dt1$chr), 4, 6)))]
    dt1$chr <- factor(dt1$chr,
                      levels = lvls)
    
    selectInput(inputId = "chr",
                label = "Select Chromosome",
                choices = levels(dt1$chr))
  })
  # Select sample----  
  output$sample <- renderUI({
    validate(need(input$fileIn != "", ""))
    dt1 <- fread(input$fileIn$datapath)
    
    setkey(dt1, start)
    lvls <- unique(dt1$chr)[order(as.numeric(substr(unique(dt1$chr), 4, 6)))]
    dt1$chr <- factor(dt1$chr,
                      levels = lvls)
    
    tmp <- droplevels(subset(dt1,
                             chr == input$chr))
    
    tmp <- melt.data.table(tmp,
                           id.vars = c("chr", "start", "end", "CpG"),
                           variable.name = "experiment",
                           value.name = "nreadns")
    tmp$experiment <- as.character(tmp$experiment)
    
    tmp$sample <- substr(tmp$experiment, 1, nchar(tmp$experiment) - 2)
    tmp$sample <- factor(tmp$sample,
                         levels = unique(tmp$sample))
    
    selectInput(inputId = "sample",
                label = "Select Sample(s)",
                choices = levels(tmp$sample),
                selected = NULL,
                multiple = TRUE)
  })
  # Select gene----  
  output$gene <- renderUI({
    gn <- droplevels(subset(dt.g,
                            chr == input$chr))
    
    selectInput(inputId = "gene",
                label = "Select Gene",
                choices = c("All", levels(gn$gene)[order(levels(gn$gene))]),
                selected = "All")
  })
  
  # Plot area selection----
  ranges <- reactiveValues(x = NULL)
  
  # Plot1: % methylated----
  output$plot1 <- renderPlot({
    tmp <- dd()
    
    ggplot(tmp) +
      facet_wrap(~ sample,
                 ncol = 1) +
      geom_rect(aes(xmin = start,
                    xmax = end,
                    ymin = ystart,
                    ymax = pct.meth),
                col = "green") +
      ggtitle(paste("Chromosome:", tmp$chr[1], ", Percent CpG Methylated")) +
      coord_cartesian(xlim = ranges$x) +
      scale_x_continuous("Starting Location of Regions")
  })
  
  # Plot2: gene regions----
  output$plot2 <- renderPlot({
    gn <- droplevels(subset(dt.g,
                            chr == input$chr))
    ggplot(gn) +
      geom_rect(aes(xmin = start,
                    xmax = stop,
                    ymin = rep(0, length(start)),
                    ymax = rep(1, length(start)),
                    fill = reg)) +
      coord_cartesian(xlim = ranges$x) +
      theme(legend.position = "bottom") +
      guides(fill = guide_legend(title = "Region",
                                 nrow = 1))
  })
  
  # Observe events----
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
  
  # When a gene selected, override x-limit to the range of that gene
  observeEvent(input$gene, {
    if(input$gene != "All") {
      gn <- subset(dt.g,
                   chr == input$chr &
                     gene == input$gene)
      
      ranges$x <- c(min(gn$start, na.rm = TRUE),
                    max(gn$stop, na.rm = TRUE))
    } else {
      ranges$x <- NULL
    }
  })
}

shinyApp(ui, server)