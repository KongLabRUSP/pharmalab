# Sources:
# Zoom: http://shiny.rstudio.com/gallery/plot-interaction-zoom.html
# Hover: http://stackoverflow.com/questions/38917101/how-do-i-show-the-y-value-on-tooltip-while-hover-in-ggplot2
###################################################
# Header----
options(shiny.maxRequestSize=200*1024^2) 
require(shiny)
require(data.table)
library(ggplot2)
library(Cairo)

# Load genes----
# C:/git_local/pharmalab/methylation seq/source/
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
lgenes <- c("<Select a gene>", levels(dt.g$gene)[order(levels(dt.g$gene))])

# UI----
ui <- fluidPage(
  fluidRow(
    column(5,
           style = "position:relative",
           fileInput(inputId = "fileIn", 
                     label = "Upload File",
                     accept = ".csv"),
           uiOutput("chr"),
           uiOutput("gene")
    )#,
    # column(5,
    #        selectizeInput(inputId = "allgenes",
    #                       label = "Gene Lookup",
    #                       choices = lgenes,
    #                       options = list(maxOptions = 10^6)),
    #        verbatimTextOutput(outputId = "selected")
    # )
  ),
  fluidRow(
    plotOutput("plot1", 
               height = 400,
               dblclick = "plot1_dblclick",
               brush = brushOpts(
                 id = "plot1_brush",
                 resetOnNew = TRUE)),
    plotOutput("plot2", 
               height = 100)
  )
  # fluidRow(
  #   dataTableOutput("dt1")
  # )
)

# Server----
server <- function(input, output) {
  # Select chromosome----
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
  
  # Select gene----  
  output$gene <- renderUI({
    gn <- droplevels(subset(dt.g,
                            chr == input$chr))
    
    selectizeInput(inputId = "gene",
                   label = "Select Gene",
                   choices = c("All", levels(gn$gene)[order(levels(gn$gene))]),
                   selected = "All",
                   options = list(maxOptions = 10^6))
  })
  
  ranges <- reactiveValues(x = NULL)
  
  # Plot1: CpG DMR----
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
    tmp$red.meth[tmp$`Control->Exptl:pval` > 0.05] <- "black"
    
    ggplot(tmp) +
      geom_rect(aes(xmin = start,
                    xmax = end,
                    ymin = ystart,
                    ymax = CpG),
                fill = tmp$red.meth,
                col = tmp$red.meth) +
      scale_x_continuous("Starting Location of Regions") +
      ggtitle(paste("Total CpG on Chromosome:", tmp$chr[1])) +
      coord_cartesian(xlim = ranges$x)
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
                    fill = reg,
                    alpha = 0.5)) +
      coord_cartesian(xlim = ranges$x) +
      theme(legend.position = "bottom") +
      guides(fill = guide_legend(title = "Region",
                                 nrow = 1))
  })
  
  # Zoom on click----
  observeEvent(input$plot1_dblclick, {
    brush <- input$plot1_brush
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      
    } else {
      ranges$x <- NULL
    }
  })
  
  # When a gene selected, override x-limit to the range of that gene----
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
  
  # # In gene lookup, when gene is selected, show chromosome and location
  # observeEvent(input$allgenes, {
  #   if(input$allgenes != "<Select a gene>") {
  #     output$selected <- renderText(paste("Gene Start:",
  #                                         as.character(min(dt.g$start[dt.g$gene == input$allgenes])),
  #                                         "\n",
  #                                         "Gene End:",
  #                                         as.character(max(dt.g$stop[dt.g$gene == input$allgenes])),
  #                                         "\n",
  #                                         "Chromosome:",
  #                                         dt.g$chr[dt.g$gene == input$allgenes][1])) 
  #   }
  # })
  
  output$dt1 <- renderDataTable({
    
  })
}

shinyApp(ui, server)