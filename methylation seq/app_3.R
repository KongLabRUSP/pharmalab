# Sources:
# Zoom: http://shiny.rstudio.com/gallery/plot-interaction-zoom.html
# Hover: http://stackoverflow.com/questions/38917101/how-do-i-show-the-y-value-on-tooltip-while-hover-in-ggplot2
options(shiny.maxRequestSize=100*1024^2) 
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

genes <- fread("data/genes.csv", header = FALSE)

dt.g <- subset(genes, select = c(1, 3:5, 9))
names(dt.g) <- c("chr", "reg", "start", "stop", "gene")

setkey(dt.g, start)
lvls <- unique(dt.g$chr)[order(as.numeric(substr(unique(dt.g$chr), 4, 6)))]
dt.g$chr <- factor(dt.g$chr,
                   levels = lvls)

dt.g$reg <- factor(dt.g$reg)

dt.g$gene <- substr(dt.g$gene, 11, nchar(dt.g$gene) - 2)
dt.g$gene <- factor(dt.g$gene, levels = unique(dt.g$gene))

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
               # height = 800,
               dblclick = "plot1_dblclick",
               hover = "plot_hover",
               brush = brushOpts(id = "plot1_brush",
                                 resetOnNew = TRUE)),
    plotOutput("plot2", 
               height = 100)
  )
)

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
  
# Plot1----
  output$plot1 <- renderPlot({
    tmp <- dd()
    
    if(input$gene == "All") {
      p1 <- ggplot(tmp) +
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
    } else {
      gn <- subset(dt.g,
                   chr == input$chr &
                     gene == input$gene)
      p1 <- ggplot(tmp) +
        facet_wrap(~ sample,
                   ncol = 1) +
        geom_rect(aes(xmin = start,
                      xmax = end,
                      ymin = ystart,
                      ymax = pct.meth),
                  col = "green") +
        ggtitle(paste("Chromosome:", tmp$chr[1], ", Percent CpG Methylated")) +
        scale_x_continuous("Starting Location of Regions",
                           limits = c(min(gn$start, na.rm = TRUE),
                                      max(gn$stop, na.rm = TRUE)))
    }
    p1
  })
  
# Plot2----
  output$plot2 <- renderPlot({
    gn <- droplevels(subset(dt.g,
                            chr == input$chr))
    if(input$gene == "All") {
      p2 <- ggplot(gn) +
        geom_rect(aes(xmin = start,
                      xmax = stop,
                      ymin = rep(0, length(start)),
                      ymax = rep(1, length(start)),
                      fill = reg)) +
        coord_cartesian(xlim = ranges$x) +
        theme(legend.position = "bottom") +
        guides(fill = guide_legend(title = "Region",
                                   nrow = 1))
    } else {
      gn <- subset(gn,
                   gene == input$gene)
      p2 <- ggplot(gn) +
        geom_rect(aes(xmin = start,
                      xmax = stop,
                      ymin = rep(0, length(start)),
                      ymax = rep(1, length(start)),
                      fill = reg)) +
        scale_x_continuous(limits = c(min(gn$start, na.rm = TRUE),
                                      max(gn$stop, na.rm = TRUE))) +
        theme(legend.position = "bottom") +
        guides(fill = guide_legend(title = "Region",
                                   nrow = 1))
    }
    p2
  })

  # output$hover_info <- renderUI({
  #   tmp <- dd()
  #   
  #   hover <- input$plot_hover
  #   point <- nearPoints(tmp,
  #                       hover,
  #                       threshold = 5,
  #                       maxpoints = 1,
  #                       addDist = TRUE)
  #   if (nrow(point) == 0) return(NULL)
  #   
  #   # calculate point position INSIDE the image as percent of total dimensions
  #   # from left (horizontal) and from top (vertical)
  #   left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
  #   top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
  #   
  #   # calculate distance from left and bottom side of the picture in pixels
  #   left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
  #   top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
  #   
  #   # create style property fot tooltip
  #   # background color is set so tooltip is a bit transparent
  #   # z-index is set so we are sure are tooltip will be on top
  #   style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
  #                   "left:", left_px + 2, "px; top:", top_px + 2, "px;")
  #   
  #   # actual tooltip created as wellPanel
  #   wellPanel(
  #     style = style,
  #     p(HTML(paste0("<b> Car: </b>", rownames(point), "<br/>",
  #                   "<b> mpg: </b>", point$mpg, "<br/>",
  #                   "<b> hp: </b>", point$hp, "<br/>",
  #                   "<b> Distance from left: </b>", left_px, "<b>, from top: </b>", top_px))))
  # })

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