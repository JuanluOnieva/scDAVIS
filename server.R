library(shiny)
library(shinyFiles)
library(Seurat)
library(dplyr)
library(RColorBrewer)
library(shinycssloaders)
library(Matrix)
library(DT)
library(shinyjs)

options(shiny.maxRequestSize = 500*1024^2)

source("function.R")
ExamplePath <- "examples"

shinyServer(function(input, output, session) {
  
  
  # Variable where save useful objects
  rv <- reactiveValues(sc=NULL) 
  
  
  # Create seuratSC and fill basic info
  observeEvent(input$loadSC, {
    if(input$radioLoadSC==1)
      sc <- paste(ExamplePath, "sc.rds", sep = "/")
    else if(input$radioLoadSC==2)
      sc <- paste(ExamplePath, "seuratSC.rds", sep = "/")
    else
      sc <- input$file$datapath
    withProgress({
      
    rv$sc <- readRDS(sc)
    
    rv$geneset <- rownames(rv$sc@raw.data)
    rv$identset <- colnames(rv$sc@meta.data)
    cat <- unlist(lapply(colnames(rv$sc@meta.data), function(x){ is.character(rv$sc@meta.data[, x]) | is.factor(rv$sc@meta.data[, x]) | is.logical(rv$sc@meta.data[, x]) }))
    rv$categorical <- colnames(rv$sc@meta.data)[cat]
    rv$numerical <- colnames(rv$sc@meta.data)[!cat]
    
    if(is.null(rv$gene))
      rv$gene <- rv$geneset[1]
    
    if(is.null(rv$ident))
      rv$ident <- ifelse("res.1" %in% rv$identset, "res.1", rv$identset[1])
    
    updateSelectInput(session = session, inputId = "selectIdent", choices = rv$identset, selected = rv$ident)
    updateSelectInput(session = session, inputId = "selectGene", choices = rv$geneset, selected = rv$gene)
    updateSelectInput(session = session, inputId = "selectCategorical", choices = rv$categorical, selected = rv$categorical[1])
    updateSelectInput(session = session, inputId = "selectNumerical", choices = rv$numerical, selected = rv$numerical[1])

    # Histogram about genes
    output$histGenes <- renderPlot({
      validate(
        need(input$loadSC, "Please select Seurat object"),
        need(rv$sc@meta.data$nGene, "Please select Seurat object with info about nGene")
      )
      ggplot() + geom_histogram(data = rv$sc@meta.data, aes(nGene), col="red", fill="blue")
    })
    
    # Histogram about UMIs
    output$histUMIs <- renderPlot({
      validate(
        need(input$loadSC, "Please select Seurat object"),
        need(rv$sc@meta.data$nUMI, "Please select Seurat object with info about nUMI")
      )
      ggplot() + geom_histogram(data = rv$sc@meta.data, aes(nUMI), col="red", fill="blue")
    })
   }, message = "Reading object")
  })
  
  # Info number of cells 
  output$infoNumberCells <- renderText({
    validate(
      need(input$loadSC, "Please select Seurat object")
    )
    ncol(rv$sc@data)
  })
  
  # Info number of reads 
  output$infoReadsPerCells <- renderText({
    validate(
      need(input$loadSC, "Please select Seurat object")
    )
    round(mean(colSums(rv$sc@raw.data)))
  })
  
  # Info number of genes
  output$infoNumberGenes <- renderText({
    validate(
      need(input$loadSC, "Please select Seurat object")
    )
    
    nrow(rv$sc@data)
  })
  
  # Info vln plots
  infoVln <- reactive(
    VlnPlot(object = rv$sc, features.plot = input$selectNumerical, 
            group.by = input$selectCategorical, do.return = T, point.size.use = input$checkboxPoint, y.lab.rot = T)
  )
  
  # Info vln plots
  output$InfoVln <- renderPlot({
    validate(
      need(input$loadSC, "Please select Seurat object")
    )
    infoVln()
  })
  
  # Get cluster plot
  clusterPlot <- reactive({
    nms <- rownames(rv$sc@dr[[input$radioRepresentation]]@cell.embeddings[, 1:2])
    col <- rv$sc@meta.data[nms, input$selectIdent]
    type <- rv$sc@meta.data[, input$selectIdent]
    if(is.character(type) |
       is.factor(type) |
       is.logical(type)
    )
      p <- getPlot(sc = rv$sc, rep = input$radioRepresentation, key = nms, col = col) + geom_point(aes(col=col))
    
    if(is.numeric(type))
      p <- getPlot(sc = rv$sc, rep = input$radioRepresentation, key = nms, col = col, scale = T) + geom_point(aes(col=col))
    return(p)
  }
  )
  
  
  # Render ClusterPlot for cluster tab
  output$ClusterPlot <- renderPlotly({
    validate(
      need(rv$sc, "Please select a seurat object"),
      need(input$radioRepresentation, "Please select a representation")
    )
    ggplotly(clusterPlot()) %>% layout(dragmode = "lasso")
    
  })
  
  # Fill cluster A with ident info
  output$selectA <- renderUI({
    validate(
      need(rv$sc, "Load seurat"),
      need(input$selectIdent, "Select ident")
    )
    selectInput("Ident1",  label = h5("Cluster A"), choices = levels(as.factor(rv$sc@meta.data[, input$selectIdent])))
  })
  
  # Fill cluster B with ident info
  output$selectB <- renderUI({
    validate(
      need(rv$sc, "Load seurat"),
      need(input$selectIdent, "Select ident")
    )
    selectInput("Ident2",  label = h5("Cluster B"), choices = levels(as.factor(rv$sc@meta.data[, input$selectIdent])))
  })
  
  # Fill cluster A with ident info
  output$select1VsALL <- renderUI({
    validate(
      need(rv$sc, "Load seurat"),
      need(input$selectIdent, "Select ident")
    )
    selectInput("Ident1VsRest",  label = h5("Cluster A"), choices = levels(as.factor(rv$sc@meta.data[, input$selectIdent])))
  })
  
  
  # FindMks for cluster tab one vs one
  observeEvent(input$actionMks, {
      withProgress({
        rv$ident <- input$selectIdent
        rv$sc <- SetIdent(rv$sc, ident.use = rv$sc@meta.data[, input$selectIdent])
        df <- runFindMks(sc = rv$sc, ident1 = input$Ident1, ident2 = input$Ident2)
        rv$df.cluster.all <- df[["all"]]
        rv$df.cluster <- df[["filter"]]
      })
    })
  
  # FindMks for cluster tab one vs all
  observeEvent(input$actionMks1, {
    withProgress({
      rv$ident <- input$selectIdent
      rv$sc <- SetIdent(rv$sc, ident.use = rv$sc@meta.data[, input$selectIdent])
      df <- runFindMks(sc = rv$sc, ident1 = input$Ident1VsRest)
      rv$df.cluster.all <- df[["all"]]
      rv$df.cluster <- df[["filter"]]
    })
  })
  
  
  # FindMks for cluster tab all vs all
  observeEvent(input$actionMksAll, {
    withProgress({
    rv$ident <- input$selectIdent
    rv$sc <- SetIdent(rv$sc, ident.use = rv$sc@meta.data[, input$selectIdent])
    df <- runFindAllMks(rv$sc)
    rv$df.cluster.all <- df[["all"]]
    rv$df.cluster <- df[["filter"]]
    })
  })
  
  
  # Fill table for mks cluster
  output$df.mks <- DT::renderDataTable({
    validate(
      need(rv$sc, "Please select a seurat object"),
      need(rv$df.cluster, "Please run mks")
    )
    DT::datatable(
      rv$df.cluster,
      rownames = F,
      options = list(pageLength = 10,lengthChange=FALSE)) %>% formatRound(columns = c(2:ncol(rv$df.cluster)), digits = 2) %>% formatRound(columns = 1, digits = 3) %>% formatStyle(columns = c(1:ncol(rv$df.cluster)), 'text-align' = 'center') 
  })
  
  # Identify a group cluster
  output$brushCluster <- renderText({
    d <- event_data("plotly_selected")
    if (!is.null(d)){
      rv$cellCluster <- d$key
      return(length(d$key))
    }
  })
  
  # Save cell1Cluster
  observeEvent(input$saveCell1Cluster, {
    rv$cell1Cluster <- unlist(rv$cellCluster)
  })
  
  # Save cell2Cluster
  observeEvent(input$saveCell2Cluster, {
    rv$cell2Cluster <- unlist(rv$cellCluster)
  })
  
  # Render cell1 cluster
  output$cell1Cluster <- renderText({
    length(rv$cell1Cluster)
  })
  
  # Render cell2 cluster
  output$cell2Cluster <- renderText({
    length(rv$cell2Cluster)
  })
  
  # Get new plot for cluster tab with selected cells
  plotNewClusterCluster <- reactive({
    validate(
      need(rv$sc, "Please select a seurat object"),
      need(rv$cell1Cluster, "Please select save cell1"),
      need(rv$cell2Cluster, "Please select save cell1")
    )
    df <- as.data.frame(rv$sc@dr[[input$radioRepresentation]]@cell.embeddings[, 1:2])
    normCounts <- rv$sc@scale.data
    df$val <- "Other"
    df[rv$cell1Cluster, "val"] = "cell1"
    df[rv$cell2Cluster, "val"] = "cell2"
    my.x <- colnames(df)[1]
    my.y <- colnames(df)[2]
    p <- ggplot(data=df, aes(x=df[,my.x], y=df[,my.y])) + geom_point(aes(col=df$val))+
      theme(line = element_blank(),
            text = element_blank(),
            title = element_blank())
    return(p)
  })
  
  # Render plot of new cluser in cluster tab
  observeEvent(input$showNewIdentCluster, {
    
    output$plotNewClusterCluster <- renderPlotly({  
      ggplotly(plotNewClusterCluster())
    })
  })
  
  # Run button to find markers for identified cells
  observeEvent(input$actionMksCellsCluster, {
    rv$ident <- input$selectIdent
    rv$gene <- input$selectGene
    rv$sc <- SetIdent(rv$sc, ident.use = rv$sc@meta.data[, input$selectIdent])
    withProgress({
      seuratSC.V2 <- SubsetData(rv$sc, cells.use = c(rv$cell1Cluster, rv$cell2Cluster))
      seuratSC.V2@meta.data[rv$cell1Cluster, "my.res"] = "cell1"
      seuratSC.V2@meta.data[rv$cell2Cluster, "my.res"] = "cell2"
      seuratSC.V2 <- SetIdent(seuratSC.V2, ident.use = seuratSC.V2@meta.data[, "my.res"])
      df <- runFindMks(sc = seuratSC.V2, ident1 = "cell1", ident2 = "cell2")
      rv$df.cluster <- df[["filter"]]
      rv$df.cluster.all <- df[["all"]]
      
    }, message = "Running mks")
  })
  
  # Save new Ident cluster
  observeEvent(input$saveIdentCluster, {
    withProgress({
      rv$ident <- input$selectIdent
      mix <- intersect(rv$cell1Cluster, rv$cell2Cluster)
      cell1 <- rv$cell1Cluster[!rv$cell1Cluster %in% mix]
      cell2 <- rv$cell2Cluster[!rv$cell2Cluster %in% mix]
      mix <- c(mix, rv$sc@cell.names[!rv$sc@cell.names %in% c(cell1, cell2)])
      if(input$checkboxKeepIdentCluster){
        rv$sc <- SetIdent(rv$sc, ident.use = rv$sc@meta.data[, input$nameClusteringIdentCluster])
        mix <- WhichCells(rv$sc, ident = "Mix") 
        mix <- mix[!mix %in% c(cell1, cell2)]
      }
      rv$sc@meta.data[cell1, input$nameClusteringIdentCluster] = input$nameClusteringCell1Cluster
      rv$sc@meta.data[cell2, input$nameClusteringIdentCluster] = input$nameClusteringCell2Cluster
      rv$sc@meta.data[mix, input$nameClusteringIdentCluster] = "Mix"
      rv$identset <- colnames(rv$sc@meta.data)
      updateSelectInput(session = session, inputId = "selectIdent", choices = rv$identset, selected = rv$ident)
      
    })
  })
  
  # Plot gene for tab cells
  plotGene <- reactive({
    nms <- rownames(rv$sc@dr[[input$radioRepresentation]]@cell.embeddings[, 1:2])
    normCounts <- rv$sc@scale.data
    val <- input$selectGene
    if(length(unlist(val))>1){
      col <- as.factor(color.gradient(as.numeric(colSums(normCounts[val,colnames(rv$sc@scale.data)])),dscale = quantile(x = as.numeric(colSums(normCounts[val,colnames(rv$sc@scale.data)])),probs = c(input$sliderScale[1]/100,input$sliderScale[2]/100))))
    }else
      col <- as.factor(color.gradient(as.numeric(normCounts[val,colnames(rv$sc@scale.data)]),dscale = quantile(x = as.numeric(normCounts[val,colnames(rv$sc@scale.data)]),probs = c(input$sliderScale[1]/100,input$sliderScale[2]/100))))
    
    p <- getPlot(sc = rv$sc, rep = input$radioRepresentation, col = col, key = nms) + geom_point(col=col)
    return(p)
  })
  
  
  # Render gene plot for cell tab
  output$plotGene <- renderPlotly({
    validate(
      need(rv$sc, "Please select a seurat object"),
      need(input$selectGene, "Please select a gene")
    )
    ggplotly(plotGene()) %>% layout(dragmode = "lasso")
  })
  
  # Identify a group cell
  output$brush <- renderText({
    d <- event_data("plotly_selected")
    if (!is.null(d)){
      rv$cell <- d$key
      return(length(d$key))
    }
  })
  
  # Save cell1
  observeEvent(input$saveCell1, {
    rv$cell1 <- unlist(rv$cell)
  })
  
  # Save cell2
  observeEvent(input$saveCell2, {
    rv$cell2 <- unlist(rv$cell)
  })
  

  # Render cell1
  output$cell1 <- renderText({
    length(rv$cell1)
  })
  
  # Render cell2
  output$cell2 <- renderText({
    length(rv$cell2)
  })
  
  # Get new plot for gene/cell tab with selected cells
  plotNewClusterCell <- reactive({
    validate(
      need(rv$sc, "Please select a seurat object"),
      need(rv$cell1, "Please select save cell1"),
      need(rv$cell2, "Please select save cell1")
    )    
    df <- as.data.frame(rv$sc@dr[[input$radioRepresentation]]@cell.embeddings[, 1:2])
    df$val <- "Other"
    df[rv$cell1, "val"] = "cell1"
    df[rv$cell2, "val"] = "cell2"
    my.x <- colnames(df)[1]
    my.y <- colnames(df)[2]
    p <- ggplot(data=df, aes(x=df[,my.x], y=df[,my.y])) + geom_point(aes(col=df$val))+
      theme(line = element_blank(),
            text = element_blank(),
            title = element_blank())
    return(p)
  })
  
  # Render plot of new cluser in cells tab
  observeEvent(input$showNewIdentCell, {
    
    output$plotNewClusterCell <- renderPlotly({
      ggplotly(plotNewClusterCell())
    })
  })
  
  # Run button to find markers for identified cells
  observeEvent(input$actionMksCells, {
    rv$ident <- input$selectIdent
    rv$gene <- input$selectGene
    rv$sc <- SetIdent(rv$sc, ident.use = rv$sc@meta.data[, input$selectIdent])
    withProgress({
      seuratSC.V2 <- SubsetData(rv$sc, cells.use = c(rv$cell1, rv$cell2))
      seuratSC.V2@meta.data[rv$cell1, "my.res"] = "cell1"
      seuratSC.V2@meta.data[rv$cell2, "my.res"] = "cell2"
      seuratSC.V2 <- SetIdent(seuratSC.V2, ident.use = seuratSC.V2@meta.data[, "my.res"])
      df <- runFindMks(sc = seuratSC.V2, ident1 = "cell1", ident2 = "cell2")
      rv$df.cells <- df[["filter"]]
      rv$df.cells.all <- df[["all"]]
    
    }, message = "Running mks")
  })
  
  # Fill mks table for costum cluster
  output$df.mks.cell <- DT::renderDataTable({
    validate(
        need(rv$sc, "Please select a seurat object"),
        need(rv$df.cells, "Please run mks")
    )
    
    DT::datatable(
      rv$df.cells,
      rownames = F,
      options = list(pageLength = 10,lengthChange=FALSE)) %>% formatRound(columns = c(2:ncol(rv$df.cells)), digits = 2) %>% formatRound(columns = 1, digits = 3) %>% formatStyle(columns = c(1:ncol(rv$df.cells)), 'text-align' = 'center') 
  })
  
  # Save new Ident
  observeEvent(input$saveIdent, {
    rv$gene <- input$selectGene
    mix <- intersect(rv$cell1, rv$cell2)
    cell1 <- rv$cell1[!rv$cell1 %in% mix]
    cell2 <- rv$cell2[!rv$cell2 %in% mix]
    mix <- c(mix, rv$sc@cell.names[!rv$sc@cell.names %in% c(cell1, cell2)])
    if(input$checkboxKeepIdent){
      rv$sc <- SetIdent(rv$sc, ident.use = rv$sc@meta.data[, input$nameClusteringIdent])
      mix <- WhichCells(rv$sc, ident = "Mix") 
      mix <- mix[!mix %in% c(cell1, cell2)]
    }
    rv$sc@meta.data[cell1, input$nameClusteringIdent] = input$nameClusteringCell1
    rv$sc@meta.data[cell2, input$nameClusteringIdent] = input$nameClusteringCell2
    rv$sc@meta.data[mix, input$nameClusteringIdent] = "Mix"
    rv$identset <- colnames(rv$sc@meta.data)
    updateSelectInput(session = session, inputId = "selectIdent", choices = rv$identset, selected = rv$ident)
    
  })
  

  
  # Button for download mks in cluster tab
  output$downloadDataMks <- downloadHandler(
    filename = function() {
      paste(paste(gsub("[.]", "_", input$selectIdent), 
                  "mks", sep="_"), ".txt", sep = "")
    },
    content = function(file) {
      write.table(rv$df.cluster.all, file, col.names = TRUE, quote = F, row.names = F)
    }
  )
  
  # Button for download plot in cluster tab
  output$downloadPlotCluster <- downloadHandler(
    filename = function() {
      paste(paste(gsub("[.]", "_", input$selectIdent), 
                  "plot", sep="_"), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file)
      plot(clusterPlot())
      dev.off()
    }
  )
  
  # Button for download mks in cells tab
  output$downloadDataMksCells <- downloadHandler(
    filename = function() {
      paste(paste(gsub("[.]", "_", input$saveIdent), 
                  "mks", "cells", sep="_"), ".txt", sep = "")
    },
    content = function(file) {
      write.table(rv$df.cells.all, file, col.names = TRUE, quote = F, row.names = F)
    }
  )
  
  # Button for download plot in cells tab
  output$downloadPlotCells <- downloadHandler(
    filename = function() {
      paste(paste(gsub("[.]", "_", input$saveIdent), 
                  "plot", "cells", sep="_"), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file)
      plot(rv$plot.cells)
      dev.off()
    }
  )
  
  # Button for download plot of new cluster in tab cells
  output$downloadPlotGene <- downloadHandler(
    filename = function() {
      paste(paste(gsub("[.]", "_", input$selectGene), 
                  "plot", "cells", sep="_"), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file)
      plot(plotGene())
      dev.off()
    }
  )
  
  
  # Vln plot
  vlnPlot <- reactive({
    validate(
      need(rv$sc, "Please upload seurat object"),
      need(input$selectIdent, "Please select an ident")
    )
    
    type <- rv$sc@meta.data[, input$selectIdent]
    if(is.character(type) | is.factor(type) | is.logical(type))
      rv$doVln <- T
    else
      rv$doVln <- NULL
    validate({
      need(rv$sc, "Please select seutat Object")
      need(rv$doVln, "Please the ident must be categorical")
    }
    )
    point <- ifelse(input$checkboxShowPoint, 1, 0)
    plot <- VlnPlot(object = rv$sc, features.plot = input$selectGene, 
                    group.by = input$selectIdent, do.return = T, point.size.use = point, y.lab.rot = T)
    return(plot)
  })
  
  # Get vln plot
  output$plotVln <- renderPlot(vlnPlot())
  
  # Button for download vln plots
  output$downloadVlnPlot <- downloadHandler(
    filename = function() {
      paste(paste("vlnPlot", input$selectGene, input$selectIdent, sep = "_"), "pdf", sep = ".")
    },
    content = function(file) {
      pdf(file)
      plot(vlnPlot())
      dev.off()
    }
  )
  
  
  # Run mks of heatMap
  hmap.top <- reactive({
    rv$sc <- SetIdent(rv$sc, ident.use = rv$sc@meta.data[, input$selectIdent])
    withProgress(df <- FindAllMarkers(rv$sc), message = "Running mks")
    df <- df %>% group_by(cluster) %>% top_n(min(nrow(df), input$sliderNumMksHmap), avg_logFC)
    return(df$gene)
  })
  
  # Which gene use in heatMap
  geneSelected <- reactive({
    if(input$radioHMap==0)
      gene <- input$selectGene
    else
      gene <- hmap.top()
    return(gene)
  })
  
  # Action button heatMap
  observeEvent(input$actionMksHMap, {
    rv$ident <- input$selectIdent
    rv$gene <- input$selectGene
    rv$HMPlot <- DoHeatmap(rv$sc, genes.use = geneSelected(), group.by = input$selectIdent, group.label.rot = T, do.plot = F,
                           remove.key = F, slim.col.label = T)
  })
  
  output$HeatMapPlot <- renderPlot({
    validate(
      need(rv$HMPlot, "Need to select an option")
    )
    rv$HMPlot
  })
  
  # Button for download heatMap plots
  output$downloadHMap <- downloadHandler(
    filename = function() {
      paste(paste("HeatMap", input$selectIdent, sep = "_"), "pdf", sep = ".")
    },
    content = function(file) {
      pdf(file)
      plot(rv$HMPlot)
      dev.off()
    }
  )
  
  # Button for download vln info plots
  output$downloadVlnInfo <- downloadHandler(
    filename = function() {
      paste(paste("VlnInfo", input$selectNumerical, input$selectCategorical, sep = "_"), "pdf", sep = ".")
    },
    content = function(file) {
      pdf(file)
      plot(infoVln())
      dev.off()
    }
  )
  
  # Button for download new plot cells plots
  output$downloadPlotCells <- downloadHandler(
    filename = function() {
      paste(paste("newCluster", input$selectNumerical, input$selectCategorical, sep = "_"), "pdf", sep = ".")
    },
    content = function(file) {
      pdf(file)
      plot(plotNewClusterCell())
      dev.off()
    }
  )
  
  # Button for download new plot cluster plots
  output$downloadPlotNewCluster <- downloadHandler(
    filename = function() {
      paste(paste("newCluster", input$selectNumerical, input$selectCategorical, sep = "_"), "pdf", sep = ".")
    },
    content = function(file) {
       pdf(file)
       plot(plotNewClusterCluster())
       dev.off()
    }
  )
  
  # Button for download seurat object
  output$downloadSC <- downloadHandler(

    filename = function() {
      paste(input$nameSC, "rds", sep = ".")
    },
    content = function(file) {
      withProgress(
      saveRDS(rv$sc, file))
    }
  )
  
  
})