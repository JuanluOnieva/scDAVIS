library(shiny)
library(plotly)
library(shinyWidgets)
library(shinydashboard)
library(shinyFiles)
library(shinyjs)
library(shinycssloaders)

shinyUI(dashboardPage(
  dashboardHeader(
    title = "scDAVIS: single-cell Data Analysis VISualization"
    , titleWidth = 350
    ),
  dashboardSidebar(
    width = 350,
    collapsed = T,
    sidebarMenu(
      menuItem("Load seurat Object"
               , tabName = "info"
               , fileInput(inputId = "file", label = "Chose directory", accept = 'rds')
               , hr()
               , radioButtons("radioLoadSC", label = h3("Select seurat object")
                              , choices = list("My sc" = 0, "Mph example" = 1, "All cell example" = 2)
                              , selected = 0)
               , actionButton(
                 inputId = "loadSC",
                 label = "Load SC object"
               )
      )
      , menuItem("Visualization"
                 , selectInput("selectIdent", label = h3("Select clustering")
                               , choices = NULL, selected = NULL)
                 , radioButtons("radioRepresentation", label = h3("Representation")
                                , choices = list("TSNE" = "tsne", "PCA" = "pca")
                                , selected = "tsne")
                 , hr()
      
      )
      , menuItem( "Map Gene"
                  , selectInput("selectGene", label = h3("Select gene")
                                , choices = NULL, selected = NULL, multiple = T)
      )
      , menuItem( "Download seurat object"
                  , fluidRow(
                    column(9
                           , textInput("nameSC", label = "Name of sc", value = "scAnalysis")
                                       
                    ),
                    column(3, offset = 9
                           , downloadButton("downloadSC", "")
                    )
                  )
      )
    )
  )
  , dashboardBody(
    navbarPage("",
               tabPanel("Info",
                        fluidPage(
                          fluidRow(
                            box(title="Basic info"
                                , valueBox(subtitle = "Number of cells",value = textOutput("infoNumberCells"), width = 6)
                                , valueBox(subtitle = "Number of reads per cells",value = textOutput("infoReadsPerCells"), width = 6)
                                , valueBox(subtitle = "Number of genes",value = textOutput("infoNumberGenes"), width = 6)
                                , collapsible = T,  collapsed = T
                                , width = 12
                            )
                            # , box(title="Cell type info"
                            #       , tableOutput("infoCelltype")
                            #       , plotOutput("infoCelltypePlot") %>% withSpinner()
                            #       , collapsible = T,  collapsed = T
                            #       , width = 12
                            # )
                            , box(title="Distribution of genes"
                                  , plotOutput("histGenes") %>% withSpinner()
                                  , collapsible = T, collapsed = T
                                  , width = 12
                            )
                            , box(title="Distribution of reads"
                                  , plotOutput("histUMIs") %>% withSpinner()
                                  , collapsible = T, collapsed = T
                                  , width = 12
                            )
                            , box(title="QC and metadata info"
                                  , fluidRow(
                                    column(4
                                           , selectInput("selectCategorical", label = h3("Select categorical")
                                                         , choices = NULL, selected = NULL)                                           
                                    ),
                                    column(4
                                           , selectInput("selectNumerical", label = h3("Select numerical")
                                                         , choices = NULL, selected = NULL, multiple = T)
                                    ),
                                    column(2
                                           , br()
                                    ),
                                    column(2
                                           ,downloadButton("downloadVlnInfo", "")  
                                           )
                                  )
                                  
                                 , plotOutput("InfoVln") %>% withSpinner()
                                 , fluidRow(
                                   column(2, offset = 10
                                          , dropdownButton(
                                            checkboxInput("checkboxPoint", label = "Show point", value = F)
                                            , icon = icon(name = "cog")
                                          )
                                                                                 
                                   )
                                 )
                                  , collapsible = T, collapsed = T
                                  , width = 12
                            )
                          )
                        )
               )
               , tabPanel("Cluster",
                          fluidRow(
                            tabBox(title = ""
                                   , width = 8
                                   , tabPanel("Graph"
                                              , fluidRow(
                                                column(2, offset = 10
                                                       , downloadButton("downloadPlotCluster", "")
                                                )
                                              )
                                              , plotlyOutput("ClusterPlot") %>% withSpinner()

                                    )
                                   , tabPanel("new Graph"
                                              , fluidRow(
                                                column(2, offset = 10
                                                       , downloadButton("downloadPlotNewCluster", "")
                                                )
                                              )
                                              , plotlyOutput("plotNewClusterCluster") %>% withSpinner()
                                              
                                    )      
                                   ,  tabPanel('Table'
                                               , fluidRow(
                                                 column(2, offset = 10
                                                        , downloadButton("downloadDataMks", "")
                                                 )
                                               )
                                               ,br()
                                               , DT::dataTableOutput('df.mks')
                                               
                                    )
                                   ),
                            tabBox(
                              side = "right", width = 4,
                              selected = "Select cluster",
                              tabPanel("Select cluster", 
                                       h3("Two condition contrast")
                                       , uiOutput("selectA")
                                       , uiOutput("selectB")
                                       # , textInput("Ident1", label = h5("Cluster A"), value = "Enter text...")
                                       # , textInput("Ident2", label = h5("Cluster B"), value = "Enter text...")
                                       , actionButton("actionMks", label = "Run mks")
                                       , hr()
                                       , h3("One vs all contrast")
                                       , uiOutput("select1VsALL")
                                       # , selectInput("Ident1VsRest", label = h5("Cluster A")
                                       #                , choices = NULL, selected = NULL)
                                       #, textInput("Ident1VsRest", label = h5("Cluster A"), value = "Enter text...")
                                       , actionButton("actionMks1", label = "Run mks")
                                       , hr()
                                       , h3("All vs all contrast")
                                       , actionButton("actionMksAll", label = "Run mks")
                                       )
                              ,tabPanel("Costum cluster",
                                       h4("Selected cells")
                                       , textOutput("brushCluster")
                                       , hr()
                                       , fluidRow(
                                         column(4
                                                , actionButton("saveCell1Cluster", label = "Save Cell1")

                                         ),
                                         column(4
                                                , actionButton("saveCell2Cluster", label = "Save Cell2")
                                         )
                                       )
                                       , br()
                                       , fluidRow(
                                         column(6
                                                , textOutput("cell1Cluster")
                                                #valueBox(subtitle = h4(""),value = h5(textOutput("cell1")), width = 6, color = "red")
                                         )
                                         , column(6,
                                                  textOutput("cell2Cluster")
                                                  #valueBox(subtitle = h4(""),value = h5(textOutput("cell2")), width = 6, color = "aqua")
                                         )
                                       )
                                       , br()
                                       , actionButton("showNewIdentCluster", label = "Show new cluster")
                                       , hr()
                                       , textInput("nameClusteringIdentCluster", label = h3("Name of ident"), value = "(Ej Arg1Cluster)")
                                       , textInput("nameClusteringCell1Cluster", label = h6("Name of cell1"), value = "(Ej Arg1+)")
                                       , textInput("nameClusteringCell2Cluster", label = h6("Name of cell2"), value = "(Ej Arg1-)")
                                       , actionButton("saveIdentCluster", label = "Save ident")
                                       , checkboxInput("checkboxKeepIdentCluster", label = "Keep old ident", value = F)
                                       , hr()
                                       , actionButton("actionMksCellsCluster", label = "Run mks")

                                       #, downloadButton("downloadPlotCluster", "Download plot as pdf")
                                       )
                              )
                          )
               ),
               tabPanel("Gene",
                        fluidRow(
                          tabBox(title = ""
                                 , width = 8
                                 , tabPanel("Graph"
                                            , fluidRow(
                                              column(2, offset = 10
                                                     , downloadButton("downloadPlotGene", "")
                                              )
                                            )
                                            , plotlyOutput("plotGene") %>% withSpinner()
                                            # box(title="Distribution of reads",
                                            #       textOutput("brush")
                                            #)
                                            , fluidRow(
                                              column(2, offset = 10
                                                     , dropdownButton(
                                                       sliderInput("sliderScale", label = h3("Scale"), min = 0, 
                                                                   max = 100, value = c(01, 99))
                                                       , icon = icon(name = "cog"), up = F, right = T
                                                     )
                                              )
                                            )
                                            
                                 )
                                 , tabPanel("new Graph"
                                            , fluidRow(
                                              column(2, offset = 10
                                                     , downloadButton("downloadPlotCells", "")
                                              )
                                            )
                                            , plotlyOutput("plotNewClusterCell") %>% withSpinner()
                                            
                                            
                                 )        
                                 ,  tabPanel('Table'
                                             , fluidRow(
                                               column(2, offset = 10
                                                      , downloadButton("downloadDataMksCells", "")
                                               )
                                             )
                                             , br()
                                             , DT::dataTableOutput('df.mks.cell')
                                             
                                 )
                                 
                                 ),
                          box(
                            side = "right",
                            width = 4,
                            h4("Selected cells")
                            , textOutput("brush")
                            , hr()
                            , fluidRow(
                              column(4
                                     , actionButton("saveCell1", label = "Save Cell1")
                                     
                              ),
                              column(4
                                     , actionButton("saveCell2", label = "Save Cell2")
                              )
                            )
                            , br()
                            , fluidRow(
                              column(6
                                     , textOutput("cell1")
                                     #valueBox(subtitle = h4(""),value = h5(textOutput("cell1")), width = 6, color = "red")
                              )
                              , column(6,
                                       textOutput("cell2")
                                       #valueBox(subtitle = h4(""),value = h5(textOutput("cell2")), width = 6, color = "aqua")
                              )
                            )
                            , hr()
                            , actionButton("showNewIdentCell", label = "Show new cluster")
                            , hr()
                                     , textInput("nameClusteringIdent", label = h3("Name of ident"), value = "(Ej Arg1Cluster)")
                                     , textInput("nameClusteringCell1", label = h6("Name of cell1"), value = "(Ej Arg1+)")
                                     , textInput("nameClusteringCell2", label = h6("Name of cell2"), value = "(Ej Arg1-)")
                                     , actionButton("saveIdent", label = "Save ident")
                                     , checkboxInput("checkboxKeepIdent", label = "Keep old ident", value = F)
                                     , hr()                               
                                     , h3("New constrast") 
                                     , actionButton("actionMksCells", label = "Run mks")
                          )
                        )
               )
               , navbarMenu("More plots"
                            , tabPanel("Violin Plots"
                                                         ,box(width = 10, 
                                                         fluidRow(
                                                           column(2, offset = 10
                                                                  , downloadButton("downloadVlnPlot", "")
                                                           )
                                                         )
                                                         , plotOutput("plotVln") %>% withSpinner()
                                                         , fluidRow(
                                                           column(2, offset = 10
                                                                  , dropdownButton(
                                                                    checkboxInput("checkboxShowPoint", label = "Show point", value = F)
                                                                    
                                                                    , icon = icon(name = "cog"), up = T, right = F
                                                                  )
                                                           )
                                                         )
                                                         )
                             )
                            , "----"
                            , tabPanel("HeatMaps"
                                       , sidebarLayout(position = "right"
                                                       , sidebarPanel(
                                                         radioButtons("radioHMap", label = h3("Which genes use")
                                                                      , choices = list("My geneset" = 0, "Top mks" = 1)
                                                                      , selected = "tsne")
                                                         , sliderInput("sliderNumMksHmap", label = h3("Number of Mks"), min = 0, max = 20, value = 10)
                                                         , actionButton("actionMksHMap", label = "Show HeatMap")
                                                         , hr()
                                                       )
                                                       , mainPanel(
                                                        box(width = 10, 
                                                        fluidRow(
                                                           column(2, offset = 10
                                                                  , downloadButton("downloadHMap", "")
                                                           )
                                                         )
                                                         , plotOutput("HeatMapPlot") %>% withSpinner()
                                                       )
                                                       )
                                       )
                            )
               )
               # , navbarMenu("Tutorial"
               #              , tabPanel("Tutorial1"
               #                         , fluidPage(
               #                                      fluidRow(
               #                                          box(title = "prueba")
               #                                         )
               #                         )
               #              )
               #              , "----"
               #              , tabPanel("Tutorial2"
               #                         , fluidPage(
               #                           fluidRow(
               #                           )
               #                         )
               #              )
               # )
               , tabPanel("About SCADV",
                          fluidPage(
                            fluidRow(
                              box(h3("single-cell Data Analysis and VISualization")
                              , width = 8
                              , br()
                              ,p("The aim of this tool is to help scientific who are studying single cell and have a previous analysis.")
                              , br()
                              , p("
                                Firstly, scDAVIS allows to visualize metadata and qc (quality control) information about the experiment; and plot it in different graph such as histogram and
                                violin plots. Moreover, the previous analysis should be get tsne and/or pca information in order to be able to plot the different cluster that 
                                are performed. Unless, this tool permit to create your own cluster by selecting manually the cells.
                                ")
                              , br()
                              , p("A different category that scDAVIS offers is the differencial expression analysis. It can be do with predefined clustering or with a custom clustering.
                                Then it can be visualized in tsne or pca graph by representing the gene expression or in violin plot as well as HeatMap.")
                              , br()
                              , p("
                                 This tool is created by Juan Luis Onieva and Carlos Torroja with the collaboration of CNIC"))
                            )
                          )
                          
               )
    )
  )
  
)
)



