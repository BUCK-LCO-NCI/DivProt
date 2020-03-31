#DivProt networks Shiny App

library(shiny)
library(shinythemes)
library(pheatmap)
library(igraph)
#THIS IS THE UI OBJECT
# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  theme = shinytheme("flatly"),
  # App title ----
  titlePanel("DivProt Networks"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
     
       # Input: Select a file 
      fileInput("file1", "Upload your matrix:",
                multiple = TRUE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv", ".tsv", ".txt")),
      
      # Horizontal line ----
      tags$hr(),
      
      # Input: Slider for the number of bins ----
      sliderInput(inputId = "edge_weight",
                  label = "Edge weight scaling variable:",
                  min = 0.0, #0% scaling
                  max = 7.0, #1000% scaling
                  value = 1.0,
                  step = 0.1),
      
      sliderInput(inputId = "avg_aa_len",
                  label = "Average amino acid length of input seqs:",
                  min = 20.0, #20 aa blah
                  max = 1000.0, #1000 aa len
                  value = 150.0,
                  step = 25.0),
      
      # Horizontal line ----
      tags$hr(),
      
      # Input: Checkbox for delete nodes if edge = 0 ----
      checkboxInput("checkbox_results", "Delete nodes with no edges", value=FALSE)
      
    ),
    
    
  
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Network
      plotOutput(outputId = "plot")
      
    )
  )
)


 
#THIS IS THE SERVER OBJECT
server <- function(input, output) {

  output$plot <- renderPlot({

  Ts <- as.numeric(input$edge_weight)
  
  Tem_custom <- ((3.657* input$avg_aa_len) + 260.1) + (((3.657* input$avg_aa_len) + 260.1) * Ts) 
  
  #matrix_confirm_temp = read.csv("../expanded_koonin_out/final_time/Final_3_align_matrix.csv", header = TRUE, row.names = 1, check.names = FALSE)
  req(input$file1) #nned this not to have error for file when app is first launched
  matrix_confirm_temp = read.csv(input$file1$datapath,
                                 header = TRUE, 
                                 row.names = 1, 
                                 check.names = FALSE)
 
   #1. Network construction
  #formatting
  matrix_confirm <- as.matrix(matrix_confirm_temp)
  mode(matrix_confirm) <- "numeric"   
  
  #actual igraph, two for different layout options
  library("igraph")
  ig <- graph.adjacency(matrix_confirm, mode="undirected", weighted=TRUE, diag = FALSE, add.colnames = NULL)

  
    if (input$checkbox_results) {
      
      te <- delete_edges(ig, E(ig)[weight< Tem_custom])
      dv <- delete.vertices(te, V(te)[degree(te)==0])
      community_clustering <- multilevel.community(dv)
      cluster_colors <- rainbow(max(membership(community_clustering)), alpha = 0.5)
      lld <- layout_with_graphopt(dv, niter=1000,charge = 0.01) #this is BY FAR the best layout here
      
      #plot with edge weight applied and zero edge weight nodes deleted
      
      plot(dv,
          layout=lld, 
          edge.arrow.size=0.5, 
          vertex.label.cex=1.0, 
          vertex.label.family="Helvetica",
          vertex.label.font=1.5,
          vertex.label.dist=1,
          vertex.shape="circle", 
          vertex.size=3,
          vertex.color=cluster_colors[membership(community_clustering)], 
          vertex.label.color="black", 
          edge.width=0.4
          )

      } else {
  #output$net <- renderPlot({
        
        community_clustering <- fastgreedy.community(ig)
        cluster_colours <- rainbow(max(membership(community_clustering)), alpha = 0.5)
        
        ll <- layout.fruchterman.reingold(ig, niter=1000)
        #2. Apply edgeweight cutoff,
        te <- delete_edges(ig, E(ig)[weight< Tem_custom])
        community_clustering <- fastgreedy.community(te)
        cluster_colours <- rainbow(max(membership(community_clustering)), alpha = 0.5)
        te2 <- layout.fruchterman.reingold(te, niter=1000)
        
            plot(te,
            layout=te2,
            edge.arrow.size=0.5,
            vertex.label.cex=1.0,
            vertex.label.family="Helvetica",
            vertex.label.font=1.5,
            vertex.label.dist=1,
            vertex.shape="circle",
            vertex.size=3,
            vertex.color=cluster_colours[membership(community_clustering)],
            vertex.label.color="black",
            edge.width=0.3
            ) }
    
    
    }, height = 900, width = 900) 
    
 
}


#THIS PUTS IT TOGETHER
shinyApp(ui = ui, server = server)


