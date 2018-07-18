library(here)
library(reshape2)
library(corrplot)

setwd(here())
#taxa_mat<- read.csv('shiny/order_composition/order_mat.csv', row.names = 1) 


all_ecology <- read.csv('data/output_data/all_bats/sitewise_all_individual_info.csv')
all_ecology <- all_ecology[,c(seq(24,40), 59)]


unique(all_ecology$Species)



# Define UI for miles per gallon app ----
ui <- fluidPage(
  
  # App title ----
  headerPanel("Minimum number of assigned MOTU per bat (across species)"),
  
  # Sidebar panel for inputs ----
  sidebarPanel(
    # Input: Selector for variable to plot against mpg ----
    sliderInput("min_nodes", "Minimum number of MOTU consumed:",
                min = 0, max = 27,
                value = 1),
  selectInput("species", "Species:", 
              choices=c('all', as.character(unique(all_ecology$Species))))
  ),
  
  # Main panel for displaying outputs ----
  mainPanel("",
            fluidRow(plotOutput("corr")
                     
            )
            
            
            
            
            
  ) 
)



# Define server logic to plot various variables against mpg ----
server <- function(input, output) {
  
  #output$caption <- renderText({
  #  paste("Minimum of", input$min_nodes, 'nodes', sep =' ')
  #})
  # 
   dat <- reactive({
     min <- as.numeric(input$min_nodes)
     sp <- as.character(input$species)
     
     print(min)
     print(sp)
     if(sp =='all' | sp== ''){
       taxa_mat <- all_ecology
       taxa_mat$Species <- NULL
     }else{
       taxa_mat <- all_ecology[which(all_ecology$Species== sp),]
       taxa_mat$Species <- NULL
     }
     taxa_mat
     
     for_bigmat <- t(taxa_mat[,which(colSums(taxa_mat)>min)])
     if(0 %in% colSums(for_bigmat)){
       for_bigmat <- for_bigmat[,-which(colSums(for_bigmat)==0)]
     }
     
     big_cor <- for_bigmat
     bigcormat <- round(cor(big_cor),2)
     #print('made bigcormat')
     resbig <- cor.mtest(for_bigmat) 
     #print('made resbig')
     #print(resbig)
     #molten <- melt(for_bigmat)
     #colnames(molten) <- c('sample', 'order', 'value')
     bigcormat
     
   })
   
  
  # 
  # dat <- reactive({
  #   
  # 
  #   min <- as.numeric(input$min_nodes)
  #                         sp <- as.character(input$species)
  #                         
  #                         print(input)
  #                         print(sp)
  #                         if(sp =='all' | sp== ''){
  #                           taxa_mat <- all_ecology
  #                           taxa_mat$Species <- NULL
  #                         }else{
  #                           taxa_mat <- all_ecology[which(all_ecology$Species== sp),]
  #                           taxa_mat$Species <- NULL
  #                         }
  #                         taxa_mat
  #                         
  #                         for_bigmat <- t(taxa_mat[,which(colSums(taxa_mat)>min)])
  #                         if(0 %in% colSums(for_bigmat)){
  #                           for_bigmat <- for_bigmat[,-which(colSums(for_bigmat)==0)]
  #                         }
  #                         
  #                         big_cor <- for_bigmat
  #                         bigcormat <- round(cor(big_cor),2)
  #                         #print('made bigcormat')
  #                         resbig <- cor.mtest(for_bigmat) 
  #                         #print('made resbig')
  #                         #print(resbig)
  #                         #molten <- melt(for_bigmat)
  #                         #colnames(molten) <- c('sample', 'order', 'value')
  #                         print(bigcormat)
  #                         
  # 
  # 
  #   # renderPlot({corrplot(bigcormat, method = "circle", p.mat = resbig$p, sig.level = .05, type = 'upper', order = 'AOE',
  #   #                                  tl.col = "black", tl.srt = 45, insig = 'blank',
  #   #                                  bg = "black")})
  # 
  # })  
  # 
  output$corr <- renderPlot({corrplot(dat(), method = "circle", p.mat = resbig$p, sig.level = .05, type = 'upper', order = 'AOE',
                       tl.col = "black", tl.srt = 45, insig = 'blank',
                       bg = "black")})
  
  
}


shinyApp(ui, server)

