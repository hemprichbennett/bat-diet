library(here)
library(reshape2)
library(corrplot)

setwd(here())
taxa_mat<- read.csv('shiny/order_composition/order_mat.csv', row.names = 1) 



# Define UI for miles per gallon app ----
ui <- fluidPage(
  
  # App title ----
  headerPanel("Minimum number of assigned MOTU per bat (across species)"),
  
  # Sidebar panel for inputs ----
  sidebarPanel(
    # Input: Selector for variable to plot against mpg ----
    sliderInput("min_nodes", "Minimum number of MOTU consumed:",
                min = 0, max = 27,
                value = 1)
    
  ),
  
  # Main panel for displaying outputs ----
  mainPanel("",
            fluidRow(plotOutput("corr")
              #splitLayout(cellWidths = c("60%", "40%"), plotOutput("corr"), plotOutput("bars"))
            )
    
    
    #a <- renderPrint({cat('the maximum number of samples available for each network is', max_samples)}),
    #a()
    
    
    
  ) 
)



# Define server logic to plot various variables against mpg ----
server <- function(input, output) {
  
  output$caption <- renderText({
    paste("Minimum of", input$min_nodes, 'nodes', sep =' ')
  })
  
  dat <- reactive({
    min <- as.numeric(input$min_nodes)
    
    print(input)
    
    for_bigmat <- t(taxa_mat[,which(colSums(taxa_mat)>min)])
    if(0 %in% colSums(for_bigmat)){
      for_bigmat <- for_bigmat[,-which(colSums(for_bigmat)==0)]
    }
    
    big_cor <- for_bigmat
    bigcormat <- round(cor(big_cor),2)
    resbig <- cor.mtest(for_bigmat) 
    
    molten <- melt(for_bigmat)
    colnames(molten) <- c('sample', 'order', 'value')
    print(bigcormat)
    
  })
  
  
  
  mel <- reactive({
    min <- as.numeric(input$min_nodes)
    
    print(input)
    
    for_bigmat <- t(taxa_mat[,which(colSums(taxa_mat)>min)])
    if(0 %in% colSums(for_bigmat)){
      for_bigmat <- for_bigmat[,-which(colSums(for_bigmat)==0)]
    }
    
    
    molten <- melt(for_bigmat)
    colnames(molten) <- c('sample', 'order', 'value')
    
    molten[which(molten$order %in% c('Lepidoptera', 'Diptera', 'Coleoptera')),]
    
  })
  
  
  output$corr <- renderPlot(corrplot(dat(), method = "circle", p.mat = resbig$p, sig.level = .05, type = 'upper', order = 'AOE',
                                     tl.col = "black", tl.srt = 45, insig = 'blank',
                                     bg = "black"))
  
  output$dots <- renderPlot(ggplot(mel(), aes(x=sample, y = value, colour=order))+ geom_point(alpha= 0.7)+
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
      geom_line()+labs(y='N MOTU'))
  
  output$bars <- renderPlot(ggplot(mel()) + geom_bar(aes(y = value, x = sample, fill = order), data = molten[which(molten$order %in% c('Lepidoptera', 'Diptera', 'Coleoptera')),],
                      stat="identity")+
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    )
  
  
  
}

shinyApp(ui, server)

