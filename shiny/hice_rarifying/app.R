library(shiny)
library(ggplot2)
library(here)
library(ggridges)
library(reshape2)
#library(EnvStats)
#setwd(here())
hice_data<- read.csv('hice_rarified.csv') 
hice_data$bat.ID <- as.factor(hice_data$bat.ID)
hice_data <- hice_data[,-1]

hice_data$habitat_type <- rep(NA, nrow(hice_data))
hice_data[grep('SAFE', hice_data$Network), 'habitat_type'] <- 'Logged'
hice_data[grep('SBE', hice_data$Network), 'habitat_type'] <- 'Logged, replanted'
hice_data[grep('Danum', hice_data$Network), 'habitat_type'] <- 'Primary'
hice_data[grep('Maliau', hice_data$Network), 'habitat_type'] <- 'Primary'

max_samples <- c()
for(i in 1:length(unique(hice_data$Network))){
  m <- max(hice_data$N_bats[hice_data$Network== unique(hice_data$Network)[i]])
  names(m) <- unique(hice_data$Network)[i]
  #print(m)
  max_samples <- c(max_samples, m)
}


# Define UI for miles per gallon app ----
ui <- fluidPage(
  
  # App title ----
  headerPanel("Hipposideros cervinus rarifying"),
  
  # Sidebar panel for inputs ----
  sidebarPanel(
    # Input: Selector for variable to plot against mpg ----
    sliderInput("max_nodes", "maximum number of samples:",
                min = 20, max = 118,
                value = 1),
    
     sliderInput("iteration", "iteration:",
                 min = 1, max = 100,
                 value = 1)
    
    
    
  ),
  
  # Main panel for displaying outputs ----
  mainPanel(
    # Output: Formatted text for caption ----
    h3(textOutput("caption")),
    
    # Output: Plot of the requested variable against mpg ----
    tabsetPanel(
      tabPanel("Plot", plotOutput("facet_ridge"))#,
      #tabPanel("Table", tableOutput('table'))
    )
    
    
    #a <- renderPrint({cat('the maximum number of samples available for each network is', max_samples)}),
    #a()
    
    
    
  )
)



# Define server logic to plot various variables against mpg ----
server <- function(input, output) {
  
  output$caption <- renderText({
    paste("Maximum of", input$max_nodes, 'nodes', sep =' ')
  })
  
  dat <- reactive({
    max <- as.numeric(input$max_nodes)
    it <- as.numeric(input$iteration)
    #max <- 15
    print(input)
    
    
    # Return the formula text for printing as a caption ----

    
    for_plot <- hice_data[hice_data$iteration==it,]
    
    to_filter <- names(which(max_samples >= max))
    use_max <- max_samples[which(max_samples < max)]
    
    df <- matrix(nrow= 0, ncol = ncol(for_plot))
    for(i in 1:length(use_max)){
      biggest <- for_plot[for_plot$Network == names(use_max)[i] & for_plot$N_bats == use_max[i],]
      df <- rbind(df, biggest)
    }
    
    
    
    for_plot <- rbind(for_plot[which(for_plot$N_bats == max & for_plot$Network %in% to_filter),], df)
    
    
    
    melted <- melt(for_plot, id.vars = c('Network', 'N_bats', 'habitat_type', 'bat.ID', 'iteration'))
    
    melted
    
  })
  
  output$facet_ridge <- renderPlot({ggplot(dat(), aes (y=Network, x =value, fill=habitat_type)) + 
    geom_density_ridges(scale= 0.5)+ #The scale determines the space between the rows
    theme_ridges()+ #This changes the theme to make it more aesthetically pleasing
    scale_fill_cyclical(values = c("#85d7da","#cfb4de","#d0ca9f"), guide = 'legend', name = 'Habitat type')+
    scale_x_continuous(expand = c(0.01, 0)) + #Make the space between the labels and plot smaller
    scale_y_discrete(expand = c(0.01, 0))+ #Make it so the top series actually fits in the plot
    ylab(NULL)+ xlab(NULL)+
    facet_wrap( ~ variable, ncol=1, scales = 'free_x', strip.position = 'bottom')+ #free_x is required so that the x-axes aren't all constrained to showing the same thing
    theme(strip.background = element_rect(fill="white"), strip.placement = "outside", panel.spacing = unit(0.8, "lines"),#strip stuff sorts the facet labels, spacing adjusts the space between facets
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          text = element_text(size=12))+
    theme(legend.text =element_text(size = 10))},height = 600,width = 400) #Trying to standardise the sizes across both plots
          #stat_n_text() +
  #facet_ridge
  #This is giving us some non-finite warnings, so something isn't right here
  
  output$table <- renderTable(melted[match(unique(melted$Network), melted$Network),c('Network', 'N_bats')])
  
}

shinyApp(ui, server)

