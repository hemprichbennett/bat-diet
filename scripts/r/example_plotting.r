#### Header ####
## Project: 
## Script purpose: 
## Date: 
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################
options(warnPartialMatchArgs = TRUE)

#
library(tidyverse)
library(ggplotify)
library(igraph)
library(gridExtra)
library(grid)

# Setup -------------------------------------------------------------------


# Load in the external data
nest <- read.csv("data/dummy_data/nested_2.csv")
mod <- read.csv("data/dummy_data/modularity.csv")



# Functions ---------------------------------------------------------------

# A function which turns the input matrix into a nicely formatted edgelist,
# which ggplot will happily play with
matrix_formatting <- function(mat) {
  out <- mat %>%
    # Rename the prey species column
    rename(`Prey species` = X) %>%
    # make data long format, instead of wide
    gather(-`Prey species`, key = "Predator species", value = "Interacts") %>%
    # turn the interaction column into T/F instead of binary
    # mutate(Interacts = ifelse(Interacts == 1, T, F)) %>%
    # Reorder the factors, as ggplot is plotting them wrong
    mutate(`Prey species` = fct_rev(`Prey species`))


  return(out)
}


# A function to make a random edgelist with the exact same number of nodes and
# connectance as the input edgelist
non_fun <- function(x) {
  
  # Decide which rows should have interactions
  to_fill <- sample(
    seq(1, length(x$Interacts)),
    sum(x$Interacts)
  )
  
  output <- x %>%
    # Make all the interactions empty before adding the new ones
    mutate(Interacts = 0)
  
  # Add the new interactions
  output$`Interacts`[to_fill] <- 1
  
  return(output)
}


# Bipartite data manipulation ---------------------------------------------

# set up the object to be used for plotting modularity and nestedness

Nested <- matrix_formatting(nest)
Modular <- matrix_formatting(mod)

set.seed(3)

`Non-nested` <- non_fun(Nested)
`Non-modular` <- non_fun(Modular)

outlist <- list()
nets <- list(Modular, `Non-modular`, Nested, `Non-nested`)
names(nets) <- c('Modular', 'Non-modular', 'Nested', 'Non-nested')

for (i in 1:length(nets)) {
  nets[[i]]$dataset <- names(nets)[i]
  
}

all_df <- do.call(rbind, nets) %>%
  # Turn the interaction column from a binary variable to a T/F
  mutate(Interacts = ifelse(Interacts == 1, T, F)) %>%
  # Now turn it into a factor, where T is before F
  mutate(Interacts = fct_rev(as.factor(Interacts))) %>%
  # Rename the predators and prey for the plotting
  mutate(`Predator species` = gsub('Pred', 'Predator ', `Predator species`),
         `Prey species` = gsub('Bug', 'Prey ', `Prey species`)) %>%
  # Reverse the factor orders for prey species, as ggplot plots them from bottom
  # to top, not vice versa as would be helpful here
  mutate(`Prey species` = fct_rev(`Prey species`)) %>%
  mutate(dataset = gsub('Modular', 'A) Modular', dataset),
        dataset = gsub('Non-modular', 'B) Non-modular', dataset),
        dataset = gsub('Nested', 'C) Nested', dataset),
        dataset = gsub('Non-nested', 'D) Non-nested', dataset))


# Plot the bipartite networks ---------------------------------------------

bipartite_nets <- ggplot(all_df, aes(x = `Predator species`, y = `Prey species`, fill = Interacts)) +
  # Add the tiles. 'Colour = black' adds gridlines
  geom_tile(colour = "black") +
  scale_fill_manual(values = c("black", "#f2f2f2"),
                    labels = c('True', 'False')) +
  facet_wrap(. ~ dataset, strip.position = 'bottom') +
  theme_classic() +
  theme(
    # Put the legend at the bottom
    legend.position = "bottom",
    # # Get rid of the axis labels
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = -0.1),
    # Get rid of the background colours on facet labels
    strip.background = element_rect(
      color="black", fill="white", size=1, linetype="solid"),
    panel.spacing = unit(2, "lines"),
    text = element_text(size=20)
  ) +
  scale_x_discrete(position = "top") 

bipartite_nets
ggsave('example_images/biparite_networks.png', bipartite_nets,
       width = 9.5, height = 9)

# Unipartite network ------------------------------------------------------

# # Make network
# close_graph <- graph(edges = c(
#   1, 3,
#   2, 3,
#   4, 3,
#   5, 3
# ), n = 5, directed = F)
# 
# # Print network stats
# igraph::closeness(close_graph)
# igraph::betweenness(close_graph)
# 
# # Colour the nodes, with node 3 being black
# V(close_graph)$color <- ifelse(igraph::betweenness(close_graph)==0, 'yellow', 'blue')
# # Remove the vertex labels
# V(close_graph)$label <- NA
# # Make the edges wider for nicer plotting
# E(close_graph)$width <- 10
# 
# # Save image
# pdf('example_images/centrality.pdf', width = 15, height = 15)
# plot(close_graph, pch = 19, vertex.size = 25)
# dev.off()

far_graph <- graph(edges = c(
  1, 2,
  2, 3,
  4, 3,
  5, 3
), n = 5, directed = F)

# Make a vector of colours, then use it to colour the nodes
colrs <- c("gray50", "tomato", "gold", 'blue')
colscheme <- colrs[as.factor(closeness(far_graph))]
V(far_graph)$color <- colscheme
# Remove the labels
V(far_graph)$label <- ""
E(far_graph)$label <- ""
# Make the edges wider
E(far_graph)$width <- 10

plot(far_graph)


# Store the igraph plot as a grob, so we can use it in gridextra
graph_grob <- ggplotify::base2grob(~plot(far_graph))


# Make a graphical object of the table of centrality values --------------------



betweenness_vals <- igraph::betweenness(far_graph)
closeness_vals <- round(igraph::closeness(far_graph), 2)
colscheme

# bipartite::specieslevel(matrix(c(1,1,1,1), nrow = 1))
# bipartite::specieslevel(matrix(c(1,1,0,
#          1,0,1), ncol = 3), index = c('betweenness', 'closeness'))

unipartite_df <- data.frame('Node' = colscheme,
                            'Betweenness centrality' = betweenness_vals,
                            'Closeness centrality' = closeness_vals)

# store an alternate version of unipartite_df, with no values in the first 
# column. This is because when plotting the table, we want those cells to be 
# empty, so we can just fill it with the appropriate colours. 
# Also, rename the columns

temp_df <- unipartite_df %>%
  mutate(Node = '') %>%
  rename(`Betweenness centrality` = Betweenness.centrality,
         `Closeness centrality` = Closeness.centrality)


table_for_plot <- tableGrob(
  temp_df,
  rows = NULL,
  theme = ttheme_minimal()
)

# # Alter the formatting further, using adding in column separators
# separators <- replicate(ncol(table_for_plot) - 1,
#                         segmentsGrob(x1 = unit(0, "npc"), gp=gpar(lty=2)),
#                         simplify=FALSE)
# 
# table_for_plot <- gtable::gtable_add_grob(table_for_plot, grobs = separators,
#                              t = 2, b = nrow(table_for_plot), l = seq_len(ncol(table_for_plot)-2)+2)

# The following code is awful, but it works. It sets the colours of each
# node cell, so that they match the colour scheme in the rest of the figure.
# Code modified from the 'Accessing existing grobs in the table' section of 
# https://cran.r-project.org/web/packages/gridExtra/vignettes/tableGrob.html
find_cell <- function(table, row, col, name="core-fg"){
  l <- table$layout
  which(l$t==row & l$l==col & l$name==name)
}

ind1 <- find_cell(table_for_plot, 2, 1, "core-bg")
table_for_plot$grobs[ind1][[1]][["gp"]] <- gpar(fill="gray50", lwd=5)

ind2 <- find_cell(table_for_plot, 3, 1, "core-bg")
table_for_plot$grobs[ind2][[1]][["gp"]] <- gpar(fill="gold",  lwd=5)

ind3 <- find_cell(table_for_plot, 4, 1, "core-bg")
table_for_plot$grobs[ind3][[1]][["gp"]] <- gpar(fill="blue",  lwd=5)

ind4 <- find_cell(table_for_plot, 5, 1, "core-bg")
table_for_plot$grobs[ind4][[1]][["gp"]] <- gpar(fill="tomato",  lwd=5)

ind5 <- find_cell(table_for_plot, 6, 1, "core-bg")
table_for_plot$grobs[ind5][[1]][["gp"]] <- gpar(fill="tomato",  lwd=5)

#grid.draw(table_for_plot)



#table_for_plot$grobs[10][[1]][["gp"]] <- gpar(fill = 'blue')

# Store the layout for the plots to be arranged in
layout <- rbind(c(1,2),
                c(1,3))

# Put all the objects together into a grid object
my_grid <- grid.arrange(
  bipartite_nets, 
  graph_grob, 
  table_for_plot,
             layout_matrix = layout)
# Save it
ggsave('example_images/example_grid.pdf', my_grid, width = 16)
