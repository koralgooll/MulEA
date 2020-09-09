# Eszter's first example.
rm(list = ls())
root_of_cpp_project <- 'D:/projects/science/enrichment-analysis'
root_of_cpp_project <- '/home/cezary/Cezary/projects/science/enrichment-analysis'
setwd(root_of_cpp_project)


convert_model_list_to_model_df <- function(model_list) {
  model_df <- plyr::adply(.data = names(DB), .margins=1, .fun = function(list_name){
    
    data.frame('ontologyId' = list_name, 'ontologyName' = list_name, 
               'listOfValues' = I(DB[list_name]), stringsAsFactors = FALSE)
  })[c('ontologyId', 'ontologyName', 'listOfValues')]
  model_df
}


# Test Case 01
source("src/set-based-enrichment-test.R") # It contains the R source of enrichment analizis and it compiles the C++ part.
load(file="test-parameters-01.RData") # It loas a dummy databaseinto variabla 'DB', and a 'pool' and 'selet' lists; The DB contains only disjoint sets.  One of the sets is empty.
number_of_steps <- 100

your_res_01 <- set.based.enrichment.test(steps=number_of_steps, pool=pool, 
                                         select=select, DB=DB)

model_df <- convert_model_list_to_model_df(DB)
mulea_ora_01 <- MulEA::ORA(gmt = model_df, testData = select,
                           pool = pool, adjustMethod = "PT",
                           numberOfPermutations = number_of_steps)
mulea_res_01 <- MulEA::runTest(mulea_ora_01)


# Merge model dataframe with mulea results.
model_with_res <- merge(x = model_df, y = mulea_res_01, 
                        by.x = "ontologyId", by.y = "DB_names", all = TRUE)
names(model_df)
names(mulea_res_01)

model_with_res[[1,'listOfValues']]

enrichplot::cnetplot()

install.packages("ggnetwork")


# Plot Processing.
library(ggplot2)
library(dplyr)


bar_plot_mulea <- function(result_data_frame, selection_vector=1:10, 
                           categories_names='DB_names', probabilities_values='FDR') {
  result_data_frame[selection_vector,] %>% 
    ggplot( aes_string(x=categories_names, y=probabilities_values, fill=probabilities_values)) +
    geom_bar(stat="identity") +
    scale_fill_gradient2(mid='darkgreen', high='red') +
    coord_flip()
}

bar_plot_mulea(result_data_frame = mulea_res_01, selection_vector = c(1, 7, 15, 24, 27, 51, 61, 84, 86))




# Thwo variables encoded.
library(reshape2)

mulea_res_01.long <- melt(mulea_res_01[selection_vector,])
mulea_res_01.long$DB_names <- as.character(mulea_res_01.long$DB_names)
mulea_res_01.long$variable <- as.character(mulea_res_01.long$variable)

mulea_res_01.long %>% filter(variable %in% c('P', 'FDR')) %>% 
  ggplot(aes(DB_names, value, alpha=variable, fill=value)) +
  geom_bar(stat="identity", position = "dodge") +
  scale_alpha_ordinal(range = c(0.5, 0.9)) +
  scale_fill_gradient2(mid='red', high='darkgreen', space='Lab') +
  coord_flip() + 
  xlab("") +
  theme_bw()






# Graph's plots
install.packages(c("tidyverse", "tidygraph", "network", "igraph"))
library(tidyverse)
edge_list <- tibble(from = c(1, 2, 2, 3, 4), to = c(2, 3, 4, 2, 1))
node_list <- tibble(id = 1:4)

edge_list

library(network)
?network()
edges <- data.frame(from = c(1, 2, 2, 3, 4), to = c(2, 3, 4, 2, 1), weight = c(1, 7, 4, 7, 12))
nodes <- data.frame(id = c(1, 2, 3, 4), label = c("l-one", "l-two", "l-third", "l-fourth"), stringsAsFactors = FALSE)
# routes_network <- network(edges, vertex.attr = nodes, matrix.type = "edgelist", ignore.eval = FALSE)
routes_network <- network(edges, matrix.type = "edgelist", ignore.eval = FALSE)
plot(routes_network, vertex.cex = 3)
plot(routes_network, vertex.cex = 3, mode = "circle")

library(igraph)
routes_igraph <- graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)
plot(routes_igraph, edge.arrow.size = 0.2)
plot(routes_igraph, layout = layout_with_graphopt, edge.arrow.size = 0.2)


library(tidygraph)
library(ggraph)
routes_tidy <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
routes_igraph_tidy <- as_tbl_graph(routes_igraph)

ggraph(routes_tidy) + geom_edge_link() + geom_node_point() + theme_graph()
ggraph(routes_tidy, layout = "graphopt") + 
  geom_node_point() +
  geom_edge_link(aes(width = weight), alpha = 0.8) + 
  scale_edge_width(range = c(0.2, 2)) +
  geom_node_text(aes(label = label), repel = TRUE) +
  labs(edge_width = "Overlaped genes") +
  theme_graph()

ggraph(routes_igraph, layout = "linear") + 
  geom_edge_arc(aes(width = weight), alpha = 0.8) + 
  scale_edge_width(range = c(0.2, 2)) +
  geom_node_text(aes(label = label)) +
  labs(edge_width = "Letters") +
  theme_graph()

