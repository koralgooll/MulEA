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

# Mock data
model_df <- data.frame('ontologyId'=c("CAT_g0001", "CAT_g0002", "CAT_g0003", "CAT_g0004", "CAT_g0005"), 
                       'ontologyName'=c("CAT_g0001", "CAT_g0002", "CAT_g0003", "CAT_g0004", "CAT_g0005"), 
                       'listOfValues'=I(list(
                         c("g_001", "g_002", "g_003", "g_004"), 
                         c("g_003", "g_004"), 
                         c("g_003", "g_004", "g_005", "g_006"), 
                         c("g_004", "g_005", "g_006", "g_007"),
                         c("g_001", "g_002", "g_007"))))
mulea_ora_01 <- MulEA::ORA(gmt = model_df, testData = c("g_002", "g_003", "g_004", "g_005", "g_008", "g_009"),
                           pool = c("g_001", "g_002", "g_003", "g_004", "g_005", "g_006", "g_007"), adjustMethod = "PT",
                           numberOfPermutations = 100)
mulea_res_01 <- MulEA::runTest(mulea_ora_01)

# Merge model dataframe with mulea results.
model_with_res <- merge(x = model_df, y = mulea_res_01, 
                        by.x = "ontologyId", by.y = "DB_names", all = TRUE)

# Create relaxed dataframe from our structure.
model_with_res_dt <- data.table::setDT(model_with_res)
model_with_res_dt_size = 0
for (i in 1:nrow(model_with_res_dt[,1])) {
  model_with_res_dt_size = model_with_res_dt_size + length(model_with_res_dt[[i, 'listOfValues']])
}

model_with_res_dt_relaxed <- data.table::data.table(
  ontologyId=rep('a',length.out=model_with_res_dt_size), 
  gen_in_ontology=rep('a',length.out=model_with_res_dt_size),
  ontology_p_stat=rep(1.0,length.out=model_with_res_dt_size))

model_with_res_dt_relaxed_counter = 1
for (i in 1:nrow(model_with_res_dt[,1])) {
  category_name <- model_with_res_dt[[i, 'ontologyId']]
  category_p_stat <- model_with_res_dt[[i, 'P']]
  for (item_name in model_with_res_dt[[i, 'listOfValues']]) {
    model_with_res_dt_relaxed[model_with_res_dt_relaxed_counter, 
                              c("ontologyId", "gen_in_ontology", "ontology_p_stat"):=list(category_name, item_name, category_p_stat)]
    model_with_res_dt_relaxed_counter = model_with_res_dt_relaxed_counter + 1
  }
}
model_with_res_dt_relaxed


relax_model_with_results <- function(model_with_res) {
  
}




ontologies <-unique(model_with_res_dt_relaxed[,'ontologyId'])
ontologies_graph_edges_num <- sum(1:(nrow(ontologies)-1))
ontologies_graph_edges <- data.table::data.table(
  from=rep('a', length.out=ontologies_graph_edges_num), 
  to=rep('a', length.out=ontologies_graph_edges_num),
  weight=rep(0, length.out=ontologies_graph_edges_num))


ontologies_graph_edges_counter <- 1
for (i in 1:(nrow(ontologies)-1)) {
  ontology_name_i <- ontologies[i, ontologyId]
  print(ontology_name_i)
  genes_in_ontology_i <- model_with_res_dt_relaxed[ontologyId==ontology_name_i, gen_in_ontology]
  print(genes_in_ontology_i)
  
  for (j in (i+1):nrow(ontologies)) {
    ontology_name_j <- ontologies[j, ontologyId]
    genes_in_ontology_j <- model_with_res_dt_relaxed[ontologyId==ontology_name_j, gen_in_ontology]
    genes_in_ontology_i_j_intersection_num <- length(intersect(genes_in_ontology_i, genes_in_ontology_j))
    print(genes_in_ontology_i_j_intersection_num)
    ontologies_graph_edges[ontologies_graph_edges_counter,
                           c('from', 'to', 'weight'):=list(ontology_name_i, ontology_name_j,
                                                           genes_in_ontology_i_j_intersection_num)]
    ontologies_graph_edges_counter = ontologies_graph_edges_counter + 1
  }
}
ontologies_graph_edges


nodes_ids <- c("CAT_g0001", "CAT_g0002", "CAT_g0003", "CAT_g0004", "CAT_g0005")
nodes_p_stat <- c(0.6, 0.4, 0.6, 0.8, 0.1)
ontologies_graph_nodes <- data.table::data.table(
  id=nodes_ids, 
  label=nodes_ids,
  p_stat=nodes_p_stat)


nodes <- data.frame(id = c("CAT_g0001", "CAT_g0002", "CAT_g0003", "CAT_g0004"), 
                    label = c("CAT_g0001", "CAT_g0002", "CAT_g0003", "CAT_g0004"), 
                    p_stat = c(0.6, 0.4, 0.6, 0.8, 0.1),
                    stringsAsFactors = FALSE)


routes_tidy <- tidygraph::tbl_graph(nodes = ontologies_graph_nodes, 
                                    edges = ontologies_graph_edges, directed = TRUE)
library(ggraph)

# ggraph(routes_tidy, layout = "graphopt") + 
#   geom_node_point() +
#   geom_edge_link(aes(width = weight), alpha = 0.8) + 
#   scale_edge_width(range = c(0.2, 2)) +
#   geom_node_text(aes(label = label), repel = TRUE) +
#   labs(edge_width = "Overlaped genes") +
#   theme_graph()
# 
# routes_tidy[,p_stat]
library(ggraph)
library(ggforce)
ggraph(routes_tidy, layout = "linear", circular = TRUE) +
  # geom_edge_arc(aes(width = weight, colour = weight), alpha = 0.6) +
  # scale_edge_color_gradient2(mid='black', high='brown') +
  geom_edge_arc(aes(width = weight), alpha = 0.5) +
  scale_edge_width(range = c(0, 3)) +
  geom_node_point(aes(color=p_stat)) +
  geom_node_point(aes(color=p_stat, size=(1-p_stat)), show.legend = FALSE) +
  scale_size_area(max_size = 10) +
  # guides(fill = guide_legend(reverse = TRUE)) +
  # guide
  # guides(size = guide_legend(override.aes = list(size = seq(0,10, by=2)) ) ) +
  scale_color_gradient2(mid='darkgreen', high='red') +
  geom_node_text(aes(label = label), repel = TRUE, fonts = "mono") +
  
  # theme_void()
  # theme_no_axes()
  # theme_minimal()
  # theme_light()
  # theme_classic()
  # theme_gray()
  # theme_bw()
  # theme_dark()
  # TODO : Ask Eszter about theme.
  theme_graph(base_family = "mono")

# scale_fill_gradient2(mid='darkgreen', high='red')

ggraph(routes_tidy, layout = 'kk', maxiter = 100) + 
  geom_edge_link(aes(colour = factor(weight))) + 
  geom_node_point()

ggraph(routes_tidy, layout = 'linear') + 
  geom_edge_arc(aes(colour = factor(weight)))

ggraph(routes_tidy, layout = 'linear', circular = TRUE) + 
  geom_edge_arc(aes(colour = factor(weight)))

# Add metadata to node
ggraph(routes_tidy, 'partition') + 
  geom_node_tile(aes(fill = p_stat), size = 0.25)


library(tidygraph)
library(ggraph)

edges <- data.frame(from = c(1, 2, 2, 3, 4), to = c(2, 3, 4, 2, 1), weight = c(1, 7, 4, 7, 12))
nodes <- data.frame(id = c(1, 2, 3, 4), label = c("l-one", "l-two", "l-third", "l-fourth"), stringsAsFactors = FALSE)

routes_tidy <- tidygraph::tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
# routes_igraph_tidy <- as_tbl_graph(routes_igraph)

ggraph::ggraph(routes_tidy) + ggraph::geom_edge_link() + ggraph::geom_node_point() + ggraph::theme_graph()

ggraph(routes_tidy) + 
  geom_node_point() +
  geom_edge_link(aes(width = weight), alpha = 0.8) + 
  scale_edge_width(range = c(0.2, 2)) +
  geom_node_text(aes(label = label), repel = TRUE) +
  labs(edge_width = "Overlaped genes") +
  theme_graph()






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
bar_plot_mulea(result_data_frame = mulea_res_01)



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

ggraph::ggraph(routes_tidy) + ggraph::geom_edge_link() + ggraph::geom_node_point() + ggraph::theme_graph()
ggraph(routes_tidy) + 
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



# Heatmap implementation.
library(ggplot2)

# Dummy data
x <- LETTERS[1:20]
y <- paste0("var", seq(1,20))
data <- expand.grid(X=x, Y=y)
data$Z <- runif(400, 0, 5)

# Heatmap---- 
ggplot(data, aes(X, Y, fill= Z)) + 
  geom_tile()


# Order categories by p-value, top is big p-values 
ggplot(model_with_res_dt_relaxed, aes(gen_in_ontology, ontologyId, fill= ontology_p_stat)) + 
  scale_fill_gradient2(mid='darkgreen', high='red', midpoint = 0.5) +
  geom_tile() +
  # theme_void()
  # theme_no_axes()
  # theme_minimal()
  theme_light()
  # theme_classic()
  # theme_gray()
  # theme_bw()
  # theme_dark()
  # TODO : Ask Eszter about theme.
  # theme_graph(base_family = "mono")


