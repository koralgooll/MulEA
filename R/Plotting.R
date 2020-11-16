
# PUBLIC API (Plotting)
#' @description
#' \code{relaxModelAndResults}
#'
#' \code{relaxModelAndResults} merge model and model relsuts into 
#' one relaxed datatable for easy resutls graphical interpretation.
#'
#' @param mulea_model MulEA object represents model. For example created by MulEA::ORA.
#' @param mulea_model_resuts Results from model, in most cases it is returned by MulEA::runTest generic method.
#' 
#' 
#' @title Input/Output Functions
#' @name  InputOutputFunctions
#' @export
#'
#' @return Return relaxed datatable where model and results are merged for plotting purposes. 
relaxModelAndResults <- function(mulea_model=NULL, mulea_model_resuts=NULL) {
  
  model_with_res <- merge(x = mulea_model@gmt, y = mulea_model_resuts, 
                          by.x = "ontologyId", by.y = "DB_names", all = TRUE)
  # Create relaxed dataframe from our structure.
  print(model_with_res)
  model_with_res_dt <- data.table::setDT(model_with_res)
  model_with_res_dt_size = 0
  print(model_with_res_dt[,1])
  for (i in 1:nrow(model_with_res_dt[,1])) {
    model_with_res_dt_size <- model_with_res_dt_size + length(model_with_res_dt[[i, 'listOfValues']])
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
}


# PUBLIC API (Plotting)
#' @description
#' \code{plotGraph}
#'
#' \code{plotGraph} merge model and model relsuts into 
#' one relaxed datatable for easy resutls graphical interpretation.
#'
#' @param mulea_relaxed_resuts data.table in relaxed form.
#' 
#' 
#' @title Input/Output Functions
#' @name  InputOutputFunctions
#' @export
#'
#' @return Return plot. 
plotGraph <- function(mulea_relaxed_resuts) {
  model_with_res_dt_relaxed <- mulea_relaxed_resuts
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
  
  
  nodes_ids <- model_with_res_dt_relaxed[,ontologyId]
  # nodes_ids <- c("CAT_g0001", "CAT_g0002", "CAT_g0003", "CAT_g0004", "CAT_g0005")
  nodes_p_stat <- model_with_res_dt_relaxed[,ontology_p_stat]
  # nodes_p_stat <- c(0.6, 0.4, 0.6, 0.8, 0.1)
  ontologies_graph_nodes <- data.table::data.table(
    id=nodes_ids, 
    label=nodes_ids,
    p_stat=nodes_p_stat)
  
  ontologies_graph_nodes <- unique(ontologies_graph_nodes)
  
  
  routes_tidy <- tidygraph::tbl_graph(nodes = ontologies_graph_nodes, 
                                      edges = ontologies_graph_edges, directed = TRUE)
  # library(ggraph)
  
  ggraph(routes_tidy, layout = "linear", circular = TRUE) +
    geom_edge_arc(aes(width = weight), alpha = 0.5) +
    scale_edge_width(range = c(0, 3)) +
    geom_node_point(aes(color=p_stat)) +
    geom_node_point(aes(color=p_stat, size=(1-p_stat)), show.legend = FALSE) +
    scale_size_area(max_size = 10) +
    scale_color_gradient2(mid='darkgreen', high='red') +
    geom_node_text(aes(label = label), repel = TRUE, fonts = "mono") +
    theme_graph(base_family = "mono")
}



