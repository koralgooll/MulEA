
filterRelaxedResultsForPlotting <- function(mulea_relaxed_resuts, 
                                            statistics_value_colname='ontology_p_stat', 
                                            statistics_value_cutoff=0.05) {
  
  mulea_relaxed_resuts_filtered_na <- mulea_relaxed_resuts[!is.na(mulea_relaxed_resuts[[statistics_value_colname]]),]
  mulea_relaxed_resuts_filtered_cutoff <- mulea_relaxed_resuts_filtered_na[mulea_relaxed_resuts_filtered_na[[statistics_value_colname]]<=statistics_value_cutoff,]
  mulea_relaxed_resuts_filtered_cutoff
}

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
relaxModelAndResults <- function(mulea_model=NULL, mulea_model_resuts=NULL, 
                                 mulea_model_ontology_col_name="ontologyId", 
                                 mulea_model_resuts_ontology_col_name='DB_names') {
  
  model_with_res <- merge(x = mulea_model@gmt, y = mulea_model_resuts, 
                          by.x = mulea_model_ontology_col_name, 
                          by.y = mulea_model_resuts_ontology_col_name, all = TRUE)
  # Create relaxed dataframe from our structure.
  model_with_res_dt <- data.table::setDT(model_with_res)
  model_with_res_dt_size = 0
  for (i in 1:nrow(model_with_res_dt)) {
    model_with_res_dt_size <- model_with_res_dt_size + length(model_with_res_dt[[i, 'listOfValues']])
  }
  model_with_res_dt_relaxed <- data.table::data.table(
    ontologyId=rep('a',length.out=model_with_res_dt_size), 
    gen_in_ontology=rep('a',length.out=model_with_res_dt_size),
    ontology_p_stat=rep(1.0,length.out=model_with_res_dt_size))
  
  model_with_res_dt_relaxed_counter = 1
  for (i in 1:nrow(model_with_res_dt)) {
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
plotGraph <- function(mulea_relaxed_resuts, edge_weight_cutoff=0, 
                      statistics_value_colname='ontology_p_stat', 
                      statistics_value_cutoff=0.05) {
  
  model_with_res_dt_relaxed <- MulEA:::filterRelaxedResultsForPlotting(mulea_relaxed_resuts=mulea_relaxed_resuts,
                                          statistics_value_colname=statistics_value_colname,
                                          statistics_value_cutoff=statistics_value_cutoff)
  
  ontologies <-unique(model_with_res_dt_relaxed[,'ontologyId'])
  ontologies_graph_edges_num <- sum(1:(nrow(ontologies)-1))
  ontologies_graph_edges <- data.table::data.table(
    from=rep('a', length.out=ontologies_graph_edges_num), 
    to=rep('a', length.out=ontologies_graph_edges_num),
    weight=rep(0, length.out=ontologies_graph_edges_num))
  
  if (0 == ontologies_graph_edges_num) {
    stop('No edges at all. Wrong data.table or manipulate statistics_value_cutoff please.')
  }
  
  ontologies_graph_edges_counter <- 1
  for (i in 1:(nrow(ontologies)-1)) {
    ontology_name_i <- ontologies[i, ontologyId]
    genes_in_ontology_i <- model_with_res_dt_relaxed[ontologyId==ontology_name_i, gen_in_ontology]
    
    for (j in (i+1):nrow(ontologies)) {
      ontology_name_j <- ontologies[j, ontologyId]
      genes_in_ontology_j <- model_with_res_dt_relaxed[ontologyId==ontology_name_j, gen_in_ontology]
      genes_in_ontology_i_j_intersection_num <- length(intersect(genes_in_ontology_i, genes_in_ontology_j))
      if (edge_weight_cutoff < genes_in_ontology_i_j_intersection_num) {
        ontologies_graph_edges[ontologies_graph_edges_counter,
                               c('from', 'to', 'weight'):=list(ontology_name_i, ontology_name_j,
                                                               genes_in_ontology_i_j_intersection_num)]
        ontologies_graph_edges_counter = ontologies_graph_edges_counter + 1
      }
    }
  }
  ontologies_graph_edges <- ontologies_graph_edges[0:(ontologies_graph_edges_counter-1), ]
  
  nodes_ids <- model_with_res_dt_relaxed[,ontologyId]
  nodes_p_stat <- model_with_res_dt_relaxed[,ontology_p_stat]
  ontologies_graph_nodes <- data.table::data.table(
    id=nodes_ids, 
    label=nodes_ids,
    p_stat=nodes_p_stat)
  
  ontologies_graph_nodes <- unique(ontologies_graph_nodes)
  
  routes_tidy <- tidygraph::tbl_graph(nodes = ontologies_graph_nodes,
                                      edges = ontologies_graph_edges, 
                                      directed = TRUE)
  
  graph_plot <- ggraph(routes_tidy, layout = "linear", circular = TRUE)
  
  if (0 != nrow(routes_tidy %>% tidygraph::activate(edges) %>% tidygraph::as_tibble())) {
    graph_plot <- graph_plot + geom_edge_arc(aes(width = weight), alpha = 0.5) 
  }
  
  graph_plot <- graph_plot + scale_edge_width(range = c(0, 3)) +
    geom_node_point(aes(color=p_stat)) +
    geom_node_point(aes(color=p_stat, size=(1-p_stat)), show.legend = FALSE) +
    scale_size_area(max_size = 10) +
    scale_color_gradient2(mid='darkgreen', high='red', limits=c(0.0, 1.0)) +
    geom_node_text(aes(label = label), repel = TRUE) +
    theme_graph(base_family = "mono")
  graph_plot
}
?continuous_scale

# PUBLIC API (Plotting)
#' @description
#' \code{plotBarplot}
#'
#' \code{plotBarplot} barplot of p-values.
#'
#' @param mulea_relaxed_resuts  data.table in relaxed form.
#' @param selection_vector vector for selecting variables to plot. 
#' 
#' 
#' @title Input/Output Functions
#' @name  InputOutputFunctions
#' @export
#'
#' @return Return plot. 
plotBarplot <- function(mulea_relaxed_resuts, selection_vector=NULL, 
                        categories_names='ontologyId', 
                        statistics_value_colname='ontology_p_stat', 
                        statistics_value_cutoff=0.05) {
  
  mulea_relaxed_resuts <- MulEA:::filterRelaxedResultsForPlotting(
    mulea_relaxed_resuts=mulea_relaxed_resuts,
    statistics_value_colname=statistics_value_colname,
    statistics_value_cutoff=statistics_value_cutoff)
  
  if (is.null(selection_vector)) {
    selection_vector <- 1:nrow(mulea_relaxed_resuts)
  }
  unique_mulea_relaxed_resuts <- unique(
    mulea_relaxed_resuts[selection_vector, c(..categories_names, ..statistics_value_colname)]) 
  unique_mulea_relaxed_resuts <- unique_mulea_relaxed_resuts %>% 
    dplyr::arrange(dplyr::desc((!!as.name(statistics_value_colname))))

  unique_mulea_relaxed_resuts_df <- as.data.frame(unique_mulea_relaxed_resuts)
  unique_mulea_relaxed_resuts_df[, 1] <- factor(unique_mulea_relaxed_resuts_df[[1]], 
                                                levels = unique_mulea_relaxed_resuts_df[[1]])
  mulea_gg_plot <- ggplot(unique_mulea_relaxed_resuts_df, 
                          aes_string(x=categories_names, y=statistics_value_colname, 
                                     fill=statistics_value_colname)) +
    geom_bar(stat="identity") +
    scale_fill_gradient2(mid='darkgreen', high='red', limits=c(0.0, 1.0)) +
    coord_flip() +
    theme_light()
  mulea_gg_plot
}


# PUBLIC API (Plotting)
#' @description
#' \code{plotHeatmap}
#'
#' \code{plotHeatmap} merge model and model relsuts into 
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
plotHeatmap <- function(mulea_relaxed_resuts, 
                        statistics_value_colname='ontology_p_stat', 
                        statistics_value_cutoff=0.05) {
  
  model_with_res_dt_relaxed <- MulEA:::filterRelaxedResultsForPlotting(
    mulea_relaxed_resuts=mulea_relaxed_resuts,
    statistics_value_colname=statistics_value_colname,
    statistics_value_cutoff=statistics_value_cutoff)
  
  model_with_res_dt_relaxed_sort_pval <- model_with_res_dt_relaxed %>% dplyr::arrange(., desc(ontology_p_stat), .by_group = FALSE)
  model_with_res_dt_relaxed_sort_pval[,1] <- factor(model_with_res_dt_relaxed_sort_pval[[1]], 
                                                    levels = unique(model_with_res_dt_relaxed_sort_pval[[1]]))
  model_with_res_dt_relaxed_sort_pval[,2] <- factor(model_with_res_dt_relaxed_sort_pval[[2]], 
                                                    levels = unique(rev(model_with_res_dt_relaxed_sort_pval[[2]])))
  ggplot(model_with_res_dt_relaxed_sort_pval, aes(gen_in_ontology, ontologyId, fill= ontology_p_stat)) + 
    scale_fill_gradient2(mid='darkgreen', high='red', limits=c(0.0, 1.0)) +
    geom_tile() +
    theme_light() +
    theme(axis.text.x = element_text(angle = 90))
}
