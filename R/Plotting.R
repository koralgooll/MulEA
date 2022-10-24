validate_column_names_and_function_args <- function(data, ...) {
  arguments <- list(...)
  if (!all(unlist(arguments) %in% names(data))) {
    stop('Wrongly set data column names.')
  }
}

filterRelaxedResultsForPlotting <- function(reshaped_results,
                                            p_value_type_colname = 'ontologyStatValue',
                                            p_value_max_threshold = 0.05) {
  include <- !is.na(reshaped_results[[p_value_type_colname]])
  reshaped_results_filtered_na <-
    reshaped_results[include,]
  include <-
    reshaped_results_filtered_na[[p_value_type_colname]] <= p_value_max_threshold
  reshaped_results_filtered_cutoff <-
    reshaped_results_filtered_na[include, ]
  reshaped_results_filtered_cutoff
}

#' Reshape Results
#' 
#' This function merges model and model results into a single data frame.
#'
#' @param model a MulEA model, created e.g. by ora().
#' @param model_results Results from model, returned by run_test().
#' @param model_ontology_col_name character
#' @param ontology_id_colname character
#' @param p_value_type_colname character
#' @param p_value_max_threshold logical
#' @seealso \code{\link{plot_graph}}, \code{\link{plot_barplot}},
#' \code{\link{plot_heatmap}}
#' @importFrom data.table :=
#' @export
#'
#' @return Return detailed and relaxed datatable where model and results are
#' merged for plotting purposes.
reshape_results <-
  function(model = NULL,
           model_results = NULL,
           model_ontology_col_name = 'ontologyId',
           ontology_id_colname = 'ontologyId',
           p_value_type_colname = 'adjustedPValue',
           p_value_max_threshold = TRUE) {
    genIdInOntology <- NULL
    model_with_res <-
      merge(
        x = model@gmt,
        y = model_results,
        by.x = model_ontology_col_name,
        by.y = ontology_id_colname,
        all = TRUE
      )
    model_with_res_dt <- data.table::setDT(model_with_res)
    model_with_res_dt_size = 0
    for (i in 1:nrow(model_with_res_dt)) {
      model_with_res_dt_size <-
        model_with_res_dt_size + length(model_with_res_dt[[i, 'listOfValues']])
    }
    model_with_res_dt_relaxed <- data.table::data.table(
      ontologyId = rep('a', length.out = model_with_res_dt_size),
      genIdInOntology = rep('a', length.out = model_with_res_dt_size),
      ontologyStatValue = rep(1.0, length.out = model_with_res_dt_size)
    )
    
    model_with_res_dt_relaxed_counter = 1
    for (i in 1:nrow(model_with_res_dt)) {
      category_name <- model_with_res_dt[[i, 'ontologyId']]
      category_p_stat <-
        model_with_res_dt[[i, p_value_type_colname]]
      for (item_name in model_with_res_dt[[i, 'listOfValues']]) {
        # TODO: THE LINE BELOW DOES NOT UPDATE THE OBJECT IS THIS INTENTIONAL?
        model_with_res_dt_relaxed[model_with_res_dt_relaxed_counter,
                                  c("ontologyId", "genIdInOntology", "ontologyStatValue") :=
                                    list(category_name, item_name, category_p_stat)]
        model_with_res_dt_relaxed_counter = model_with_res_dt_relaxed_counter + 1
      }
    }
    if (p_value_max_threshold) {
      model_with_res_dt_relaxed <-
        model_with_res_dt_relaxed[genIdInOntology %in% model@element_names]
    }
    names(model_with_res_dt_relaxed) <-
      c("ontologyId", "genIdInOntology", p_value_type_colname)
    model_with_res_dt_relaxed
  }


#' Plot Graph (Network)
#' 
#' Plots graph representation of enrichment results.
#'
#' @param reshaped_results data.table in relaxed form.
#' @param shared_elements_min_threshold numeric
#' @param p_value_type_colname character
#' @param ontology_id_colname character
#' @param ontology_element_colname numeric
#' @param p_value_max_threshold numeric
#' @return Return a graph.
#' @importFrom data.table :=
#' @importFrom rlang .data
#' @seealso \code{\link{reshape_results}}
#' @export
plot_graph <- function(reshaped_results,
                      shared_elements_min_threshold = 0,
                      p_value_type_colname = 'adjustedPValue',
                      ontology_id_colname = 'ontologyId',
                      ontology_element_colname = 'genIdInOntology',
                      p_value_max_threshold = 0.05) {
  ontologyId <- NULL
  edges <- NULL
  validate_column_names_and_function_args(
    data = reshaped_results,
    p_value_type_colname,
    ontology_id_colname,
    ontology_element_colname
  )
  reshaped_results <- data.table::setDT(reshaped_results)
  model_with_res_dt_relaxed <-
    filterRelaxedResultsForPlotting(
      reshaped_results = reshaped_results,
      p_value_type_colname = p_value_type_colname,
      p_value_max_threshold = p_value_max_threshold
    )
  
  ontologies <-
    unique(model_with_res_dt_relaxed[, ontology_id_colname, with = FALSE])
  ontologies_graph_edges_num <- sum(1:(nrow(ontologies) - 1))
  ontologies_graph_edges <- data.table::data.table(
    from = rep('a', length.out = ontologies_graph_edges_num),
    to = rep('a', length.out = ontologies_graph_edges_num),
    weight = rep(0, length.out = ontologies_graph_edges_num)
  )
  
  if (0 == ontologies_graph_edges_num) {
    stop('No edges at all. Wrong data.table or manipulate p_value_max_threshold please.')
  }
  
  ontologies_graph_edges_counter <- 1
  for (i in 1:(nrow(ontologies) - 1)) {
    ontology_name_i <- ontologies[i, ontologyId]
    genes_in_ontology_i <-
      model_with_res_dt_relaxed[ontologyId == ontology_name_i, ][[ontology_element_colname]]
    
    for (j in (i + 1):nrow(ontologies)) {
      ontology_name_j <- ontologies[j, ontologyId]
      genes_in_ontology_j <- model_with_res_dt_relaxed[ontologyId == ontology_name_j, ][[ontology_element_colname]]
      genes_in_ontology_i_j_intersection_num <-
        length(intersect(genes_in_ontology_i, genes_in_ontology_j))
      if (shared_elements_min_threshold < genes_in_ontology_i_j_intersection_num) {
        ontologies_graph_edges[ontologies_graph_edges_counter,
                               c('from', 'to', 'weight') := list(ontology_name_i,
                                                                 ontology_name_j,
                                                                 genes_in_ontology_i_j_intersection_num)]
        ontologies_graph_edges_counter = ontologies_graph_edges_counter + 1
      }
    }
  }
  ontologies_graph_edges <-
    ontologies_graph_edges[0:(ontologies_graph_edges_counter - 1), ]
  
  nodes_ids <- model_with_res_dt_relaxed[, ontologyId]
  nodes_p_stat <-
    model_with_res_dt_relaxed[[p_value_type_colname]]
  ontologies_graph_nodes <- data.table::data.table(id = nodes_ids,
                                                   label = nodes_ids,
                                                   p_stat = nodes_p_stat)
  
  ontologies_graph_nodes <- unique(ontologies_graph_nodes)
  
  routes_tidy <-
    tidygraph::tbl_graph(nodes = ontologies_graph_nodes,
                         edges = ontologies_graph_edges,
                         directed = TRUE)
  
  graph_plot <-
    ggraph::ggraph(routes_tidy, layout = "linear", circular = TRUE)
  
  
  if (0 != nrow(tibble::as_tibble(tidygraph::activate(routes_tidy, edges)))) {
    graph_plot <-
      graph_plot + ggraph::geom_edge_arc(aes(width = .data$weight), alpha = 0.5)
  }
  
  graph_plot <- graph_plot + ggraph::scale_edge_width(range = c(0, 3)) +
    ggraph::geom_node_point(aes(color = .data$p_stat)) +
    ggraph::geom_node_point(aes(color = .data$p_stat, size = (1 - .data$p_stat)), show.legend = FALSE) +
    scale_size_area(max_size = 10) +
    scale_color_gradient2(
      mid = 'darkgreen',
      high = 'red',
      limits = c(0.0, 1.0),
      name = p_value_type_colname
    ) +
    ggraph::geom_node_text(aes(label = .data$label), repel = TRUE) +
    ggraph::theme_graph(base_family = "mono")
  graph_plot
}


#' Plot Barplot
#' 
#' Plots barplot of p-values.
#'
#' @param reshaped_results  data.table in relaxed form.
#' @param selected_rows_to_plot numeric; which rows of the reshaped results data
#' frame should be included in the plot?
#' @param ontology_id_colname character
#' @param p_value_type_colname character
#' @param p_caue_max_threshold numeric
#' @importFrom magrittr %>%
#' @import ggplot2
#' @seealso \code{\link{reshape_results}}
#' @export
#'
#' @return Return a barplot.
plot_barplot <-
  function(reshaped_results,
           selected_rows_to_plot = NULL,
           ontology_id_colname = 'ontologyId',
           p_value_type_colname = 'adjustedPValue',
           p_value_max_threshold = 0.05) {
    validate_column_names_and_function_args(data = reshaped_results,
                                                    p_value_type_colname, ontology_id_colname)
    reshaped_results <- filterRelaxedResultsForPlotting(
      reshaped_results = reshaped_results,
      p_value_type_colname = p_value_type_colname,
      p_value_max_threshold = p_value_max_threshold
    )
    
    if (is.null(selected_rows_to_plot)) {
      selected_rows_to_plot <- 1:nrow(reshaped_results)
    }
    unique_reshaped_results <-
      unique(reshaped_results[selected_rows_to_plot, c(ontology_id_colname, p_value_type_colname), with = FALSE])
    unique_reshaped_results <- unique_reshaped_results %>%
      dplyr::arrange(dplyr::desc((!!as.name(
        p_value_type_colname
      ))))
    
    unique_reshaped_results_df <-
      as.data.frame(unique_reshaped_results)
    unique_reshaped_results_df[, 1] <-
      factor(unique_reshaped_results_df[[1]],
             levels = unique_reshaped_results_df[[1]])
    mulea_gg_plot <- ggplot(
      unique_reshaped_results_df,
      aes_string(x = ontology_id_colname, y = p_value_type_colname,
                 fill = p_value_type_colname)
    ) +
      geom_bar(stat = "identity") +
      scale_fill_gradient2(mid = 'darkgreen',
                           high = 'red',
                           limits = c(0.0, 1.0)) +
      coord_flip() +
      theme_light()
    mulea_gg_plot
  }


#' Plot Heatmap
#' 
#' Plots heatmap of enriched terms and obtained p-values.
#'
#' @param reshaped_results data.table in relaxed form.
#' @param p_value_type_colname character
#' @param ontology_element_colname character
#' @param p_value_max_threshold numeric
#' @importFrom magrittr %>%
#' @import ggplot2
#' @seealso \code{\link{reshape_results}}
#' @export
#'
#' @return Return a heatmap.
plot_heatmap <- function(reshaped_results,
                        p_value_type_colname = 'adjustedPValue',
                        ontology_element_colname = 'genIdInOntology',
                        p_value_max_threshold = 0.05) {
  validate_column_names_and_function_args(data = reshaped_results,
                                                  p_value_type_colname,
                                                  ontology_element_colname)
  model_with_res_dt_relaxed <-
    filterRelaxedResultsForPlotting(
      reshaped_results = reshaped_results,
      p_value_type_colname = p_value_type_colname,
      p_value_max_threshold = p_value_max_threshold
    )
  
  model_with_res_dt_relaxed_sort_pval <-
    model_with_res_dt_relaxed %>%
    dplyr::arrange(dplyr::desc((
      !!rlang::sym(p_value_type_colname)
    )), .by_group = FALSE)
  model_with_res_dt_relaxed_sort_pval[, 1] <-
    factor(
      model_with_res_dt_relaxed_sort_pval[[1]],
      levels = unique(model_with_res_dt_relaxed_sort_pval[[1]])
    )
  model_with_res_dt_relaxed_sort_pval[, 2] <-
    factor(model_with_res_dt_relaxed_sort_pval[[2]],
           levels = unique(rev(model_with_res_dt_relaxed_sort_pval[[2]])))
  ggplot(
    model_with_res_dt_relaxed_sort_pval,
    aes(
      !!rlang::sym(ontology_element_colname),
      .data$ontologyId,
      fill = !!rlang::sym(p_value_type_colname)
    )
  ) +
    scale_fill_gradient2(mid = 'darkgreen',
                         high = 'red',
                         limits = c(0.0, 1.0)) +
    geom_tile() +
    theme_light() +
    theme(axis.text.x = element_text(angle = 90))
}
