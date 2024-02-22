#' Run enrichment analysis procedure
#' 
#' This is a generic function that chooses an enrichment analysis procedure
#' based on the model class and runs the analysis.
#' @param model An S4 object which represents one of mulea's tests (ORA or GSEA). See details
#' for more information.
#' @details The function requires the definition of a model. Models currently
#' implemented in mulea include Gene Set Enrichment Analysis (GSEA) and
#' Over-Representation Analysis (ORA). These models must be defined through
#' their specific functions which are provided in this package. 
#' @seealso \code{\link{gsea}}, \code{\link{ora}}
#' @export
#' @examples
#' library(mulea)
#' library(tidyverse)
#' geo2r_result_tab <- read_tsv(file = system.file(package="mulea", "extdata", "GSE55662.table_wt_non_vs_cipro.tsv"))
#' geo2r_result_tab %<>% 
#' # extracting the first gene symbol from the Gene.symbol column
#' mutate(Gene.symbol = str_remove(string = Gene.symbol,
#'                                 pattern = "\\/.*")) %>% 
#'  # removing rows where Gene.symbol is NA
#'  filter(!is.na(Gene.symbol)) %>% 
#'  # ordering by logFC
#'  arrange(desc(logFC))
#'  
#'  sign_genes <- geo2r_result_tab %>% 
#' # filtering for adjusted p-value < 0.05 and logFC > 1
#' filter(adj.P.Val < 0.05 & logFC > 1) %>% 
#'  # selecting the Gene.symbol column
#'  select(Gene.symbol) %>% 
#'  # converting the tibble to a vector
#'  pull() %>% 
#'  # removing duplicates
#'  unique()
#'  
#'  background_genes <- geo2r_result_tab %>% 
#' # selecting the Gene.symbol column
#' select(Gene.symbol) %>% 
#'  # convertin the tibble to a vector
#'  pull() %>% 
#'  # removing duplicates
#'  unique()
#'  
#' tf_gmt <- read_gmt(file = system.file(package="mulea", "extdata", "Transcription_factor_RegulonDB_Escherichia_coli_GeneSymbol.gmt"))
#' tf_gmt_filtered <- filter_ontology(gmt = tf_gmt,
#'                           min_nr_of_elements = 3,
#'                           max_nr_of_elements = 400)
#'
#' # creating the ORA model using the GMT variable
#' ora_model <- ora(gmt = tf_gmt_filtered, 
#'                 # the test set variable
#'                 element_names = sign_genes, 
#'                 # the background set variable
#'                 background_element_names = background_genes, 
#'                 # the p-value adjustment method
#'                 p_value_adjustment_method = "eFDR", 
#'                 # the number of permutations
#'                 number_of_permutations = 10000,
#'                 # the number of processor threads to use
#'                  nthreads = 4)
#' # running the ORA
#' ora_results <- run_test(ora_model)
#' 
#' @return Results in form of data frame. Structure of data frame depends on
#' object processed by this generic method.
#' @importFrom methods new
setGeneric("run_test", function(model)
  standardGeneric("run_test"))
