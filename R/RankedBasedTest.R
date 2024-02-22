#' Gene Set Enrichment Analysis (GSEA)
#' 
#' An S4 class to represent a ranked based tests in mulea.
#'
#' @slot gmt  A data.frame representing the GMT representation of a model.
#' @slot element_names A vector of elements names to include in the analysis,
#' ordered by their scores.
#' @slot element_scores A vector of element_scores per element_names.
#' @slot gsea_power A power of weight. Default value is 1.
#' @slot element_score_type Defines the GSEA score type. Only positive
#' element_scores - "pos", only negative element_scores - "neg" and mixed
#' (standard) - "std".
#' @slot number_of_permutations The number of permutations used in KS test.
#' Default value is 1000.
#' @slot test character
#' @return GSEA object. This object represents ranked based tests.
#' @export
#' @examples
#' library(mulea)
#' library(tidyverse)
#' geo2r_result_tab <- read_tsv("GSE55662.table_wt_non_vs_cipro.tsv")
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
#' tf_gmt <- read_gmt("Transcription_factor_RegulonDB_Escherichia_coli_GeneSymbol.gmt")
#' tf_gmt_filtered <- filter_ontology(gmt = tf_gmt,
#'                           min_nr_of_elements = 3,
#'                           max_nr_of_elements = 400)
#' # creating the GSEA model using the GMT variable
#' gsea_model <- gsea(gmt = tf_gmt_filtered,
#'                    # the names of elements to test
#'                    element_names = geo2r_result_tab_filtered$Gene.symbol,
#'                    # the logFC-s of elements to test
#'                    element_scores = geo2r_result_tab_filtered$logFC,
#'                    # consider elements having positive logFC values only
#'                    element_score_type = "pos",
#'                    # the number of permutations
#'                    number_of_permutations = 10000)

gsea <- setClass(
  "gsea",
  slots = list(
    gmt = "data.frame",
    element_names = "character",
    element_scores = "numeric",
    gsea_power = "numeric",
    element_score_type = "character",
    number_of_permutations = "numeric",
    test = "function"
  )
)

setMethod("initialize", "gsea",
          function(.Object,
                   gmt = data.frame(),
                   element_names = character(),
                   element_scores = numeric(),
                   gsea_power = 1,
                   element_score_type = "std",
                   number_of_permutations = 1000,
                   test = NULL,
                   ...) {
            .Object@gmt <- gmt
            .Object@element_names <- element_names
            .Object@element_scores <- element_scores
            .Object@gsea_power <- gsea_power
            .Object@element_score_type <- element_score_type
            
            .Object@number_of_permutations <- number_of_permutations
            
            .Object@test <- function(rankedBasemodel) {
              rankedTestRes <- NULL
              
              subramanianTest <- SubramanianTest(
                gmt = rankedBasemodel@gmt,
                element_names = rankedBasemodel@element_names,
                element_scores = rankedBasemodel@element_scores,
                gsea_power = rankedBasemodel@gsea_power,
                element_score_type = rankedBasemodel@element_score_type
              )
              rankedTestRes <- run_test(subramanianTest)
              
              rankedTestRes
            }
            
            .Object
            
          })

#' @describeIn gsea runs test calculations.
#' @param model Object of S4 class representing the mulea test.
#' @return run_test method for GSEA object. Returns results of
#' the enrichment analysis.
#' @examples
#' #' library(mulea)
#' library(tidyverse)
#' geo2r_result_tab <- read_tsv("GSE55662.table_wt_non_vs_cipro.tsv")
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
#' tf_gmt <- read_gmt("Transcription_factor_RegulonDB_Escherichia_coli_GeneSymbol.gmt")
#' tf_gmt_filtered <- filter_ontology(gmt = tf_gmt,
#'                           min_nr_of_elements = 3,
#'                           max_nr_of_elements = 400)
#' # creating the GSEA model using the GMT variable
#' gsea_model <- gsea(gmt = tf_gmt_filtered,
#'                    # the names of elements to test
#'                    element_names = geo2r_result_tab_filtered$Gene.symbol,
#'                    # the logFC-s of elements to test
#'                    element_scores = geo2r_result_tab_filtered$logFC,
#'                    # consider elements having positive logFC values only
#'                    element_score_type = "pos",
#'                    # the number of permutations
#'                    number_of_permutations = 10000)
#'                    
#'                    # running the GSEA
#' gsea_results <- run_test(gsea_model)


setMethod("run_test",
          signature(model = "gsea"),
          function(model) {
            model@test(model)
          })
