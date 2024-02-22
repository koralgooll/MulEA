#' An S4 class to represent a set based tests in mulea.
#'
#' @slot method A method from set based methods to count results. Possible
#' values: "Hypergeometric", "SetBasedEnrichment".
#' @slot gmt A data.frame representing the GMT representation of a model.
#' @slot element_names Data from an experiment to analyse across model, e.g. differentially expressed genes.
#' @slot background_element_names Background set for the test, e.g. all genes in the experiment.
#' @slot p_value_adjustment_method A type of algorithm used to adjust values.
#' Possible values: "eFDR", and all options from p.adjust {stats} documentation.
#' @slot number_of_permutations A number of permutations used in set based
#' enrichment test. Default value is 10000.
#' @slot nthreads Number of processor's threads used in calculations.
#' @return ora object. This object represents set based tests in mulea.
#' @export ora
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
#'tf_gmt <- read_gmt(file = system.file(package="mulea", "extdata", 
#' "Transcription_factor_RegulonDB_Escherichia_coli_GeneSymbol.gmt")) 
#' 
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
#'                 nthreads = 4) 

ora <- setClass(
  "ora",
  slots = list(
    gmt = "data.frame",
    element_names = "character",
    background_element_names = "character",
    p_value_adjustment_method = "character",
    number_of_permutations = "numeric",
    nthreads = "numeric",
    test = "function"
  )
)

setMethod("initialize", "ora",
          function(.Object,
                   gmt = data.frame(),
                   element_names = character(),
                   background_element_names = character(),
                   p_value_adjustment_method = "eFDR",
                   number_of_permutations = 10000,
                   test = NULL,
                   nthreads = 4,
                   ...) {
            adjustMethod <- NULL
            .Object@gmt <- gmt
            .Object@element_names <- element_names
            .Object@background_element_names <- background_element_names
            .Object@p_value_adjustment_method <- p_value_adjustment_method
            .Object@number_of_permutations <- number_of_permutations
            .Object@nthreads <- nthreads
            
            .Object@test <- function(setBasemodel) {
              setBasedTestRes <- NULL
              
              if (!identical(setBasemodel@p_value_adjustment_method, character(0)) &&
                  setBasemodel@p_value_adjustment_method == "eFDR") {
                muleaSetBaseEnrichmentTest <-
                  SetBasedEnrichmentTest(
                    gmt = setBasemodel@gmt,
                    element_names = setBasemodel@element_names,
                    pool = setBasemodel@background_element_names,
                    number_of_permutations = setBasemodel@number_of_permutations,
                    nthreads = setBasemodel@nthreads
                  )
                
                muleaSetBaseEnrichmentTest <-
                  run_test(muleaSetBaseEnrichmentTest)
                
                muleaSetBaseEnrichmentTest <- merge(
                  setBasemodel@gmt[c('ontology_id', 'ontology_name')],
                  muleaSetBaseEnrichmentTest,
                  by.x = "ontology_id",
                  by.y = "DB_names",
                  all = TRUE
                )
                
                for (i in 1:length(muleaSetBaseEnrichmentTest$FDR)) {
                  if (!is.nan(muleaSetBaseEnrichmentTest$FDR[i])
                      && muleaSetBaseEnrichmentTest$FDR[i] > 1.0) {
                    muleaSetBaseEnrichmentTest$FDR[i] <- 1.0e+00
                  }
                }
                
                names(muleaSetBaseEnrichmentTest) <-
                  c(
                    'ontology_id',
                    'ontology_name',
                    'nr_common_with_tested_elements',
                    'nr_common_with_backgound_elements',
                    'Genes_in_DB',
                    'p_value',
                    'P_adj_Bonf',
                    'adjustedPValue',
                    'R_obs',
                    'R_exp',
                    'eFDR'
                  )
                
                setBasedTestRes <-
                  muleaSetBaseEnrichmentTest[, !names(muleaSetBaseEnrichmentTest) %in%
                                               c('Genes_in_DB', 'P_adj_Bonf',
                                                 'R_obs', 'R_exp', 'adjustedPValue')]
              } else {
                MuleaHypergeometricTest <-
                  MuleaHypergeometricTest(
                    gmt = setBasemodel@gmt,
                    element_names = setBasemodel@element_names,
                    pool = setBasemodel@background_element_names,
                    nthreads = setBasemodel@nthreads
                  )
                setBasedTestRes <- run_test(MuleaHypergeometricTest)
                
                muleaSetBaseEnrichmentTest <- merge(
                  setBasemodel@gmt[c('ontology_id', 'ontology_name')],
                  setBasedTestRes,
                  by.x = "ontology_id",
                  by.y = "ontology_name",
                  all = TRUE
                )
                
                names(muleaSetBaseEnrichmentTest) <-
                  c(
                    'ontology_id',
                    'ontology_name',
                    'list_of_values',
                    'p_value'
                  )
                if (!identical(setBasemodel@p_value_adjustment_method, character(0)) &&
                    setBasemodel@p_value_adjustment_method != "eFDR") {
                  muleaSetBaseEnrichmentTest <-
                    data.frame(
                      muleaSetBaseEnrichmentTest,
                      "adjusted_p_value" = stats::p.adjust(muleaSetBaseEnrichmentTest$p_value, method = setBasemodel@p_value_adjustment_method)
                    )
                  setBasedTestRes <-
                    muleaSetBaseEnrichmentTest[, !names(muleaSetBaseEnrichmentTest) %in%
                                                 c('list_of_values')]
                }
              }
              
              setBasedTestRes
            }
            
            .Object
            
          })

#' @describeIn Runs ORA test.
#' @param model Object of S4 class representing the mulea test.
#' @return run_test method for ora object. Returns the results of the overrepresentation analysis.
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

setMethod("run_test",
          signature(model = "ora"),
          function(model) {
            model@test(model)
          })
