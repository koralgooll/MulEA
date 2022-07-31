#' An S4 class to represent a set based tests in Mulea.
#'
#' @slot method A method from set based methods to count results. Possible values: "Hypergeometric", "SetBaseEnrichment".
#' @slot gmt A data.frame representing GMT's reprezentation of model.
#' @slot element_names A data from expeciment to analize accross model.
#' @slot background_element_names A background data to count test.
#' @slot p_value_adjustment_method A type of algorithm used to adjust values. Possible values: "PT" and all from p.adjust {stats} documentation.
#' @slot number_of_permutations A number of permutations used in set base enrichment test. Default vlue is 10000.
#' @slot number_of_cpu_threads Number of processor's threads used in calculations.
#' @return ora object. This object represents set based tests in Mulea.
#' @export ora
#' @examples
#' modelDfFromFile <- MulEA::read_gmt(file = system.file(package="MulEA", "extdata", "model.gmt"))
#' dataFromExperiment <- c("FBgn0004407", "FBgn0010438", "FBgn0003742", "FBgn0029709", "FBgn0030341", "FBgn0037044", "FBgn0002887", "FBgn0028434", "FBgn0030170", "FBgn0263831")
#' dataFromExperimentPool <- unique(c(c("FBgn0033690", "FBgn0261618", "FBgn0004407", "FBgn0010438", "FBgn0032154", "FBgn0039930", "FBgn0040268", "FBgn0013674", "FBgn0037008", "FBgn0003116", "FBgn0037743", "FBgn0035401", "FBgn0037044", "FBgn0051005", "FBgn0026737", "FBgn0026751", "FBgn0038704", "FBgn0002887", "FBgn0028434", "FBgn0030170", "FBgn0263831", "FBgn0000579"), c("FBgn0066666", "FBgn0000000", "FBgn0099999", "FBgn0011111", "FBgn0022222", "FBgn0777777", "FBgn0333333", "FBgn0003742", "FBgn0029709", "FBgn0030341")))
#' setBasedTest <- ora(gmt = modelDfFromFile, element_names = dataFromExperiment, 
#'                     number_of_cpu_threads = 2)
#' setBasedTestWithPool <- ora(gmt = modelDfFromFile, element_names = dataFromExperiment, background_element_names = dataFromExperimentPool, number_of_cpu_threads = 2)
#' setBasedTestWithPoolAndAdjust <- ora(gmt = modelDfFromFile, element_names = dataFromExperiment, background_element_names = dataFromExperimentPool, p_value_adjustment_method = "BH", number_of_cpu_threads = 2)
#' setBaseTestWithPermutationTestAdjustment <- ora(gmt = modelDfFromFile, element_names = dataFromExperiment, p_value_adjustment_method = "PT", number_of_cpu_threads = 2)
ora <- setClass(
  "ora",
  slots = list(
    gmt = "data.frame",
    element_names = "character",
    background_element_names = "character",
    p_value_adjustment_method = "character",
    number_of_permutations = "numeric",
    number_of_cpu_threads = "numeric",
    test = "function"
  )
)

setMethod("initialize", "ora",
          function(.Object,
                   gmt = data.frame(),
                   element_names = character(),
                   background_element_names = character(),
                   p_value_adjustment_method = character(),
                   number_of_permutations = 10000,
                   test = NULL,
                   number_of_cpu_threads = 4,
                   ...) {
            .Object@gmt <- gmt
            .Object@element_names <- element_names
            .Object@background_element_names <- background_element_names
            .Object@p_value_adjustment_method <- 'PT'
            .Object@number_of_permutations <- number_of_permutations
            .Object@number_of_cpu_threads <- number_of_cpu_threads
            
            .Object@test <- function(setBasemodel) {
              setBasedTestRes <- NULL
              
              if (!identical(setBasemodel@p_value_adjustment_method, character(0)) &&
                  setBasemodel@p_value_adjustment_method == "PT") {
                muleaSetBaseEnrichmentTest <-
                  SetBasedEnrichmentTest(
                    gmt = setBasemodel@gmt,
                    element_names = setBasemodel@element_names,
                    pool = setBasemodel@background_element_names,
                    number_of_permutations = setBasemodel@number_of_permutations,
                    nthreads = setBasemodel@number_of_cpu_threads
                  )
                
                muleaSetBaseEnrichmentTest <-
                  run_test(muleaSetBaseEnrichmentTest)
                
                muleaSetBaseEnrichmentTest <- merge(
                  setBasemodel@gmt[c('ontologyId', 'ontologyName')],
                  muleaSetBaseEnrichmentTest,
                  by.x = "ontologyId",
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
                    'ontologyId',
                    'ontologyName',
                    'nrCommonGenesOntologySet',
                    'nrCommonGenesOntologyBackground',
                    'Genes_in_DB',
                    'pValue',
                    'P_adj_Bonf',
                    'adjustedPValue',
                    'R_obs',
                    'R_exp',
                    'adjustedPValueEmpirical'
                  )
                
                setBasedTestRes <-
                  muleaSetBaseEnrichmentTest[, !names(muleaSetBaseEnrichmentTest) %in%
                                               c('Genes_in_DB', 'P_adj_Bonf',
                                                 'R_obs', 'R_exp')]
              } else {
                muleaHypergeometricTest <-
                  MuleaHypergeometricTest(
                    gmt = setBasemodel@gmt,
                    element_names = setBasemodel@element_names,
                    pool = setBasemodel@pool
                  )
                setBasedTestRes <- run_test(muleaHypergeometricTest)
                if (!identical(setBasemodel@p_value_adjustment_method, character(0)) &&
                    setBasemodel@p_value_adjustment_method != "PT") {
                  setBasedTestRes <-
                    data.frame(
                      setBasedTestRes,
                      "q.value" = p.adjust(setBasedTestRes$p.value, method = adjustMethod)
                    )
                }
              }
              
              setBasedTestRes
            }
            
            .Object
            
          })

#' @describeIn ora runs test calculations.
#' @param model Object of s4 class represents Mulea Test.
#' @return run_test method for ora object. Returns results of counting using methods from set based area.
#' @examples
#' modelDfFromFile <- MulEA::read_gmt(file = system.file(package="MulEA", "extdata", "model.gmt"))
#' dataFromExperiment <- c("FBgn0004407", "FBgn0010438", "FBgn0003742", "FBgn0029709", "FBgn0030341", "FBgn0037044", "FBgn0002887", "FBgn0028434", "FBgn0030170", "FBgn0263831")
#' dataFromExperimentPool <- unique(c(c("FBgn0033690", "FBgn0261618", "FBgn0004407", "FBgn0010438", "FBgn0032154", "FBgn0039930", "FBgn0040268", "FBgn0013674", "FBgn0037008", "FBgn0003116", "FBgn0037743", "FBgn0035401", "FBgn0037044", "FBgn0051005", "FBgn0026737", "FBgn0026751", "FBgn0038704", "FBgn0002887", "FBgn0028434", "FBgn0030170", "FBgn0263831", "FBgn0000579"), c("FBgn0066666", "FBgn0000000", "FBgn0099999", "FBgn0011111", "FBgn0022222", "FBgn0777777", "FBgn0333333", "FBgn0003742", "FBgn0029709", "FBgn0030341")))
#' setBasedTest <- ora(gmt = modelDfFromFile, element_names = dataFromExperiment, 
#'                     number_of_cpu_threads = 2)
#' setBasedTestWithPool <- ora(gmt = modelDfFromFile, 
#'                             element_names = dataFromExperiment, 
#'                             background_element_names = dataFromExperimentPool, number_of_cpu_threads = 2)
#' setBasedTestWithPoolAndAdjust <- ora(gmt = modelDfFromFile, 
#'                                      element_names = dataFromExperiment, 
#'                                      background_element_names = dataFromExperimentPool, 
#'                                      p_value_adjustment_method = "BH", number_of_cpu_threads = 2)
#' setBaseTestWithPermutationTestAdjustment <- ora(gmt = modelDfFromFile, element_names = dataFromExperiment, p_value_adjustment_method = "PT", number_of_cpu_threads = 2)
#' setBasedTestRes <- MulEA::run_test(setBasedTest)
#' setBasedTestWithPoolRes <- MulEA::run_test(setBasedTestWithPool)
#' setBasedTestWithPoolAndAdjustRes <- MulEA::run_test(setBasedTestWithPoolAndAdjust)
#' setBaseTestWithPermutationTestAdjustmentRes <- MulEA::run_test(setBaseTestWithPermutationTestAdjustment)
setMethod("run_test",
          signature(model = "ora"),
          function(model) {
            model@test(model)
          })
