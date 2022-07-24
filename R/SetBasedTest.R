#' An S4 class to represent a set based tests in Mulea.
#'
#' @slot method A method from set based methods to count results. Possible values: "Hypergeometric", "SetBaseEnrichment".
#' @slot gmt A data.frame representing GMT's reprezentation of model.
#' @slot testData A data from expeciment to analize accross model.
#' @slot pool A background data to count test.
#' @slot adjustMethod A type of algorithm used to adjust values. Possible values: "PT" and all from p.adjust {stats} documentation.
#' @slot numberOfPermutations A number of permutations used in set base enrichment test. Default vlue is 10000.
#' @slot nthreads Number of processor's threads used in calculations.
#' @return ORA object. This object represents set based tests in Mulea.
#' @export "ORA"
#' @examples
#' modelDfFromFile <- MulEA::readGmtFileAsDataFrame(gmtFilePath = system.file(package="MulEA", "extdata", "model.gmt"))
#' dataFromExperiment <- c("FBgn0004407", "FBgn0010438", "FBgn0003742", "FBgn0029709", "FBgn0030341", "FBgn0037044", "FBgn0002887", "FBgn0028434", "FBgn0030170", "FBgn0263831")
#' dataFromExperimentPool <- unique(c(c("FBgn0033690", "FBgn0261618", "FBgn0004407", "FBgn0010438", "FBgn0032154", "FBgn0039930", "FBgn0040268", "FBgn0013674", "FBgn0037008", "FBgn0003116", "FBgn0037743", "FBgn0035401", "FBgn0037044", "FBgn0051005", "FBgn0026737", "FBgn0026751", "FBgn0038704", "FBgn0002887", "FBgn0028434", "FBgn0030170", "FBgn0263831", "FBgn0000579"), c("FBgn0066666", "FBgn0000000", "FBgn0099999", "FBgn0011111", "FBgn0022222", "FBgn0777777", "FBgn0333333", "FBgn0003742", "FBgn0029709", "FBgn0030341")))
#' setBasedTest <- ORA(gmt = modelDfFromFile, testData = dataFromExperiment, 
#'                     nthreads = 2)
#' setBasedTestWithPool <- ORA(gmt = modelDfFromFile, testData = dataFromExperiment, pool = dataFromExperimentPool, nthreads = 2)
#' setBasedTestWithPoolAndAdjust <- ORA(gmt = modelDfFromFile, testData = dataFromExperiment, pool = dataFromExperimentPool, adjustMethod = "BH", nthreads = 2)
#' setBaseTestWithPermutationTestAdjustment <- ORA(gmt = modelDfFromFile, testData = dataFromExperiment, adjustMethod = "PT", nthreads = 2)
ORA <- setClass(
  "ORA",
  slots = list(
    gmt = "data.frame",
    testData = "character",
    pool = "character",
    adjustMethod = "character",
    numberOfPermutations = "numeric",
    nthreads = "numeric",
    test = "function"
  )
)

setMethod("initialize", "ORA",
          function(.Object,
                   gmt = data.frame(),
                   testData = character(),
                   pool = character(),
                   adjustMethod = character(),
                   numberOfPermutations = 10000,
                   test = NULL,
                   nthreads = 4,
                   ...) {
            .Object@gmt <- gmt
            .Object@testData <- testData
            .Object@pool <- pool
            .Object@adjustMethod <- 'PT'
            .Object@numberOfPermutations <- numberOfPermutations
            .Object@nthreads <- nthreads
            
            .Object@test <- function(setBasemodel) {
              setBasedTestRes <- NULL
              
              if (!identical(setBasemodel@adjustMethod, character(0)) &&
                  setBasemodel@adjustMethod == "PT") {
                muleaSetBaseEnrichmentTest <-
                  SetBasedEnrichmentTest(
                    gmt = setBasemodel@gmt,
                    testData = setBasemodel@testData,
                    pool = setBasemodel@pool,
                    numberOfPermutations = setBasemodel@numberOfPermutations,
                    nthreads = setBasemodel@nthreads
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
                    testData = setBasemodel@testData,
                    pool = setBasemodel@pool
                  )
                setBasedTestRes <- run_test(muleaHypergeometricTest)
                if (!identical(setBasemodel@adjustMethod, character(0)) &&
                    setBasemodel@adjustMethod != "PT") {
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

#' @describeIn ORA runs test calculations.
#' @param model Object of s4 class represents Mulea Test.
#' @return run_test method for ORA object. Returns results of counting using methods from set based area.
#' @examples
#' modelDfFromFile <- MulEA::readGmtFileAsDataFrame(gmtFilePath = system.file(package="MulEA", "extdata", "model.gmt"))
#' dataFromExperiment <- c("FBgn0004407", "FBgn0010438", "FBgn0003742", "FBgn0029709", "FBgn0030341", "FBgn0037044", "FBgn0002887", "FBgn0028434", "FBgn0030170", "FBgn0263831")
#' dataFromExperimentPool <- unique(c(c("FBgn0033690", "FBgn0261618", "FBgn0004407", "FBgn0010438", "FBgn0032154", "FBgn0039930", "FBgn0040268", "FBgn0013674", "FBgn0037008", "FBgn0003116", "FBgn0037743", "FBgn0035401", "FBgn0037044", "FBgn0051005", "FBgn0026737", "FBgn0026751", "FBgn0038704", "FBgn0002887", "FBgn0028434", "FBgn0030170", "FBgn0263831", "FBgn0000579"), c("FBgn0066666", "FBgn0000000", "FBgn0099999", "FBgn0011111", "FBgn0022222", "FBgn0777777", "FBgn0333333", "FBgn0003742", "FBgn0029709", "FBgn0030341")))
#' setBasedTest <- ORA(gmt = modelDfFromFile, testData = dataFromExperiment, 
#'                     nthreads = 2)
#' setBasedTestWithPool <- ORA(gmt = modelDfFromFile, 
#'                             testData = dataFromExperiment, 
#'                             pool = dataFromExperimentPool, nthreads = 2)
#' setBasedTestWithPoolAndAdjust <- ORA(gmt = modelDfFromFile, 
#'                                      testData = dataFromExperiment, 
#'                                      pool = dataFromExperimentPool, 
#'                                      adjustMethod = "BH", nthreads = 2)
#' setBaseTestWithPermutationTestAdjustment <- ORA(gmt = modelDfFromFile, testData = dataFromExperiment, adjustMethod = "PT", nthreads = 2)
#' setBasedTestRes <- MulEA::run_test(setBasedTest)
#' setBasedTestWithPoolRes <- MulEA::run_test(setBasedTestWithPool)
#' setBasedTestWithPoolAndAdjustRes <- MulEA::run_test(setBasedTestWithPoolAndAdjust)
#' setBaseTestWithPermutationTestAdjustmentRes <- MulEA::run_test(setBaseTestWithPermutationTestAdjustment)
setMethod("run_test",
          signature(model = "ORA"),
          function(model) {
            model@test(model)
          })
