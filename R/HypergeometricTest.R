#' PRIVATE class : An S4 class to represent a Hypergeometric tests in Mulea.
#'
#' @slot gmt A data.frame representing GMT's reprezentation of model.
#' @slot testData A data from expeciment to analize accross model.
#' @slot pool A background data to count test.
#' @return MuleaHypergeometricTest object. Used as private function.
#' @examples
#' \dontrun{
#' #It is a private s4 object. Look at ORA's examples.
#' }
MuleaHypergeometricTest <- setClass(
  "MuleaHypergeometricTest",
  slots = list(
    gmt = "data.frame",
    testData = "character",
    pool = "character",
    test = "function"
  )
)

setMethod("initialize", "MuleaHypergeometricTest",
          function(.Object,
                   gmt = data.frame(),
                   testData = character(),
                   pool = character(),
                   test = NULL,
                   ...) {
            .Object@gmt <- gmt
            .Object@testData <- testData
            .Object@pool <- pool
            
            .Object@test <- function(model) {
              model@testData <- checkIfPoolIncludeSample(model@gmt, model@testData, model@pool)
              
              muleaSetBaseEnrichmentTest <-
                SetBasedEnrichmentTest(
                  gmt = model@gmt,
                  testData = model@testData,
                  pool = model@pool,
                  only_hyper_geometric_test = TRUE
                )
              
              muleaSetBaseEnrichmentTestResult <<-
                run_test(muleaSetBaseEnrichmentTest)
              modelGlobal <<- model
              
              
              testResults <- data.frame(
                'ontologyName' = muleaSetBaseEnrichmentTestResult$DB_names,
                'listOfValues' = model@gmt$listOfValues,
                'p.value' = muleaSetBaseEnrichmentTestResult$P_val,
                row.names = NULL
              )
              
              testResults
            }
            
            .Object
            
          })

#' @describeIn MuleaHypergeometricTest runs test calculations.
#' @param model Object of s4 class represents Mulea Test.
#' @return run_test method for MuleaHypergeometricTest object. Used as private
#' function.
#' @examples
#' \dontrun{
#' #It is a private method. Look at run_test of ORA's examples.
#' }
setMethod("run_test",
          signature(model = "MuleaHypergeometricTest"),
          function(model) {
            model@test(model)
          })
