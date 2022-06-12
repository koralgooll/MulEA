
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
MuleaHypergeometricTest <- setClass("MuleaHypergeometricTest",
                             slots = list(
                               gmt = "data.frame",
                               testData = "character",
                               pool = "character",
                               test = "function"
                             ))

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

            .Object@test <- function(testObject) {
                testObject@testData <- checkIfPoolIncludeSample(
                  testObject@gmt, testObject@testData, testObject@pool)
                
                muleaSetBaseEnrichmentTest <- SetBasedEnrichmentTest(
                  gmt = testObject@gmt,
                  testData = testObject@testData,
                  pool = testObject@pool,
                  only_hyper_geometric_test=TRUE)
                
                muleaSetBaseEnrichmentTestResult <<- runTest(muleaSetBaseEnrichmentTest)
                testObjectGlobal <<- testObject

                testResults <- data.frame(
                  'ontologyName' = muleaSetBaseEnrichmentTestResult$DB_names,
                  'listOfValues' = testObject@gmt$listOfValues,
                  'p.value' = muleaSetBaseEnrichmentTestResult$P_val, row.names = NULL)
                
                testResults
            }

            .Object

          })

#' @describeIn MuleaHypergeometricTest runs test calculations.
#' @param testObject Object of s4 class represents Mulea Test.
#' @return runTest method for MuleaHypergeometricTest object. Used as private
#' function.
#' @examples
#' \dontrun{
#' #It is a private method. Look at runTest of ORA's examples.
#' }
setMethod("runTest",
          signature(testObject = "MuleaHypergeometricTest"),
          function(testObject) {
            testObject@test(testObject)
          })
