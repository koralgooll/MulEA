#' PRIVATE class : An S4 class to represent a Hypergeometric tests in Mulea.
#'
#' @slot gmt A data.frame representing GMT's reprezentation of model.
#' @slot element_names A data from expeciment to analize accross model.
#' @slot background_element_names A background data to count test.
#' @return MuleaHypergeometricTest object. Used as private function.
#' @examples
#' \dontrun{
#' #It is a private s4 object. Look at ora's examples.
#' }
MuleaHypergeometricTest <- setClass(
  "MuleaHypergeometricTest",
  slots = list(
    gmt = "data.frame",
    element_names = "character",
    pool = "character",
    test = "function"
  )
)

setMethod("initialize", "MuleaHypergeometricTest",
          function(.Object,
                   gmt = data.frame(),
                   element_names = character(),
                   pool = character(),
                   test = NULL,
                   ...) {
            .Object@gmt <- gmt
            .Object@element_names <- element_names
            .Object@pool <- pool
            
            .Object@test <- function(model) {
              model@element_names <- checkIfPoolIncludeSample(model@gmt, model@element_names, model@pool)
              
              muleaSetBaseEnrichmentTest <-
                SetBasedEnrichmentTest(
                  gmt = model@gmt,
                  element_names = model@element_names,
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
#' #It is a private method. Look at run_test of ora's examples.
#' }
setMethod("run_test",
          signature(model = "MuleaHypergeometricTest"),
          function(model) {
            model@test(model)
          })
