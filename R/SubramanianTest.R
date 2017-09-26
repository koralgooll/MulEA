
#' PRIVATE class : An S4 class to represent a ranked based tests in Mulea.
#'
#' @slot gmt A data.frame representing GMT's reprezentation of model.
#' @slot testData A data from expeciment to analize accross model.
#' @slot scores A vectore of scores per testData.
#' @slot p A power of weight.
#' @slot numberOfPermutations A number of permutations used in renked test.
#' @return KolmogorovSmirnovTest object. Used as private object.
#' @examples
#' \dontrun{
#' #It is a private s4 object. Look at RankedBasedTest's examples.
#' }
SubramanianTest <- setClass("SubramanianTest",
                            slots = list(
                              gmt = "data.frame",
                              testData = "character",
                              scores = "numeric",
                              p = "numeric",
                              numberOfPermutations = "numeric",
                              test = "function"
                            ))

setMethod("initialize", "SubramanianTest",
          function(.Object,
                   gmt = data.frame(),
                   testData = character(),
                   scores = numeric(),
                   p = 1,
                   numberOfPermutations = 1000,
                   test = NULL,
                   ...) {

            .Object@gmt <- gmt
            .Object@testData <- testData
            .Object@scores <- scores
            .Object@p <- p
            .Object@numberOfPermutations <- numberOfPermutations

            .Object@test <- function(testObject) {

              listmodelDfFromFile <- testObject@gmt$listOfValues
              names(listmodelDfFromFile) <- testObject@gmt$ontologyId

              samplesToAnalisys <- testObject@scores
              names(samplesToAnalisys) <- testObject@testData

              fgseaRes <- fgsea(pathways = listmodelDfFromFile,
                                stats = samplesToAnalisys,
                                gseaParam = testObject@p, nperm = testObject@numberOfPermutations)

              resultDf <- merge(testObject@gmt, fgseaRes, by.x = "ontologyId", by.y = "pathway", all = TRUE)[c("ontologyId", "ontologyName", "listOfValues", "pval")]
              colnames(resultDf) <- c("ontologyId", "ontologyName", "listOfValues", "p.value")
              resultDf
            }

            .Object

          })

#' @describeIn SubramanianTest runs test calculations.
#' @param testObject Object of s4 class represents Mulea Test.
#' @return runTest method for SubramanianTest object. Used as private function.
#' @examples
#' \dontrun{
#' #It is a private method. Look at runTest of RankedBasedTest's examples.
#' }
setMethod("runTest",
          signature(testObject = "SubramanianTest"),
          function(testObject) {
            testObject@test(testObject)
          })
