
#' PRIVATE class : An S4 class to represent a ranked based tests in Mulea.
#'
#' @slot gmt A data.frame representing GMT's reprezentation of model.
#' @slot testData A data from expeciment to analize accross model.
#' @slot scores A vectore of scores per testData.
#' @slot p A power of weight.
#' @slot scoreType Defines the GSEA score type. Only positive scores - "pos", only negative scores - "neg" and mixed (standard) - "std".
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
                              scoreType = "character",
                              test = "function"
                            ))

setMethod("initialize", "SubramanianTest",
          function(.Object,
                   gmt = data.frame(),
                   testData = character(),
                   scores = numeric(),
                   p = 1,
                   scoreType = "std",
                   test = NULL,
                   ...) {

            .Object@gmt <- gmt
            .Object@testData <- testData
            .Object@scores <- scores
            .Object@p <- p
            .Object@scoreType <- scoreType

            .Object@test <- function(testObject) {

              listmodelDfFromFile <- testObject@gmt$listOfValues
              names(listmodelDfFromFile) <- testObject@gmt$ontologyId

              samplesToAnalisys <- testObject@scores
              names(samplesToAnalisys) <- testObject@testData

              fgseaRes <- fgsea::fgsea(pathways = listmodelDfFromFile, 
                                       stats = samplesToAnalisys, 
                                       gseaParam = testObject@p,
                                       scoreType = testObject@scoreType)

              resultDf <- merge(testObject@gmt, fgseaRes, by.x = "ontologyId", 
                                by.y = "pathway", all = TRUE)
              resultDf <- plyr::ddply(.data = resultDf, .variables = c('ontologyId'), .fun = function(df_row){
                nrCommonGenesOntologyBackground <- length(df_row[,'listOfValues'][[1]])
                cbind(df_row, nrCommonGenesOntologyBackground)
              })[c("ontologyId", "ontologyName", 'size', 'nrCommonGenesOntologyBackground', "pval", "padj")]
              colnames(resultDf) <- c("ontologyId", "ontologyName", 
                                      'nrCommonGenesOntologySet', 'nrCommonGenesOntologyBackground', 
                                      "pValue", "adjustedPValue")
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
