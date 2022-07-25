#' PRIVATE class : An S4 class to represent a ranked based tests in Mulea.
#'
#' @slot gmt A data.frame representing GMT's reprezentation of model.
#' @slot testData A data from expeciment to analize accross model.
#' @slot element_scores A vectore of element_scores per testData.
#' @slot p A power of weight.
#' @slot element_score_type Defines the GSEA score type. Only positive element_scores - "pos", only negative element_scores - "neg" and mixed (standard) - "std".
#' @return dataframe with presented columns 'ontologyId', 'ontologyName',
#' 'nrCommonGenesOntologySet', 'nrCommonGenesOntologyBackground',
#' 'pValue', 'adjustedPValue'
#' @examples
#' \dontrun{
#' #It is a private s4 object. Look at RankedBasedTest's examples.
#' }
SubramanianTest <- setClass(
  "SubramanianTest",
  slots = list(
    gmt = "data.frame",
    testData = "character",
    element_scores = "numeric",
    p = "numeric",
    element_score_type = "character",
    test = "function"
  )
)

setMethod("initialize", "SubramanianTest",
          function(.Object,
                   gmt = data.frame(),
                   testData = character(),
                   element_scores = numeric(),
                   p = 1,
                   element_score_type = "std",
                   test = NULL,
                   ...) {
            .Object@gmt <- gmt
            .Object@testData <- testData
            .Object@element_scores <- element_scores
            .Object@p <- p
            .Object@element_score_type <- element_score_type
            
            .Object@test <- function(model) {
              listmodelDfFromFile <- model@gmt$listOfValues
              names(listmodelDfFromFile) <-
                model@gmt$ontologyId
              
              samplesToAnalisys <- model@element_scores
              names(samplesToAnalisys) <- model@testData
              
              fgseaRes <-
                fgsea::fgsea(
                  pathways = listmodelDfFromFile,
                  stats = samplesToAnalisys,
                  gseaParam = model@p,
                  scoreType = model@element_score_type
                )
              
              resultDf <-
                merge(
                  model@gmt,
                  fgseaRes,
                  by.x = "ontologyId",
                  by.y = "pathway",
                  all = TRUE
                )
              resultDf <-
                plyr::ddply(
                  .data = resultDf,
                  .variables = c('ontologyId'),
                  .fun = function(df_row) {
                    nrCommonGenesOntologyBackground <-
                      length(df_row[, 'listOfValues'][[1]])
                    cbind(df_row, nrCommonGenesOntologyBackground)
                  }
                )[c(
                  "ontologyId",
                  "ontologyName",
                  'size',
                  'nrCommonGenesOntologyBackground',
                  "pval",
                  "padj"
                )]
              colnames(resultDf) <- c(
                "ontologyId",
                "ontologyName",
                'nrCommonGenesOntologySet',
                'nrCommonGenesOntologyBackground',
                "pValue",
                "adjustedPValue"
              )
              resultDf
            }
            
            .Object
            
          })

#' @describeIn SubramanianTest runs test calculations.
#' @param model Object of s4 class represents Mulea Test.
#' @return run_test method for SubramanianTest object. Used as private function.
#' @examples
#' \dontrun{
#' #It is a private method. Look at run_test of RankedBasedTest's examples.
#' }
setMethod("run_test",
          signature(model = "SubramanianTest"),
          function(model) {
            model@test(model)
          })
