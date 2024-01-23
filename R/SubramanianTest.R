#' PRIVATE class : An S4 class to represent a ranked based tests in mulea.
#'
#' @slot gmt A data.frame representing GMT's reprezentation of model.
#' @slot element_names A data from expeciment to analize accross model.
#' @slot element_scores A vectore of element_scores per element_names.
#' @slot p A power of weight.
#' @slot element_score_type Defines the GSEA score type. Only positive element_scores - "pos", only negative element_scores - "neg" and mixed (standard) - "std".
#' @return dataframe with presented columns 'ontology_id', 'ontology_name',
#' 'nrCommonGenesOntologySet', 'nrCommonGenesOntologyBackground',
#' 'pValue', 'adjustedPValue'
#' @examples
#' \dontrun{
#' #It is a private s4 object. Look at GSAE's examples.
#' }
SubramanianTest <- setClass(
  "SubramanianTest",
  slots = list(
    gmt = "data.frame",
    element_names = "character",
    element_scores = "numeric",
    gsea_power = "numeric",
    element_score_type = "character",
    test = "function"
  )
)

#' @importFrom fgsea fgsea
setMethod("initialize", "SubramanianTest",
          function(.Object,
                   gmt = data.frame(),
                   element_names = character(),
                   element_scores = numeric(),
                   gsea_power = 1,
                   element_score_type = "std",
                   test = NULL,
                   ...) {
            .Object@gmt <- gmt
            .Object@element_names <- element_names
            .Object@element_scores <- element_scores
            .Object@gsea_power <- gsea_power
            .Object@element_score_type <- element_score_type
            
            .Object@test <- function(model) {
              listmodelDfFromFile <- model@gmt$list_of_values
              names(listmodelDfFromFile) <-
                model@gmt$ontology_id
              
              samplesToAnalisys <- model@element_scores
              names(samplesToAnalisys) <- model@element_names

              fgseaRes <-
                fgsea::fgsea(
                  pathways = listmodelDfFromFile,
                  stats = samplesToAnalisys,
                  gseaParam = model@gsea_power,
                  scoreType = model@element_score_type
                )
              
              resultDf <-
                merge(
                  model@gmt,
                  fgseaRes,
                  by.x = "ontology_id",
                  by.y = "pathway",
                  all = TRUE
                )
              resultDf <-
                plyr::ddply(
                  .data = resultDf,
                  .variables = c('ontology_id'),
                  .fun = function(df_row) {
                    nrCommonGenesOntologyBackground <-
                      length(df_row[, 'list_of_values'][[1]])
                    cbind(df_row, nrCommonGenesOntologyBackground)
                  }
                )[c(
                  "ontology_id",
                  "ontology_name",
                  'size',
                  'nrCommonGenesOntologyBackground',
                  "pval",
                  "padj"
                )]
              colnames(resultDf) <- c(
                "ontology_id",
                "ontology_name",
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
#' @param model Object of s4 class represents mulea Test.
#' @return run_test method for SubramanianTest object. Used as private function.
#' @examples
#' \dontrun{
#' #It is a private method. Look at run_test of GSEA's examples.
#' }
setMethod("run_test",
          signature(model = "SubramanianTest"),
          function(model) {
            model@test(model)
          })
