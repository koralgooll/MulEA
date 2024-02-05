#' Gene Set Enrichment Analysis (GSEA)
#' 
#' An S4 class to represent a ranked based tests in mulea.
#'
#' @slot gmt A data frame which contains the data, imported from a GMT file.
#' @slot element_names A vector of elements names to include in the analysis,
#' ordered by their scores.
#' @slot element_scores A vector of element_scores per element_names.
#' @slot gsea_power A power of weight. Default value is 1.
#' @slot element_score_type Defines the GSEA score type. Only positive
#' element_scores - "pos", only negative element_scores - "neg" and mixed
#' (standard) - "std".
#' @slot number_of_permutations The number of permutations used in KS test.
#' Default value is 1000.
#' @slot test character
#' @return GSEA object. This object represents ranked based tests.
#' @export
#' @examples
#' modelDfFromFile <- read_gmt(
#'   file = system.file(package="mulea", "extdata", "model.gmt"))
#' dataFromExperiment <- c("FBgn0004407", "FBgn0010438", "FBgn0003742",
#'                         "FBgn0029709", "FBgn0030341", "FBgn0037044",
#'                         "FBgn0002887", "FBgn0028434", "FBgn0030170",
#'                         "FBgn0263831")
#' dataFromExperimentScores <- c(0.09, 0.11, 0.15, 0.20, 0.21, 0.24, 0.28, 0.30,
#'                               0.45, 0.50)
#' GSEASubramanian <- gsea(gmt = modelDfFromFile,
#'                         element_names = dataFromExperiment,
#'                         element_scores = dataFromExperimentScores)
gsea <- setClass(
  "gsea",
  slots = list(
    gmt = "data.frame",
    element_names = "character",
    element_scores = "numeric",
    gsea_power = "numeric",
    element_score_type = "character",
    number_of_permutations = "numeric",
    test = "function"
  )
)

setMethod("initialize", "gsea",
          function(.Object,
                   gmt = data.frame(),
                   element_names = character(),
                   element_scores = numeric(),
                   gsea_power = 1,
                   element_score_type = "std",
                   number_of_permutations = 1000,
                   test = NULL,
                   ...) {
            .Object@gmt <- gmt
            .Object@element_names <- element_names
            .Object@element_scores <- element_scores
            .Object@gsea_power <- gsea_power
            .Object@element_score_type <- element_score_type
            
            .Object@number_of_permutations <- number_of_permutations
            
            .Object@test <- function(rankedBasemodel) {
              rankedTestRes <- NULL
              
              subramanianTest <- SubramanianTest(
                gmt = rankedBasemodel@gmt,
                element_names = rankedBasemodel@element_names,
                element_scores = rankedBasemodel@element_scores,
                gsea_power = rankedBasemodel@gsea_power,
                element_score_type = rankedBasemodel@element_score_type
              )
              rankedTestRes <- run_test(subramanianTest)
              
              rankedTestRes
            }
            
            .Object
            
          })

#' @describeIn gsea runs test calculations.
#' @param model Object of s4 class represents mulea Test.
#' @return run_test method for GSEA object. Returns results of
#' counting using methods from ranking based area.
#' @examples
#' GSEASubramanianRes <- run_test(GSEASubramanian)
setMethod("run_test",
          signature(model = "gsea"),
          function(model) {
            model@test(model)
          })
