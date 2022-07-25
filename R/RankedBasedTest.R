#' An S4 class to represent a ranked based tests in Mulea.
#'
#' @slot gmt A data.frame representing GMT's reprezentation of model.
#' @slot element_names A data from expeciment to analize accross model.
#' @slot element_scores A vectore of element_scores per element_names.
#' @slot p A power of weight. Default value is 1.
#' @slot element_score_type Defines the GSEA score type. Only positive element_scores - "pos",
#' only negative element_scores - "neg" and mixed (standard) - "std".
#' @slot number_of_permutations A number of permutations used in KS test. Default
#' vlue is 1000.
#' @return GSEA object. This object represents ranked based tests in
#' Mulea.
#' @export "GSEA"
#' @examples
#' modelDfFromFile <- MulEA::read_gmt(
#'   file = system.file(package="MulEA", "extdata", "model.gmt"))
#' dataFromExperiment <- c("FBgn0004407", "FBgn0010438", "FBgn0003742",
#'                         "FBgn0029709", "FBgn0030341", "FBgn0037044",
#'                         "FBgn0002887", "FBgn0028434", "FBgn0030170",
#'                         "FBgn0263831")
#' dataFromExperimentScores <- c(0.09, 0.11, 0.15, 0.20, 0.21, 0.24, 0.28, 0.30,
#'                               0.45, 0.50)
#' rankedBasedTestSubramanian <- GSEA(gmt = modelDfFromFile,
#'                                               element_names = dataFromExperiment,
#'                                               element_scores = dataFromExperimentScores)
GSEA <- setClass(
  "GSEA",
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

setMethod("initialize", "GSEA",
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

#' @describeIn GSEA runs test calculations.
#' @param model Object of s4 class represents Mulea Test.
#' @return run_test method for GSEA object. Returns results of
#' counting using methods from ranking based area.
#' @examples
#' rankedBasedTestSubramanianRes <- MulEA::run_test(rankedBasedTestSubramanian)
setMethod("run_test",
          signature(model = "GSEA"),
          function(model) {
            model@test(model)
          })
