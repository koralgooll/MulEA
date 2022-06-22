#' An S4 class to represent a ranked based tests in Mulea.
#'
#' @slot gmt A data.frame representing GMT's reprezentation of model.
#' @slot testData A data from expeciment to analize accross model.
#' @slot scores A vectore of scores per testData.
#' @slot p A power of weight. Default value is 1.
#' @slot scoreType Defines the GSEA score type. Only positive scores - "pos",
#' only negative scores - "neg" and mixed (standard) - "std".
#' @slot numberOfPermutations A number of permutations used in KS test. Default
#' vlue is 1000.
#' @return RankedBasedTest object. This object represents ranked based tests in
#' Mulea.
#' @export "RankedBasedTest"
#' @examples
#' modelDfFromFile <- MulEA::readGmtFileAsDataFrame(
#'   gmtFilePath = system.file(package="MulEA", "extdata", "model.gmt"))
#' dataFromExperiment <- c("FBgn0004407", "FBgn0010438", "FBgn0003742",
#'                         "FBgn0029709", "FBgn0030341", "FBgn0037044",
#'                         "FBgn0002887", "FBgn0028434", "FBgn0030170",
#'                         "FBgn0263831")
#' dataFromExperimentScores <- c(0.09, 0.11, 0.15, 0.20, 0.21, 0.24, 0.28, 0.30,
#'                               0.45, 0.50)
#' rankedBasedTestSubramanian <- RankedBasedTest(gmt = modelDfFromFile,
#'                                               testData = dataFromExperiment,
#'                                               scores = dataFromExperimentScores)
RankedBasedTest <- setClass(
  "RankedBasedTest",
  slots = list(
    gmt = "data.frame",
    testData = "character",
    scores = "numeric",
    p = "numeric",
    scoreType = "character",
    numberOfPermutations = "numeric",
    test = "function"
  )
)

setMethod("initialize", "RankedBasedTest",
          function(.Object,
                   gmt = data.frame(),
                   testData = character(),
                   scores = numeric(),
                   p = 1,
                   scoreType = "std",
                   numberOfPermutations = 1000,
                   test = NULL,
                   ...) {
            .Object@gmt <- gmt
            .Object@testData <- testData
            .Object@scores <- scores
            .Object@p <- p
            .Object@scoreType <- scoreType
            
            .Object@numberOfPermutations <- numberOfPermutations
            
            .Object@test <- function(rankedBaseTestObject) {
              rankedTestRes <- NULL
              
              subramanianTest <- SubramanianTest(
                gmt = rankedBaseTestObject@gmt,
                testData = rankedBaseTestObject@testData,
                scores = rankedBaseTestObject@scores,
                p = rankedBaseTestObject@p,
                scoreType = rankedBaseTestObject@scoreType
              )
              rankedTestRes <- runTest(subramanianTest)
              
              rankedTestRes
            }
            
            .Object
            
          })

#' @describeIn RankedBasedTest runs test calculations.
#' @param testObject Object of s4 class represents Mulea Test.
#' @return runTest method for RankedBasedTest object. Returns results of
#' counting using methods from ranking based area.
#' @examples
#' rankedBasedTestSubramanianRes <- MulEA::runTest(rankedBasedTestSubramanian)
setMethod("runTest",
          signature(testObject = "RankedBasedTest"),
          function(testObject) {
            testObject@test(testObject)
          })
