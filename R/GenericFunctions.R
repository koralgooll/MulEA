

#' Method num.
#' @name runTest
#' @rdname runTest-methods
#' @param testObject Object of s4 class represents one of Mulea's Tests.
#' @export
#' @examples
#' modelDfFromFile <- readGmtFileAsDataFrame(
#'   gmtFilePath = system.file(package="MulEA", "extdata", "model.gmt"))
#' dataFromExperiment <- c(
#'   "FBgn0004407", "FBgn0010438", "FBgn0003742", "FBgn0029709", "FBgn0030341",
#'   "FBgn0037044", "FBgn0002887", "FBgn0028434", "FBgn0030170", "FBgn0263831")
#' dataFromExperimentPool <- unique(
#'   c(c("FBgn0033690", "FBgn0261618", "FBgn0004407", "FBgn0010438",
#'       "FBgn0032154", "FBgn0039930", "FBgn0040268", "FBgn0013674",
#'       "FBgn0037008", "FBgn0003116", "FBgn0037743", "FBgn0035401",
#'       "FBgn0037044", "FBgn0051005", "FBgn0026737", "FBgn0026751",
#'       "FBgn0038704", "FBgn0002887", "FBgn0028434", "FBgn0030170",
#'       "FBgn0263831", "FBgn0000579"),
#'    c("FBgn0066666", "FBgn0000000", "FBgn0099999", "FBgn0011111",
#'      "FBgn0022222", "FBgn0777777", "FBgn0333333", "FBgn0003742",
#'      "FBgn0029709", "FBgn0030341")))
#' setBasedTest <- ORA(gmt = modelDfFromFile, testData = dataFromExperiment)
#' setBasedTestWithPool <- ORA(gmt = modelDfFromFile,
#'                             testData = dataFromExperiment,
#'                             pool = dataFromExperimentPool)
#' setBasedTestWithPoolAndAdjust <- ORA(gmt = modelDfFromFile,
#'                                      testData = dataFromExperiment,
#'                                      pool = dataFromExperimentPool,
#'                                      adjustMethod = "BH")
#' setBasedTestRes <- runTest(setBasedTest)
#' setBasedTestWithPoolRes <- runTest(setBasedTestWithPool)
#' setBasedTestWithPoolAndAdjustRes <- runTest(setBasedTestWithPoolAndAdjust)
#' dataFromExperimentScores <- c(0.09, 0.11, 0.15, 0.20, 0.21, 0.24, 0.28,
#'                               0.30, 0.45, 0.50)
#' rankedBasedTestKs <- RankedBasedTest(method = "KS",
#'                                      gmt = modelDfFromFile,
#'                                      testData = dataFromExperiment)
#' rankedBasedTestSubramanian <- RankedBasedTest(
#'   method = "Subramanian",
#'   gmt = modelDfFromFile,
#'   testData = dataFromExperiment,
#'   scores = dataFromExperimentScores)
#' rankedBasedTestKsRes <- MulEA::runTest(rankedBasedTestKs)
#' rankedBasedTestSubramanianRes <- MulEA::runTest(rankedBasedTestSubramanian)
#' @return Results in form of data frame. Structure of data frame depends on
#' object processed by this generic method.
setGeneric("runTest", function(testObject)
  standardGeneric("runTest"))
