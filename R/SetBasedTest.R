
#' An S4 class to represent a set based tests in Mulea.
#'
#' @slot gmt A data.frame representing GMT's reprezentation of model.
#' @slot testData A data from expeciment to analize accross model.
#' @slot pool A background data to count test.
#' @slot adjustMethod A type of algorithm used to adjust values.
#' @return SetBasedTest object. This object represents set based tests in Mulea.
#' @examples
#' modelDfFromFile <- MulEA::readGmtFileAsDataFrame(gmtFilePath = system.file(package="MulEA", "extdata", "model.gmt"))
#' dataFromExperiment <- c("FBgn0004407", "FBgn0010438", "FBgn0003742", "FBgn0029709", "FBgn0030341", "FBgn0037044", "FBgn0002887", "FBgn0028434", "FBgn0030170", "FBgn0263831")
#' dataFromExperimentPool <- unique(c(c("FBgn0033690", "FBgn0261618", "FBgn0004407", "FBgn0010438", "FBgn0032154", "FBgn0039930", "FBgn0040268", "FBgn0013674", "FBgn0037008", "FBgn0003116", "FBgn0037743", "FBgn0035401", "FBgn0037044", "FBgn0051005", "FBgn0026737", "FBgn0026751", "FBgn0038704", "FBgn0002887", "FBgn0028434", "FBgn0030170", "FBgn0263831", "FBgn0000579"), c("FBgn0066666", "FBgn0000000", "FBgn0099999", "FBgn0011111", "FBgn0022222", "FBgn0777777", "FBgn0333333", "FBgn0003742", "FBgn0029709", "FBgn0030341")))
#' setBasedTest <- SetBasedTest(gmt = modelDfFromFile, testData = dataFromExperiment)
#' setBasedTestWithPool <- SetBasedTest(gmt = modelDfFromFile, testData = dataFromExperiment, pool = dataFromExperimentPool)
#' setBasedTestWithPoolAndAdjust <- SetBasedTest(gmt = modelDfFromFile, testData = dataFromExperiment, pool = dataFromExperimentPool, adjustMethod = "BH")
SetBasedTest <- setClass("SetBasedTest",
                                    slots = list(
                                      gmt = "data.frame",
                                      testData = "character",
                                      pool = "character",
                                      adjustMethod = "character",
                                      test = "function"
                                    ))

setMethod("initialize", "SetBasedTest",
          function(.Object,
                   gmt = data.frame(),
                   testData = character(),
                   pool = character(),
                   adjustMethod = character(),
                   test = NULL,
                   ...) {

            .Object@gmt <- gmt
            .Object@testData <- testData
            .Object@pool <- pool
            .Object@adjustMethod <- adjustMethod

            .Object@test <- function(setBaseTestObject) {
              muleaHypergeometricTest <- MuleaHypergeometricTest(gmt = setBaseTestObject@gmt, testData = setBaseTestObject@testData, pool = setBaseTestObject@pool)
              muleaHypergeometricTestRes <- runTest(muleaHypergeometricTest)
              if (length(setBaseTestObject@adjustMethod) != 0) {
                if (setBaseTestObject@adjustMethod == "EszterMethod") {
                  # adjustResult <- data.frame(modelWithTestsResults, gseaPermutationTest(modelWithTestDf = modelWithTestsResults, steps = steps, sampleVector = sampleVector, poolVector = poolVector))
                } else {
                  muleaHypergeometricTestRes <- data.frame(muleaHypergeometricTestRes, "q.value" = p.adjust(muleaHypergeometricTestRes$p.value, method = adjustMethod))
                }
              }
              muleaHypergeometricTestRes
            }

            .Object

          })

#' @describeIn SetBasedTest runs test calculations.
#' @param testObject Object of s4 class represents Mulea Test.
#' @return runTest method for SetBasedTest object. Returns results of counting using methods from set based area.
#' @examples
#' modelDfFromFile <- MulEA::readGmtFileAsDataFrame(gmtFilePath = system.file(package="MulEA", "extdata", "model.gmt"))
#' dataFromExperiment <- c("FBgn0004407", "FBgn0010438", "FBgn0003742", "FBgn0029709", "FBgn0030341", "FBgn0037044", "FBgn0002887", "FBgn0028434", "FBgn0030170", "FBgn0263831")
#' dataFromExperimentPool <- unique(c(c("FBgn0033690", "FBgn0261618", "FBgn0004407", "FBgn0010438", "FBgn0032154", "FBgn0039930", "FBgn0040268", "FBgn0013674", "FBgn0037008", "FBgn0003116", "FBgn0037743", "FBgn0035401", "FBgn0037044", "FBgn0051005", "FBgn0026737", "FBgn0026751", "FBgn0038704", "FBgn0002887", "FBgn0028434", "FBgn0030170", "FBgn0263831", "FBgn0000579"), c("FBgn0066666", "FBgn0000000", "FBgn0099999", "FBgn0011111", "FBgn0022222", "FBgn0777777", "FBgn0333333", "FBgn0003742", "FBgn0029709", "FBgn0030341")))
#' setBasedTest <- SetBasedTest(gmt = modelDfFromFile, testData = dataFromExperiment)
#' setBasedTestWithPool <- SetBasedTest(gmt = modelDfFromFile, testData = dataFromExperiment, pool = dataFromExperimentPool)
#' setBasedTestWithPoolAndAdjust <- SetBasedTest(gmt = modelDfFromFile, testData = dataFromExperiment, pool = dataFromExperimentPool, adjustMethod = "BH")
#' setBasedTestRes <- MulEA::runTest(setBasedTest)
#' setBasedTestWithPoolRes <- MulEA::runTest(setBasedTestWithPool)
#' setBasedTestWithPoolAndAdjustRes <- MulEA::runTest(setBasedTestWithPoolAndAdjust)
setMethod("runTest",
          signature(testObject = "SetBasedTest"),
          function(testObject) {
            testObject@test(testObject)
          })
