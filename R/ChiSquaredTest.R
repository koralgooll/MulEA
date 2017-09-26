
setClass("muleaChiSquaredTest",
         slots = list(
           testData = "character",
           test = "function"
         ))

setMethod("initialize", "muleaChiSquaredTest",
          function(.Object,
                   testData = character(),
                   test = NULL,
                   ...) {

            .Object@testData <- testData

            .Object@test <- function(dataObject, testObject) {
              MulEA::calculateChiSquaredTest(model = dataObject@gmt, sampleVector = testObject@testData)
            }

            .Object

          })

