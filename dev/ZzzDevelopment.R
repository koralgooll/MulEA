# BiocManager::install(version = "3.7")
# BiocManager::install("fgsea", version = "3.7", update = FALSE)
# install.packages("RCurl")
# install.packages("roxygen2")


muleaPkgDir <- find.package("MulEA")
modelDfFromFile <- MulEA::readGmtFileAsDataFrame(gmtFilePath = paste(muleaPkgDir,"/extdata/model.gmt", sep = ""))
dataFromExperiment <- c("FBgn0004407", "FBgn0010438", "FBgn0003742", "FBgn0029709", "FBgn0030341", "FBgn0037044", "FBgn0002887", "FBgn0028434", "FBgn0030170", "FBgn0263831", "FBgn0261618", "FBgn0038704", "FBgn0000579")
dataFromExperimentScores <- c(0.09, 0.11, 0.15, 0.20, 0.21, 0.24, 0.28, 0.30, 0.45, 0.50, 0.53, 0.60, 0.61)
dataFromExperimentPool <- unique(c(c("FBgn0033690", "FBgn0261618", "FBgn0004407", "FBgn0010438", "FBgn0032154", "FBgn0039930", "FBgn0040268", "FBgn0013674",
                                     "FBgn0037008", "FBgn0003116", "FBgn0037743", "FBgn0035401", "FBgn0037044", "FBgn0051005", "FBgn0026737", "FBgn0026751",
                                     "FBgn0038704", "FBgn0002887", "FBgn0028434", "FBgn0030170", "FBgn0263831", "FBgn0000579"),
                                   c("FBgn0066666", "FBgn0000000", "FBgn0099999", "FBgn0011111", "FBgn0022222", "FBgn0777777", "FBgn0333333", "FBgn0003742",
                                     "FBgn0029709", "FBgn0030341")))

library(MulEA)

setBasedTestWithPoolAndAdjust <- SetBasedTest(gmt = modelDfFromFile, testData = dataFromExperiment,
                                              pool = dataFromExperimentPool, adjustMethod = "BH")
setBasedTestWithPoolAndAdjustRes <- MulEA::runTest(setBasedTestWithPoolAndAdjust)

setBasedTestWithPoolAndAdjust <- SetBasedTest(gmt = modelDfFromFile, testData = dataFromExperiment,
                                              pool = dataFromExperimentPool, adjustMethod = "PT")
setBasedTestWithPoolAndAdjustRes <- MulEA::runTest(setBasedTestWithPoolAndAdjust)





