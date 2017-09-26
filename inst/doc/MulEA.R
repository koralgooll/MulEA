## ----global_options, include=TRUE, echo=FALSE----------------------------
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE, error = TRUE)

## ---- results = 'asis'---------------------------------------------------
pathToModelGmtFile <- paste(find.package("MulEA"),"/example/model.gmt", sep = "")
modelDfFromFile <- MulEA::readGmtFileAsDF(gmtFilePath = pathToModelGmtFile)
knitr::kable(modelDfFromFile, caption = "Model Data Frame")

## ---- eval=FALSE---------------------------------------------------------
#     MulEA::saveModelFromDataFrameToGmtFile(modelDF = modelDfFromFile, gmtFilePath = pathToModelGmtFile)

## ---- eval=FALSE---------------------------------------------------------
#    MulEA::startLocalDatabase(":memory:")
#    MulEA::addModelToLocalDatabase(model = modelDfFromFile, taxonomy_id = 9001, model_source = "GO", version = 0);
#    modelDfFromLocalDB <- MulEA::getModelFromLocalDatabaseAsDf(taxonomy_id = 9001, model_source = "GO", version = 0)
#    modelListFromLocalDB <- MulEA::getModelFromLocalDatabaseAsList(taxonomy_id = 9001, model_source = "GO", version = 0)
#    MulEA::stopLocalDatabase()

## ---- results = 'asis', echo = TRUE--------------------------------------
muleaPkgDir <- find.package("MulEA")
modelDfFromFile <- MulEA::readGmtFileAsDF(gmtFilePath = paste(muleaPkgDir,"/example/model.gmt", sep = ""))
dataFromExperiment <- c("FBgn0004407", "FBgn0010438", "FBgn0003742", "FBgn0029709", "FBgn0030341", "FBgn0037044", "FBgn0002887", "FBgn0028434", "FBgn0030170", "FBgn0263831")
dataFromExperimentPool <- unique(c(c("FBgn0033690", "FBgn0261618", "FBgn0004407", "FBgn0010438", "FBgn0032154", "FBgn0039930", "FBgn0040268", "FBgn0013674",
                                   "FBgn0037008", "FBgn0003116", "FBgn0037743", "FBgn0035401", "FBgn0037044", "FBgn0051005", "FBgn0026737", "FBgn0026751",
                                   "FBgn0038704", "FBgn0002887", "FBgn0028434", "FBgn0030170", "FBgn0263831", "FBgn0000579"),
                                 c("FBgn0066666", "FBgn0000000", "FBgn0099999", "FBgn0011111", "FBgn0022222", "FBgn0777777", "FBgn0333333", "FBgn0003742",
                                   "FBgn0029709", "FBgn0030341")))

## ---- results = 'asis', echo = TRUE--------------------------------------
setBasedTest <- SetBasedTest(gmt = modelDfFromFile, testData = dataFromExperiment)
setBasedTestRes <- MulEA::runTest(setBasedTest)

## ---- results = 'asis'---------------------------------------------------
knitr::kable(setBasedTestRes, caption = "Set Based Test Result Data Frame")

## ---- results = 'asis', echo = TRUE--------------------------------------
setBasedTestWithPool <- SetBasedTest(gmt = modelDfFromFile, testData = dataFromExperiment, pool = dataFromExperimentPool)
setBasedTestWithPoolRes <- MulEA::runTest(setBasedTestWithPool)

## ---- results = 'asis'---------------------------------------------------
knitr::kable(setBasedTestWithPoolRes, caption = "Set Based Test Result Data Frame")

## ---- results = 'asis', echo = TRUE--------------------------------------
setBasedTestWithPoolAndAdjust <- SetBasedTest(gmt = modelDfFromFile, testData = dataFromExperiment, pool = dataFromExperimentPool, adjustMethod = "BH")
setBasedTestWithPoolAndAdjustRes <- MulEA::runTest(setBasedTestWithPoolAndAdjust)

## ---- results = 'asis'---------------------------------------------------
knitr::kable(setBasedTestWithPoolAndAdjustRes, caption = "Set Based Test Result Data Frame")

## ---- results = 'asis', echo = TRUE--------------------------------------
muleaPkgDir <- find.package("MulEA")
modelDfFromFile <- MulEA::readGmtFileAsDF(gmtFilePath = paste(muleaPkgDir,"/example/model.gmt", sep = ""))
dataFromExperiment <- c("FBgn0004407", "FBgn0010438", "FBgn0003742", "FBgn0029709", "FBgn0030341", "FBgn0037044", "FBgn0002887", "FBgn0028434", "FBgn0030170", "FBgn0263831")
dataFromExperimentScores <- c(0.09, 0.11, 0.15, 0.20, 0.21, 0.24, 0.28, 0.30, 0.45, 0.50)

## ---- results = 'asis', echo=TRUE----------------------------------------
rankedBasedTestKs <- RankedBasedTest(method = "KS", gmt = modelDfFromFile, testData = dataFromExperiment)
rankedBasedTestKsRes <- MulEA::runTest(rankedBasedTestKs)

## ---- results = 'asis', echo=TRUE----------------------------------------
rankedBasedTestSubramanian <- RankedBasedTest(method = "Subramanian", gmt = modelDfFromFile, testData = dataFromExperiment, scores = dataFromExperimentScores)
rankedBasedTestSubramanianRes <- MulEA::runTest(rankedBasedTestSubramanian)

## ---- results = 'asis'---------------------------------------------------
knitr::kable(rankedBasedTestKsRes, caption = "Ranked Based Test Result Data Frame")

