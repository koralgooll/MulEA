# install.packages("BiocManager")
# BiocManager::install("fgsea")
# install.packages("DBI")
# install.packages("RSQLite")
# install.packages('plyr')

# install.packages("roxygen2")
# install.packages("tidygraph")
# install.packages("ggraph")


# install.packages("RCurl")
# install.packages("devtools")
# BiocManager::install("enrichplot")
# BiocManager::install("clusterProfiler")
# extra data -> BiocManager::install("breastCancerMAINZ")


library(MulEA)

muleaPkgDir <- find.package("MulEA")
modelDfFromFile <- MulEA::readGmtFileAsDataFrame(gmtFilePath = paste(muleaPkgDir,"/extdata/model.gmt", sep = ""))
dataFromExperiment <- c("FBgn0004407", "FBgn0010438", "FBgn0003742", "FBgn0029709", "FBgn0030341", "FBgn0037044", "FBgn0002887", "FBgn0028434", "FBgn0030170", "FBgn0263831", "FBgn0261618", "FBgn0038704", "FBgn0000579")
dataFromExperimentScores <- c(0.09, 0.11, 0.15, 0.20, 0.21, 0.24, 0.28, 0.30, 0.45, 0.50, 0.53, 0.60, 0.61)
dataFromExperimentPool <- unique(c(c("FBgn0033690", "FBgn0261618", "FBgn0004407", "FBgn0010438", "FBgn0032154", "FBgn0039930", "FBgn0040268", "FBgn0013674",
                                     "FBgn0037008", "FBgn0003116", "FBgn0037743", "FBgn0035401", "FBgn0037044", "FBgn0051005", "FBgn0026737", "FBgn0026751",
                                     "FBgn0038704", "FBgn0002887", "FBgn0028434", "FBgn0030170", "FBgn0263831", "FBgn0000579"),
                                   c("FBgn0066666", "FBgn0000000", "FBgn0099999", "FBgn0011111", "FBgn0022222", "FBgn0777777", "FBgn0333333", "FBgn0003742",
                                     "FBgn0029709", "FBgn0030341")))

setBasedTestWithPool <- ORA(gmt = modelDfFromFile, testData = dataFromExperiment,
                                              pool = dataFromExperimentPool)
setBasedTestWithPoolRes <- MulEA::runTest(setBasedTestWithPool)


setBasedTestWithPoolAndAdjust <- ORA(gmt = modelDfFromFile, testData = dataFromExperiment,
                                              pool = dataFromExperimentPool, adjustMethod = "BH")
setBasedTestWithPoolAndAdjustRes <- MulEA::runTest(setBasedTestWithPoolAndAdjust)


setBasedTestWithPoolAndAdjust <- ORA(gmt = modelDfFromFile, testData = dataFromExperiment,
                                              pool = dataFromExperimentPool, adjustMethod = "PT",
                                              numberOfPermutations = 100)
setBasedTestWithPoolAndAdjustResNew <- MulEA::runTest(setBasedTestWithPoolAndAdjust)



pathToModelGmtFile <- '/home/cezary/Downloads/Caenorhabditis_elegans_interactions_SS_GMT_ncbiGeneID.gmt'
KnowledgeBaseDf <- MulEA::readGmtFileAsDataFrame(gmtFilePath = pathToModelGmtFile)






# Read inputs
input_gmt <- MulEA::readGmtFileAsDataFrame('D:/projects/science/GSEA_4.1.0/inputs/KEGG_Paths2geneSymbols.gmt')
input_rank <- read.table(file='D:/projects/science/GSEA_4.1.0/inputs/lung_FoldChange.rnk', sep='\t')
input_select <- input_rank[['V1']]
input_select_scores <- input_rank[['V2']]

# Setup params
number_of_steps <- 10

# Perform ORA
# IMPORTANT : Scores are very different based on scoreType param.
mulea_ora_model <- MulEA::ORA(
  gmt = input_gmt, testData = input_select, adjustMethod = "PT",
  numberOfPermutations = number_of_steps)
mulea_ora_results <- MulEA::runTest(mulea_ora_model)


mulea_ora_model@testData

mulea_ora_model@pool <- unique(unlist(mulea_ora_model@gmt[, 'listOfValues']))

mulea_ora_model@pool

length(mulea_ora_model@testData)
length(mulea_ora_model@pool)


all(mulea_ora_model@testData %in% mulea_ora_model@pool)

mulea_ora_model@testData <- mulea_ora_model@testData[mulea_ora_model@testData %in% mulea_ora_model@pool]
length(mulea_ora_model@testData)

length(mulea_ora_model@testData %in% mulea_ora_model@pool)

all(c('a', 'b', 'e') %in% c('a', 'b', 'c', 'd'))



rstudioapi::writeRStudioPreference(name = "console_max_lines", value = as.integer('10000'))