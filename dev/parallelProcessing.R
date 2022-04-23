library(MulEA)
muleaPkgDir <- find.package("MulEA")

modelDfFromFile <- MulEA::readGmtFileAsDataFrame(
  gmtFilePath = paste(muleaPkgDir,"/extdata/model.gmt", sep = ""))

selectDf <- read.csv2(file = './inst/extdata/selectData.csv')
select <- selectDf[['select']]
selectScores <- selectDf[['score']]

poolDf <- read.csv2(file = './inst/extdata/poolData.csv')
pool <- poolDf[['pool']]

number_of_steps <- 1000



mulea_ora_model <- MulEA::ORA(
  gmt = modelDfFromFile, testData = select, 
  pool = pool, adjustMethod = "PT",
  numberOfPermutations = number_of_steps, nthreads = 8)
mulea_ora_results <- MulEA::runTest(mulea_ora_model)

mulea_ora_model_1 <- MulEA::ORA(
  gmt = modelDfFromFile, testData = select, 
  pool = pool, adjustMethod = "PT",
  numberOfPermutations = number_of_steps, nthreads = 1)

mulea_ora_results_1 <- MulEA::runTest(mulea_ora_model_1)

mulea_ora_reshaped_results <- MulEA::reshapeResults(
  mulea_model=mulea_ora_model, 
  mulea_model_resuts=mulea_ora_results, 
  category_stat_column_name='adjustedPValueEmpirical')


mulea_ranked_model <- MulEA::RankedBasedTest(
  gmt = modelDfFromFile, testData = select, scores = selectScores)

mulea_sub_results <- MulEA::runTest(mulea_ranked_model)

mulea_sub_reshaped_results <- MulEA::reshapeResults(
  mulea_model = mulea_ranked_model, 
  mulea_model_resuts = mulea_sub_results, 
  mulea_model_resuts_ontology_col_name='ontologyId')


MulEA::plotGraph(mulea_relaxed_resuts=mulea_ora_reshaped_results, statistics_value_cutoff = 1.00, 
                 statistics_value_colname = 'adjustedPValueEmpirical')

MulEA::plotBarplot(mulea_relaxed_resuts = mulea_ora_reshaped_results, statistics_value_cutoff=1.00,
                   statistics_value_colname = 'adjustedPValueEmpirical')

MulEA::plotHeatmap(mulea_relaxed_resuts=mulea_ora_reshaped_results, statistics_value_cutoff=1.00,
                   statistics_value_colname = 'adjustedPValueEmpirical')


MulEA::plotGraph(mulea_relaxed_resuts=mulea_sub_reshaped_results, statistics_value_cutoff = 1.00)

MulEA::plotBarplot(mulea_relaxed_resuts = mulea_sub_reshaped_results, statistics_value_cutoff=1.00)

MulEA::plotHeatmap(mulea_relaxed_resuts=mulea_sub_reshaped_results, statistics_value_cutoff=1.00)




DB <- "Dupa"
library(parallel)
nthread=4
cl <- makeCluster(spec=nthread, type = "PSOCK",outfile= "log1.txt")
clusterExport(cl,  "DB")
clusterExport(cl,"list_of_all_genes")
clusterExport(cl,"pool")

clusterExport(cl,"select")
clusterExport(cl,"enrichment_test_simulation")
clusterExport(cl,"seeds_per_thread")
clusterExport(cl,"steps_per_thread")
clusterEvalQ(cl,library(Rcpp))
clusterEvalQ(cl,Sys.setenv("PKG_CXXFLAGS"="-std=c++11"))
clusterEvalQ(cl,sourceCpp("src/set-based-enrichment-test.cpp"))


