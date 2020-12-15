rm(list = ls())
root_of_cpp_project <- 'D:/projects/science/enrichment-analysis'
setwd(root_of_cpp_project)


convert_model_list_to_model_df <- function(model_list) {
  model_df <- plyr::adply(.data = names(DB), .margins=1, .fun = function(list_name){
    
    data.frame('ontologyId' = list_name, 'ontologyName' = list_name, 
               'listOfValues' = I(DB[list_name]), stringsAsFactors = FALSE)
  })[c('ontologyId', 'ontologyName', 'listOfValues')]
  model_df
}


# Test Case 01
source("src/set-based-enrichment-test.R") # It contains the R source of enrichment analizis and it compiles the C++ part.
load(file="test-parameters-01.RData") # It loas a dummy databaseinto variabla 'DB', and a 'pool' and 'selet' lists; The DB contains only disjoint sets.  One of the sets is empty.
number_of_steps <- 100

your_res_01 <- set.based.enrichment.test(steps=number_of_steps, pool=pool, 
                                      select=select, DB=DB)

model_df <- convert_model_list_to_model_df(DB)
mulea_ora_01 <- MulEA::ORA(gmt = model_df, testData = select,
                        pool = pool, adjustMethod = "PT",
                        numberOfPermutations = number_of_steps)
mulea_res_01 <- MulEA::runTest(mulea_ora_01)

# TODO : Discuss with Eszter round of FDRs bigger than 1 to 1.
your_res_01$FDR-mulea_res_01$FDR



# Test Case 02
load(file="test-parameters-01.RData") # It loas a dummy databaseinto variabla 'DB', and a 'pool' and 'selet' lists; The DB contains only disjoint sets.  One of the sets is empty.

pool = c(pool, "dummy1", "dummy2", "dummy3")
select = c(select, "s_dummy1", "s_dummy2", "s_dummy3")

number_of_steps <- 100

your_res_02 <- set.based.enrichment.test(steps=number_of_steps, pool=pool, 
                                         select=select, DB=DB)

model_df <- convert_model_list_to_model_df(DB)
mulea_ora_02 <- MulEA::ORA(gmt = model_df, testData = select,
                           pool = pool, adjustMethod = "PT",
                           numberOfPermutations = number_of_steps)
mulea_res_02 <- MulEA::runTest(mulea_ora_02)

your_res_02$FDR-mulea_res_02$FDR



# Test Case 03
load(file="test-parameters-01.RData") # It loas a dummy databaseinto variabla 'DB', and a 'pool' and 'selet' lists; The DB contains only disjoint sets.  One of the sets is empty.

pool <- c("dummy1", "dummy2", "dummy3")
select <- c("s_dummy1", "s_dummy2", "s_dummy3")

number_of_steps <- 100

your_res_03 <- set.based.enrichment.test(steps=number_of_steps, pool=pool, 
                                         select=select, DB=DB)

model_df <- convert_model_list_to_model_df(DB)
mulea_ora_03 <- MulEA::ORA(gmt = model_df, testData = select,
                           pool = pool, adjustMethod = "PT",
                           numberOfPermutations = number_of_steps)
mulea_res_03 <- MulEA::runTest(mulea_ora_03)

your_res_03$FDR-mulea_res_03$FDR


# TODO : What should catch those examples?
if(all(your_res_03$P==1)){ stop()}
if(all(your_res_03$FDR==1)){ stop()}



# Test Case 04
load(file="test-parameters-01.RData") # It loas a dummy databaseinto variabla 'DB', and a 'pool' and 'selet' lists; The DB contains only disjoint sets.  One of the sets is empty.

pool <- c(pool, "dummy1", "dummy2", "dummy3") # the DB does not contain the dummy variables
select <- c( "dummy1", "dummy2", "dummy3")

number_of_steps <- 100

your_res_04 <- set.based.enrichment.test(steps=number_of_steps, pool=pool, 
                                         select=select, DB=DB)

model_df <- convert_model_list_to_model_df(DB)
mulea_ora_04 <- MulEA::ORA(gmt = model_df, testData = select,
                           pool = pool, adjustMethod = "PT",
                           numberOfPermutations = number_of_steps)
mulea_res_04 <- MulEA::runTest(mulea_ora_04)

your_res_04$FDR-mulea_res_04$FDR



# My example from package - I am using it as a test.
# Input example data in file.
# write.csv2(x=data.frame('select'=select, 'score'=selectScores), file = './inst/extdata/selectData.csv')
# write.csv2(x=data.frame('pool'=unique(pool)), file = './inst/extdata/poolData.csv')

# Standard procedure.
# TODO : Update README.md after all.
library(MulEA)
muleaPkgDir <- find.package("MulEA")
modelDfFromFile <- MulEA::readGmtFileAsDataFrame(
  gmtFilePath = paste(muleaPkgDir,"/extdata/model.gmt", sep = ""))
selectDf <- read.csv2(file = './inst/extdata/selectData.csv')
select <- selectDf[['select']]
poolDf <- read.csv2(file = './inst/extdata/poolData.csv')
pool <- poolDf[['pool']]
number_of_steps <- 1000

# TODO : DONE : Common results across set and ranked methods.
mulea_ora_model <- MulEA::ORA(
  gmt = modelDfFromFile, testData = select, 
  pool = pool, adjustMethod = "PT",
  numberOfPermutations = number_of_steps)
mulea_ora_results <- MulEA::runTest(mulea_ora_model)

# TODO : Bug after results columns extension.
mulea_ora_reshaped_results <- MulEA::reshapeResults(
  mulea_model=mulea_ora_model, 
  mulea_model_resuts=mulea_ora_results, 
  category_stat_column_name='adjustedPValueEmpirical')


# Plotting methods tests.
mulea_ora_reshaped_def_results <- MulEA::reshapeResults(
  mulea_model=mulea_ora_model, 
  mulea_model_resuts=mulea_ora_results)

# TODO : Add trim_to_testData cutoff.
# TODO : colname is inherited from category_stat_column_name
mulea_ora_reshaped_results <- MulEA::reshapeResults(
  mulea_model=mulea_ora_model, 
  mulea_model_resuts=mulea_ora_results, 
  category_stat_column_name='adjustedPValueEmpirical')


# Plots.
# TODO : Add example with names in plot. Not worry too much. :)
MulEA::plotGraph(mulea_relaxed_resuts=mulea_ora_reshaped_results, statistics_value_cutoff = 1.00)

MulEA::plotBarplot(mulea_relaxed_resuts = mulea_ora_reshaped_results, statistics_value_cutoff=1.00)

MulEA::plotHeatmap(mulea_relaxed_resuts=mulea_ora_reshaped_results, statistics_value_cutoff=1.00)



# Plots for empirically adjusted p-values.
# Plot graph.
# TODO : You need to provide proper colname for statistics_value_colname.
colnames(mulea_relaxed_adj_pval_emp_resuts) <- c("ontologyId", "genIdInOntology", "empirivalPValue")
MulEA::plotGraph(mulea_relaxed_resuts=mulea_relaxed_adj_pval_emp_resuts, statistics_value_cutoff = 1.00,
                 statistics_value_colname = "empirivalPValue")

# Plot barplot
MulEA::plotBarplot(mulea_relaxed_resuts = mulea_relaxed_adj_pval_emp_resuts, statistics_value_cutoff=1.00)

# Plot heatmap
MulEA::plotHeatmap(mulea_relaxed_resuts=mulea_relaxed_adj_pval_emp_resuts, statistics_value_cutoff=1.00, 
                   statistics_value_colname = 'empirivalPValue')


# Subramanian plots.
# TODO : Remove the method arg as results of removal of KS. 
selectScores <- selectDf[['score']]
mulea_ranked_model <- MulEA::RankedBasedTest(
  gmt = modelDfFromFile, 
  testData = select, scores = selectScores)
mulea_sub_results <- MulEA::runTest(mulea_ranked_model)
mulea_sub_reshaped_results <- MulEA::reshapeResults(
  mulea_model = mulea_ranked_model, 
  mulea_model_resuts = mulea_sub_results, 
  mulea_model_resuts_ontology_col_name='ontologyId')

MulEA::plotGraph(mulea_relaxed_resuts=mulea_sub_reshaped_results, statistics_value_cutoff = 1.00)
MulEA::plotBarplot(mulea_relaxed_resuts = mulea_sub_reshaped_results, statistics_value_cutoff=1.00)
MulEA::plotHeatmap(mulea_relaxed_resuts=mulea_sub_reshaped_results, statistics_value_cutoff=1.00)






muleaSetBaseEnrichmentTest <- MulEA:::SetBasedEnrichmentTest(gmt = setBaseTestObject@gmt,
                                                     testData = setBaseTestObject@testData,
                                                     pool = setBaseTestObject@pool,
                                                     numberOfPermutations = setBaseTestObject@numberOfPermutations)




merge(setBaseTestObject@gmt[c('ontologyId', 'ontologyName')], muleaSetBaseEnrichmentTest, by.x = "ontologyId", 
      by.y = "DB_names", all = TRUE)
