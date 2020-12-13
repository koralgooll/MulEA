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

# TODO : Improve description and README and etc.
# Standard procedure.
library(MulEA)
muleaPkgDir <- find.package("MulEA")
modelDfFromFile <- MulEA::readGmtFileAsDataFrame(gmtFilePath = paste(muleaPkgDir,"/extdata/model.gmt", sep = ""))

selectDf <- read.csv2(file = './inst/extdata/selectData.csv')
select <- selectDf[['select']]
poolDf <- read.csv2(file = './inst/extdata/poolData.csv')
pool <- poolDf[['pool']]
number_of_steps <- 1000

mulea_ora_M <- MulEA::ORA(gmt = modelDfFromFile, testData = select,
                           pool = pool, adjustMethod = "PT",
                           numberOfPermutations = number_of_steps)
mulea_res_M <- MulEA::runTest(mulea_ora_M)


# Plotting methods tests.
mulea_relaxed_pval_resuts <- MulEA::createDetailedResults(mulea_model=mulea_ora_M, 
                                                    mulea_model_resuts=mulea_res_M)
# TODO: DONE : plot base on p.val and p.adj possibility.
mulea_relaxed_adj_pval_emp_resuts <- MulEA::createDetailedResults(mulea_model=mulea_ora_M, 
                                                     mulea_model_resuts=mulea_res_M,
                                                     category_stat_column_name='adjustedPValueEmpirical')


# Plots for p-values.
# TODO : REJECTED : Let's the user choose between names or ids on plots.
# Plot graph.
MulEA::plotGraph(mulea_relaxed_resuts=mulea_relaxed_pval_resuts, statistics_value_cutoff = 1.00)

# Plot barplot
MulEA::plotBarplot(mulea_relaxed_resuts = mulea_relaxed_pval_resuts, statistics_value_cutoff=1.00)

# Plot heatmap
MulEA::plotHeatmap(mulea_relaxed_resuts=mulea_relaxed_pval_resuts, statistics_value_cutoff=1.00)


# Plots for empirically adjusted p-values.
# Plot graph.
MulEA::plotGraph(mulea_relaxed_resuts=mulea_relaxed_adj_pval_emp_resuts, statistics_value_cutoff = 1.00)

# Plot barplot
MulEA::plotBarplot(mulea_relaxed_resuts = mulea_relaxed_adj_pval_emp_resuts, statistics_value_cutoff=1.00)

# Plot heatmap
MulEA::plotHeatmap(mulea_relaxed_resuts=mulea_relaxed_adj_pval_emp_resuts, statistics_value_cutoff=1.00)


# Subramanian plots.
# TODO : DONE : Extend or check if BH correction is done on p values.
selectScores <- selectDf[['score']]
rankedBasedTestSubramanian <- MulEA::RankedBasedTest(
  method = "Subramanian", gmt = modelDfFromFile, 
  testData = select, scores = selectScores)

# TODO : listOfValues removed, as it brakes consistency.
mulea_res_sub <- MulEA::runTest(rankedBasedTestSubramanian)

# TODO : create Detailed Results - think about function name. 
mulea_relaxed_resuts_sub <- MulEA::createDetailedResults(
  mulea_model = rankedBasedTestSubramanian, 
  mulea_model_resuts = mulea_res_sub, 
  mulea_model_resuts_ontology_col_name='ontologyId')

# Plot graph.
MulEA::plotGraph(mulea_relaxed_resuts=mulea_relaxed_resuts_sub, statistics_value_cutoff = 0.35)
MulEA::plotGraph(mulea_relaxed_resuts=mulea_relaxed_resuts_sub, statistics_value_cutoff = 1.00)

# Plot barplot
MulEA::plotBarplot(mulea_relaxed_resuts = mulea_relaxed_resuts_sub, statistics_value_cutoff = 0.70)
MulEA::plotBarplot(mulea_relaxed_resuts = mulea_relaxed_resuts_sub, statistics_value_cutoff = 1.00)

# Plot heatmap
MulEA::plotHeatmap(mulea_relaxed_resuts=mulea_relaxed_resuts_sub, statistics_value_cutoff = 1.00)


# Kolmogorov-Smirnov plots.
rankedBasedTestKS <- MulEA::RankedBasedTest(
  method = "KS", gmt = modelDfFromFile, 
  testData = select, scores = selectScores)
mulea_res_ks <- MulEA::runTest(rankedBasedTestKS)

# TODO : createDetailedResults - think about function name. 
mulea_relaxed_resuts_ks <- MulEA::createDetailedResults(
  mulea_model=rankedBasedTestKS, 
  mulea_model_resuts=mulea_res_ks,
  mulea_model_resuts_ontology_col_name='ontologyId')

# Plot graph.
MulEA::plotGraph(mulea_relaxed_resuts=mulea_relaxed_resuts_ks, statistics_value_cutoff = 1.00)
