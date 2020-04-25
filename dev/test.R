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
muleaPkgDir <- find.package("MulEA")
modelDfFromFile <- MulEA::readGmtFileAsDataFrame(gmtFilePath = paste(muleaPkgDir,"/extdata/model.gmt", sep = ""))
DB <- modelDfFromFile[, 'listOfValues']
names(DB) <- modelDfFromFile$ontologyId
select <- c("FBgn0004407", "FBgn0010438", "FBgn0003742", "FBgn0029709", 
                       "FBgn0030341", "FBgn0037044", "FBgn0002887", "FBgn0028434", 
                       "FBgn0030170", 
                       "FBgn0263831", "FBgn0261618", "FBgn0038704", "FBgn0000579")
pool <- unique(c(c("FBgn0033690", "FBgn0261618", "FBgn0004407", "FBgn0010438", "FBgn0032154", "FBgn0039930", "FBgn0040268", "FBgn0013674",
                             "FBgn0037008", "FBgn0003116", "FBgn0037743", "FBgn0035401", "FBgn0037044", "FBgn0051005", "FBgn0026737", "FBgn0026751",
                             "FBgn0038704", "FBgn0002887", "FBgn0028434", "FBgn0030170", "FBgn0263831", "FBgn0000579"),
                           c("FBgn0066666", "FBgn0000000", "FBgn0099999", "FBgn0011111", "FBgn0022222", "FBgn0777777", "FBgn0333333", "FBgn0003742",
                             "FBgn0029709", "FBgn0030341")))

number_of_steps <- 100

your_res_M <- set.based.enrichment.test(steps=number_of_steps, pool=pool, 
                                         select=select, DB=DB)

mulea_ora_M <- MulEA::ORA(gmt = modelDfFromFile, testData = select,
                           pool = pool, adjustMethod = "PT",
                           numberOfPermutations = number_of_steps)
mulea_res_M <- MulEA::runTest(mulea_ora_M)

your_res_M$FDR-mulea_res_M$FDR


