# IMPORTANT FILE!!! Generate statistics to plots!
library(MulEA)
library(readr)
library(tidyverse)

# Read and filter inputs (small).
gmtFilePath <- paste(find.package("MulEA"), 
                     "/tests/inputs/KEGG_example_c_elegans_Leila.gmt", sep = "")
input_gmt <- MulEA::readGmtFileAsDataFrame(gmtFilePath)
input_gmt_filtered <- MulEA::filterOntology(input_gmt = input_gmt)
filteredGmtFilePath <- paste(find.package("MulEA"), 
                             "/tests/outputs/KEGG_filtered.gmt", sep = "")
MulEA::saveDataFrameAsGmtFile(modelDF = input_gmt_filtered, gmtFilePath = filteredGmtFilePath)
input_gmt_filtered <- MulEA::readGmtFileAsDataFrame(gmtFilePath = filteredGmtFilePath)

# Read and filter inputs (big).
gmtFilePath <- paste(find.package("MulEA"), 
                     "/tests/inputs/Pfam_Uniprot_Marton_Homo_sapiens.gmt", sep = "")
input_gmt <- MulEA::readGmtFileAsDataFrame(gmtFilePath)
input_gmt_filtered <- MulEA::filterOntology(input_gmt = input_gmt)
filteredGmtFilePath <- paste(find.package("MulEA"), 
                             "/tests/outputs/Pfam_Uniprot_Marton_Homo_sapiens.gmt", sep = "")
MulEA::saveDataFrameAsGmtFile(modelDF = input_gmt_filtered, gmtFilePath = filteredGmtFilePath)
input_gmt_filtered <- MulEA::readGmtFileAsDataFrame(gmtFilePath = filteredGmtFilePath)


# DEBUG : Global input data.
noise_ratio = 0.5
over_repr_ratio = 0.5
under_repr_ratio = 0.05
number_of_over_representation_groups = 5
number_of_under_representation_groups = 0
number_of_samples = 1
number_of_tests = 100
number_of_steps = 5000


# GREAT IDEA : ML State Management!!!
# Always same params to functions, parameters cascading is unhealthy. 
# To understand look at params of simulateMultipleTestsWithRatioParam, 
# then go inside into MulEA:::simulateMultipleTests, 
# then go inside into MulEA::ORA, MulEA:::decorateGmtByUnderOvenAndNoise, 
# MulEA:::generateInputSamples. Real waterfall of params. :)
sim_mult_tests_res <- MulEA:::simulateMultipleTestsWithRatioParam(
  input_gmt_filtered = input_gmt_filtered, 
  over_repr_ratio = 0.85,
  number_of_tests = 3, nthreads = 16)
sim_mult_tests_res <- MulEA:::simulateMultipleTestsWithRatioParam(
  input_gmt_filtered = input_gmt_filtered, 
  over_repr_ratio = 0.85,
  number_of_tests = 10, nthreads = 16)
print(object.size(sim_mult_tests_res))
sim_mult_tests_res_a <- MulEA:::simulateMultipleTestsWithRatioParam(
  input_gmt_filtered = input_gmt_filtered, 
  over_repr_ratio = 0.75,
  number_of_tests = 30, nthreads = 16)
sim_mult_tests_res_b <- MulEA:::simulateMultipleTestsWithRatioParam(
  input_gmt_filtered = input_gmt_filtered, 
  over_repr_ratio = 0.85,
  number_of_tests = 30, nthreads = 16)
sim_mult_tests_res_c <- MulEA:::simulateMultipleTestsWithRatioParam(
  input_gmt_filtered = input_gmt_filtered, 
  over_repr_ratio = 0.95,
  number_of_tests = 30, nthreads = 16)
print(object.size(sim_mult_tests_res))

# readr::write_rds(sim_mult_tests_res, file = "dev\\new_tests_res\\sim_mult_tests_res_small_085_10.rds")
sim_mult_tests_res <- readr::read_rds("dev\\new_tests_res\\sim_mult_tests_res_small_085_3.rds")
sim_mult_tests_res <- readr::read_rds("dev\\new_tests_res\\sim_mult_tests_res_big_085_10.rds")
sim_mult_tests_res <- readr::read_rds("dev\\new_tests_res\\sim_mult_tests_res_small_085_10.rds")
sim_mult_tests_res_sum <- MulEA:::getMultipleTestsSummaryAcrossCutOff(
  tests_res=sim_mult_tests_res)
print(object.size(sim_mult_tests_res_sum))
sim_mult_tests_res_to_roc <- MulEA::getSummaryToRoc(tests_res = sim_mult_tests_res)



# WORK to PAPER.
# Was over_repr_ratio = 0.6, 
# probably because of the best performance of method in this space?
sim_mult_tests_res <- MulEA:::simulateMultipleTestsWithRatioParam(
  input_gmt_filtered=input_gmt_filtered, 
  over_repr_ratio = 0.95,
  number_of_tests = 100, 
  nthreads = 16)


# WORK to PAPER.
sim_mult_tests_res <- read_rds("dev\\new_tests_res\\sim_mult_tests_res_big_075_30.rds")
sim_mult_tests_res_sum <- MulEA:::getMultipleTestsSummaryAcrossCutOff(
  tests_res=sim_mult_tests_res)


comp_mult_tests <- sim_mult_tests_res_sum
mult_tests_03_04_sum <- sim_mult_tests_res_sum
# GOTO : Noise comparison tests (line ~80, plots)
