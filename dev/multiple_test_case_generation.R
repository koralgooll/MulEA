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
  input_gmt_filtered=input_gmt_filtered, 
  over_repr_ratio = 0.6,
  number_of_tests = 2)
sim_mult_tests_res_sum <- MulEA:::getMultipleTestsSummaryAcrossCutOff(
  tests_res=sim_mult_tests_res)

sim_mult_tests_res <- MulEA:::simulateMultipleTestsWithRatioParam(
  input_gmt_filtered=input_gmt_filtered, 
  over_repr_ratio = 0.6,
  number_of_tests = 100)

# write_rds(sim_mult_tests_res, "dev\\sim_mult_tests_small_res.rds")
sim_mult_tests_res <- readr::read_rds("dev\\sim_mult_tests_small_res.rds")

sim_mult_tests_res_sum_small <- MulEA:::getMultipleTestsSummaryAcrossCutOff(
  tests_res=sim_mult_tests_res, cut_off_range = seq(0, 1, 0.01))

# write_rds(sim_mult_tests_res, "dev\\sim_mult_tests_big_res.rds")
sim_mult_tests_res <- read_rds("dev\\sim_mult_tests_big_res.rds")


sim_mult_tests_res_sum_big <- MulEA:::getMultipleTestsSummaryAcrossCutOff(
  tests_res=sim_mult_tests_res, cut_off_range = seq(0, 1, 0.01))

sim_mult_tests_res_sum <- sim_mult_tests_res_sum_big
sim_mult_tests_res_sum <- sim_mult_tests_res_sum_small

comp_mult_tests <- sim_mult_tests_res_sum
mult_tests_03_04_sum <- sim_mult_tests_res_sum
# GOTO : Noise comparison tests (line ~80, plots)
