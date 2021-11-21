library(MulEA)

# Read and filter inputs.
gmtFilePath <- paste(find.package("MulEA"), 
                     "/tests/inputs/KEGG_example_c_elegans_Leila.gmt", sep = "")
input_gmt <- MulEA::readGmtFileAsDataFrame(gmtFilePath)
input_gmt_filtered <- MulEA::filterOntology(input_gmt = input_gmt)
filteredGmtFilePath <- paste(find.package("MulEA"), 
                             "/tests/outputs/KEGG_filtered.gmt", sep = "")
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
simulateMultipleTestsWithRatioParam <- function(
  input_gmt_filtered,
  noise_ratio_range=seq(0.1, 0.5, 0.1), 
  number_of_tests = 100, 
  over_repr_ratio = 0.5, number_of_over_representation_groups = 20, 
  number_of_steps = 5000, nthreads = 16) {
  tictoc::tic()
  sim_mult_tests <- list()
  for (noise_ratio in noise_ratio_range) {
    print("noise_ratio")
    print(noise_ratio)
    sim_mult_tests <- c(
      sim_mult_tests, 
      MulEA:::simulateMultipleTests(
        input_gmt_filtered = input_gmt_filtered, number_of_tests = number_of_tests, 
        noise_ratio = noise_ratio, over_repr_ratio = over_repr_ratio,
        number_of_over_representation_groups = number_of_over_representation_groups, 
        number_of_under_representation_groups = 0, 
        number_of_steps = number_of_steps, nthreads = nthreads)
    )
  }
  print('MulEA : ratio search calculation time:')
  tictoc::toc()
  return(sim_mult_tests)
}

sim_mult_tests_res <- simulateMultipleTestsWithRatioParam(
  input_gmt_filtered=input_gmt_filtered, 
  number_of_tests = 100)


write_rds(sim_mult_tests_res, "sim_mult_tests_res.rds")
