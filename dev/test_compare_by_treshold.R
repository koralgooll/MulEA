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

# Set seed.
set.seed(seed = 1234)

# Construct search space.
construct_comparison_space <- function(nubmer_of_over_repr_terms_range=5:10, 
                                       nubmer_of_under_repr_terms_range=5:10,
                                       sample_ratio_range=seq(0.1, 0.3, by=0.1),
                                       over_repr_ratio_range=seq(0.1, 0.5, by=0.2), 
                                       under_repr_ratio_range=seq(0.1, 0.2, by=0.1)) {
  comparison_space <- expand.grid(
    no_over_repr_terms=nubmer_of_over_repr_terms_range, 
    no_under_repr_terms=nubmer_of_under_repr_terms_range,
    sample_ratio=sample_ratio_range, 
    over_repr_ratio=over_repr_ratio_range,
    under_repr_ratio=under_repr_ratio_range,
    score=0)
  
  # TODO : Validation of space
  comparison_space <- comparison_space[comparison_space$sample_ratio >= comparison_space$under_repr_ratio,]
  comparison_space <- comparison_space[comparison_space$sample_ratio+comparison_space$over_repr_ratio <= 1,]
  
  return(comparison_space)
}

comparison_space <- construct_comparison_space()

# Call method on search space.
compare_by_treshold <- function(comparison_space, input_gmt_filtered) {
  for (i in 1:length(comparison_space$no_over_repr_terms)) {
    comparison_case <- comparison_space[i,]
    print(i)
    
    input_generated <- MulEA::generateInputData(
      input_gmt = input_gmt_filtered, 
      sample_ratio=comparison_case$sample_ratio,
      over_repr_ratio = comparison_case$over_repr_ratio, 
      under_repr_ratio = comparison_case$under_repr_ratio,
      number_of_over_representation_groups = comparison_case$no_over_repr_terms,
      number_of_under_representation_groups = comparison_case$no_under_repr_terms)
    
    # Perform ORA test.
    input_select <- input_generated$input_select
    number_of_steps <- 10000
    
    # TODO : benchmark time.
    mulea_ora_model <- MulEA::ORA(
      gmt = input_gmt_filtered, testData = input_select, adjustMethod = "PT",
      numberOfPermutations = number_of_steps, nthreads = 16)
    
    mulea_ora_results <- NULL
    mulea_ora_results <- MulEA::runTest(mulea_ora_model)
  }
}

compare_by_treshold(comparison_space = comparison_space, 
                    input_gmt_filtered = input_gmt_filtered)



# DELETE - for debug only
input_gmt=input_gmt_filtered
sample_ratio=0.3
over_repr_ratio=0.5 
under_repr_ratio=0.2
number_of_over_representation_groups=10 
number_of_under_representation_groups=10 
turn_on_log=FALSE







