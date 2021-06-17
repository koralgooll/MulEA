
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

# Generate artificial select vector.
no_over_repr_terms=3
no_under_repr_terms=2

set.seed(seed = 1234)
input_generated <- MulEA::generateInputData(
  input_gmt = input_gmt_filtered, sample_ratio=0.2,
  group_under_over_representation_ratio=0.8,
  number_of_over_representation_groups = no_over_repr_terms,
  number_of_under_representation_groups = no_under_repr_terms)


# Perform ORA test.
input_select <- input_generated$input_select
number_of_steps <- 1000
mulea_ora_model <- MulEA::ORA(
  gmt = input_gmt_filtered, testData = input_select, adjustMethod = "PT",
  numberOfPermutations = number_of_steps)
mulea_ora_results <- MulEA::runTest(mulea_ora_model)

warnings()

# Summary of results
mulea_ora_results_rank <- data.table::frank(mulea_ora_results, cols='adjustedPValueEmpirical')
names(mulea_ora_results_rank) <- mulea_ora_results$ontologyId 

select_over_repr <- mulea_ora_results$ontologyId %in% input_generated$go_change_repr_over
data.frame(mulea_ora_results[select_over_repr,], 
           'rank'=mulea_ora_results_rank[select_over_repr])

select_under_repr <- mulea_ora_results$ontologyId %in% input_generated$go_change_repr_under
data.frame(mulea_ora_results[select_under_repr,], 
           'rank'=mulea_ora_results_rank[select_under_repr])


# Visualization of results.
mulea_ora_reshaped_results <- MulEA::reshapeResults(
  mulea_model=mulea_ora_model, 
  mulea_model_resuts=mulea_ora_results, 
  category_stat_column_name='adjustedPValueEmpirical')

MulEA::plotBarplot(mulea_relaxed_resuts = mulea_ora_reshaped_results, 
                   statistics_value_colname = "adjustedPValueEmpirical",
                   statistics_value_cutoff=0.05)
MulEA::plotHeatmap(mulea_relaxed_resuts=mulea_ora_reshaped_results, 
                   statistics_value_colname = "adjustedPValueEmpirical",
                   statistics_value_cutoff=0.05)
MulEA::plotGraph(mulea_relaxed_resuts=mulea_ora_reshaped_results,
                 statistics_value_colname = "adjustedPValueEmpirical",
                 statistics_value_cutoff = 0.05)
