
gmtFilePath <- paste(
  find.package("MulEA"), 
  "/tests/inputs/Bands_Drosophila_melanogaster_Flybase_Marton.gmt", sep = "")
gmtFilePath <- paste(
  find.package("MulEA"), 
  "/tests/inputs/KEGG_example_c_elegans_Leila.gmt", sep = "")


input_gmt <- MulEA::readGmtFileAsDataFrame(gmtFilePath)


# Start of test of method.
library(MulEA)

number_of_steps <- 10

hist_data <- c()

no_over_repr_terms=3 
no_under_repr_terms=2

for (i in 1:1) {
  input_gmt_filtered <- filter_ontology(input_gmt = input_gmt)
  
  input_select <- generate_input_data(
    input_gmt = input_gmt_filtered, 
    number_of_over_representation_groups = no_over_repr_terms,
    number_of_under_representation_groups = no_under_repr_terms)
  hist_data <- c(hist_data, length(input_select))
  
  mulea_ora_model <- MulEA::ORA(
    gmt = input_gmt_filtered, testData = input_select, adjustMethod = "PT",
    numberOfPermutations = number_of_steps)
  mulea_ora_results <- MulEA::runTest(mulea_ora_model)
}

hist(hist_data)


# Other stats 
terms_sizes_log <- log(terms_sizes)

hist(terms_sizes, breaks=100)
hist(terms_sizes_log, breaks=100)

quantile(terms_sizes, probs = seq(0, 1, 0.1), type = 2, na.rm = FALSE)
quantile(terms_sizes_log, probs = seq(0, 1, 0.1), type = 2, na.rm = FALSE)
exp(1.386294)
exp(4.394144)

mean(terms_sizes)
mean(terms_sizes_log)
exp(3.056871)
median(terms_sizes)
median(terms_sizes_log)
sd(terms_sizes_log)


# Setup params


# Perform ORA
mulea_ora_model <- MulEA::ORA(
  gmt = input_gmt, testData = input_select, adjustMethod = "PT",
  numberOfPermutations = number_of_steps)
mulea_ora_results <- MulEA::runTest(mulea_ora_model)


mulea_ora_reshaped_results <- MulEA::reshapeResults(
  mulea_model=mulea_ora_model, 
  mulea_model_resuts=mulea_ora_results, 
  category_stat_column_name='adjustedPValueEmpirical')

warnings()

# Plot results
MulEA::plotBarplot(mulea_relaxed_resuts = mulea_ora_reshaped_results, 
                   statistics_value_colname = "adjustedPValueEmpirical",
                   statistics_value_cutoff=0.05)
MulEA::plotHeatmap(mulea_relaxed_resuts=mulea_ora_reshaped_results, 
                   statistics_value_colname = "adjustedPValueEmpirical",
                   statistics_value_cutoff=0.05)
MulEA::plotGraph(mulea_relaxed_resuts=mulea_ora_reshaped_results,
                 statistics_value_colname = "adjustedPValueEmpirical",
                 statistics_value_cutoff = 0.05)






for (i in 1:1) {
  input_gmt_filtered <- MulEA::filterOntology(input_gmt = input_gmt)
  
  input_select <- MulEA::generateInputData(
    input_gmt = input_gmt_filtered, 
    number_of_over_representation_groups = no_over_repr_terms,
    number_of_under_representation_groups = no_under_repr_terms)
  
  mulea_ora_model <- MulEA::ORA(
    gmt = input_gmt_filtered, testData = input_select, adjustMethod = "PT",
    numberOfPermutations = number_of_steps)
  mulea_ora_results <- MulEA::runTest(mulea_ora_model)
}




