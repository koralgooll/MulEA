
gmtFilePath <- paste(
  find.package("MulEA"), 
  "/tests/inputs/Bands_Drosophila_melanogaster_Flybase_Marton.gmt", sep = "")
gmtFilePath <- paste(
  find.package("MulEA"), 
  "/tests/inputs/KEGG_example_c_elegans_Leila.gmt", sep = "")


input_gmt <- MulEA::readGmtFileAsDataFrame(gmtFilePath)

generate_input_data <- function(input_gmt, sample_ratio=0.5,
                                group_under_over_representation_ratio=0.9,
                                number_of_over_representation_groups=1, 
                                number_of_under_representation_groups=1) {
  input_select <- plyr::llply(.data = input_gmt$listOfValues, .fun = function(term) {
    size <- length(term)
    size_of_sample <- floor(size * sample_ratio)
    input_select_indicator <- sample(1:size, size_of_sample, replace=FALSE)
    input_select_term <- term[input_select_indicator]
  })
  
  input_select <- unique(unlist(input_select))
  input_select
  
  go_size <- length(input_gmt$listOfValues)
  size_of_over_under_repr <- number_of_over_representation_groups + number_of_under_representation_groups
  go_change_repr <- sample(1:go_size, size_of_over_under_repr, replace=FALSE)
  
  go_change_repr_over <- go_change_repr[1:number_of_over_representation_groups]
  go_change_repr_under <- go_change_repr[-1:-number_of_over_representation_groups]
  
  print('***')
  
  print("Over representation terms:")
  for (term_id in go_change_repr_over) {
    term <- input_gmt$listOfValues[[term_id]]
    size <- length(term)
    size_of_sample <- floor(size * group_under_over_representation_ratio)
    input_select_indicator <- sample(1:size, size_of_sample, replace=FALSE)
    input_select_term <- term[input_select_indicator]
    input_select <- c(input_select, input_select_term)
    print(paste(term_id, ' [', size, '] ', ' -> ', ' [', length(input_select_term), ']', sep = ''))
  }
  
  print("Under representation terms:")
  for (term_id in go_change_repr_under) {
    term <- input_gmt$listOfValues[[term_id]]
    size <- length(term)
    size_of_sample <- floor(size * group_under_over_representation_ratio)
    input_select_indicator <- sample(1:size, size_of_sample, replace=FALSE)
    input_select_term <- term[input_select_indicator]
    input_select <- setdiff(input_select, input_select_term)
    print(paste(term_id, ' [', size, '] ', ' -> ', ' [', length(input_select_term), ']', sep = ''))
  }

  input_select
}


library(MulEA)

number_of_steps <- 10

hist_data <- c()
for (i in 1:100) {
  input_select <- generate_input_data(input_gmt = input_gmt)
  hist_data <- c(hist_data, length(input_select))
  
  mulea_ora_model <- MulEA::ORA(
    gmt = input_gmt, testData = input_select, adjustMethod = "PT",
    numberOfPermutations = number_of_steps)
  mulea_ora_results <- MulEA::runTest(mulea_ora_model)
}

hist(hist_data)





terms_sizes <- plyr::llply(.data = input_gmt$listOfValues, .fun = function(term) {
  length(term)
})

terms_sizes <- unlist(terms_sizes)

hist(terms_sizes, breaks=100)
quantile(terms_sizes)



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
                   statistics_value_cutoff=0.00005)
MulEA::plotHeatmap(mulea_relaxed_resuts=mulea_ora_reshaped_results, 
                   statistics_value_colname = "adjustedPValueEmpirical",
                   statistics_value_cutoff=0.00005)
MulEA::plotGraph(mulea_relaxed_resuts=mulea_ora_reshaped_results,
                 statistics_value_colname = "adjustedPValueEmpirical",
                 statistics_value_cutoff = 0.00005)









