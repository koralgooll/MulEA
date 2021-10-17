
library(MulEA)

## TODO : Question about base files.
# Read and filter inputs.
gmtFilePath <- paste(find.package("MulEA"), 
                     "/tests/inputs/KEGG_example_c_elegans_Leila.gmt", sep = "")
input_gmt <- MulEA::readGmtFileAsDataFrame(gmtFilePath)
input_gmt_filtered <- MulEA::filterOntology(input_gmt = input_gmt)



filteredGmtFilePath <- paste(find.package("MulEA"), 
                             "/tests/outputs/KEGG_filtered.gmt", sep = "")
MulEA::saveDataFrameAsGmtFile(modelDF = input_gmt_filtered, gmtFilePath = filteredGmtFilePath)


# Read GO.
input_gmt_filtered <- MulEA::readGmtFileAsDataFrame(gmtFilePath = filteredGmtFilePath)



compare_methods <- function(input_gmt_filtered, 
                            no_over_repr_terms=1, no_under_repr_terms=1,
                            sample_ratio=0.2, group_under_over_representation_ratio=0.5,
                            seed=NULL, turn_on_log=FALSE, number_of_perm_in_adj=10000) {
  if (!is.null(seed)) {
    set.seed(seed = seed)
  }
  
  input_generated <- MulEA::generateInputData(
    input_gmt = input_gmt_filtered, sample_ratio=sample_ratio,
    group_under_over_representation_ratio=group_under_over_representation_ratio,
    number_of_over_representation_groups = no_over_repr_terms,
    number_of_under_representation_groups = no_under_repr_terms, 
    turn_on_log = turn_on_log)
  
  # Perform ORA test.
  if (turn_on_log) {
    print(input_generated)
  }
  input_select <- input_generated$input_select
  
  number_of_steps <- number_of_perm_in_adj
  mulea_ora_model <- MulEA::ORA(
    gmt = input_gmt_filtered, testData = input_select, adjustMethod = "PT",
    numberOfPermutations = number_of_steps)
  mulea_ora_results <- MulEA::runTest(mulea_ora_model)
  
  if (turn_on_log) {
    debug_mulea_ora_results <<- mulea_ora_results
    print(mulea_ora_results)
  }
  # warnings()
  
  
  # Summary of results
  
  # By values
  mulea_ora_results_value <- mulea_ora_results[
    mulea_ora_results$ontologyId %in% input_generated$go_change_repr_over, 
    c('adjustedPValue', 'adjustedPValueEmpirical')]
  
  mulea_ora_results_value_better <- sum(mulea_ora_results_value$adjustedPValue > mulea_ora_results_value$adjustedPValueEmpirical)
  mulea_ora_results_value_worse <- sum(mulea_ora_results_value$adjustedPValue < mulea_ora_results_value$adjustedPValueEmpirical)
  
  if (turn_on_log) {
    print('---By value---')
    to_print <- paste0('Better -> ', mulea_ora_results_value_better, ', ', mulea_ora_results_value_worse, ' <- worse.')
    print(to_print)
  }
  
  
  # By rank
  mulea_ora_perm_results_rank <- data.table::frank(mulea_ora_results, cols='adjustedPValueEmpirical')
  names(mulea_ora_perm_results_rank) <- mulea_ora_results$ontologyId 
  
  mulea_ora_bh_results_rank <- data.table::frank(mulea_ora_results, cols='adjustedPValue')
  names(mulea_ora_bh_results_rank) <- mulea_ora_results$ontologyId 
  
  known_over_repr_perm_ranks <- mulea_ora_perm_results_rank[input_generated$go_change_repr_over]
  known_over_repr_bh_ranks <- mulea_ora_bh_results_rank[input_generated$go_change_repr_over]
  
  known_over_repr_perm_ranks_mean <- mean(known_over_repr_perm_ranks)
  known_over_repr_bh_ranks_mean <- mean(known_over_repr_bh_ranks)
  
  known_over_repr_perm_ranks_median <- median(known_over_repr_perm_ranks)
  known_over_repr_bh_ranks_median <- median(known_over_repr_bh_ranks)
  
  if (turn_on_log) {
    print('---By rank---')
    to_print <- paste0('Perm rank mean -> ', known_over_repr_perm_ranks_mean, ', ', 
                       known_over_repr_bh_ranks_mean, ' <- bh rank mean')
    print(to_print)
    to_print <- paste0('Perm rank median -> ', known_over_repr_perm_ranks_median, ', ', 
                       known_over_repr_bh_ranks_median, ' <- bh rank median')
    print(to_print)
  }
  
  return(c(mulea_ora_results_value_better >= mulea_ora_results_value_worse, 
           known_over_repr_perm_ranks_mean <= known_over_repr_bh_ranks_mean, 
           known_over_repr_perm_ranks_median <= known_over_repr_bh_ranks_median))
  
}


# Find optimal usage.

no_over_repr_terms <- 1:20
no_under_repr_terms <- 0:10
sample_ratio <- seq(0, 0.3, by=0.1)
group_under_over_representation_ratio <- seq(0.1, 0.9, by=0.2)
value_score <- 0
rank_mean_score <- 0
rank_median_score <- 0

comparison_space <- expand.grid(
  no_over_repr_terms=no_over_repr_terms, no_under_repr_terms=no_under_repr_terms,
  sample_ratio=sample_ratio, 
  group_under_over_representation_ratio=group_under_over_representation_ratio,
  value_score=value_score, rank_mean_score=rank_mean_score, rank_median_score=rank_median_score)


# for (row_i in 1:nrow(comparison_space)) {
for (row_i in 1:20) {
  comparison_res <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(comparison_res) <- c('by_value', 'by_rank_mean', 'by_rank_median')
  
  current_comp_space <- comparison_space[row_i,]
  
  no_of_iters <- 10
  for (variable in 1:no_of_iters) {
    comparison_res[nrow(comparison_res) + 1,] <- compare_methods(
      input_gmt_filtered = input_gmt_filtered, 
      no_over_repr_terms=current_comp_space$no_over_repr_terms, 
      no_under_repr_terms=current_comp_space$no_under_repr_terms,
      sample_ratio=current_comp_space$sample_ratio, 
      group_under_over_representation_ratio=current_comp_space$group_under_over_representation_ratio, 
      turn_on_log = TRUE)
  }
  
  comparison_space[row_i,]$value_score <- sum(comparison_res$by_value)/no_of_iters
  comparison_space[row_i,]$rank_mean_score <- sum(comparison_res$by_rank_mean)/no_of_iters
  comparison_space[row_i,]$rank_median_score <- sum(comparison_res$by_rank_median)/no_of_iters
}



compare_methods(input_gmt_filtered = input_gmt_filtered, turn_on_log = TRUE)



comparison_res <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(comparison_res) <- c('by_value', 'by_rank_mean', 'by_rank_median')

no_of_iters <- 1000
for (variable in 1:no_of_iters) {
  comparison_res[nrow(comparison_res) + 1,] = compare_methods(input_gmt_filtered = input_gmt_filtered)
}

comparison_res

sum(comparison_res$by_value)/no_of_iters
sum(comparison_res$by_rank_mean)/no_of_iters
sum(comparison_res$by_rank_median)/no_of_iters



input_generated_serilize <- input_generated





input_generated <- MulEA::generateInputData(
  input_gmt = input_gmt_filtered, sample_ratio=0.1,
  group_under_over_representation_ratio=0.5,
  number_of_over_representation_groups = 3,
  number_of_under_representation_groups = 1, 
  turn_on_log = FALSE)

# Here is sample input.
print(input_generated)
input_select <- input_generated$input_select

# This arg is responsible for number of permutations.
number_of_steps <- 1000
mulea_ora_model <- MulEA::ORA(
  gmt = input_gmt_filtered, testData = input_select, adjustMethod = "PT",
  numberOfPermutations = number_of_steps)
mulea_ora_results <- MulEA::runTest(mulea_ora_model)

# Here will be zeros. If not re run the procedure.
mulea_ora_results



selectFilePath <- paste(find.package("MulEA"), 
                        "/tests/outputs/input_select.rds", sep = "")
saveRDS(object = input_select, file = selectFilePath)
input_select_from_file <- readRDS(file = selectFilePath)

selectCsvFilePath <- paste(find.package("MulEA"), 
                           "/tests/outputs/input_select.csv", sep = "")
write.csv(input_select, file = selectCsvFilePath)




