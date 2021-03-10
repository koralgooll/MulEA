
library(MulEA)

# Read inputs
input_gmt <- MulEA::readGmtFileAsDataFrame('D:/projects/science/GSEA_4.1.0/inputs/KEGG_Paths2geneSymbols.gmt')
input_rank <- read.table(file='D:/projects/science/GSEA_4.1.0/inputs/lung_FoldChange.rnk', sep='\t')
input_select <- input_rank[['V1']]
input_select_scores <- input_rank[['V2']]

# Setup params
number_of_steps <- 10


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
