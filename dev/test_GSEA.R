
library(MulEA)

# Read inputs
input_gmt <- MulEA::readGmtFileAsDataFrame('D:/projects/science/GSEA_4.1.0/inputs/KEGG_Paths2geneSymbols.gmt')
input_rank <- read.table(file='D:/projects/science/GSEA_4.1.0/inputs/lung_FoldChange.rnk', sep='\t')
input_select <- input_rank[['V1']]
input_select_scores <- input_rank[['V2']]

# Setup params
number_of_steps <- 10

# Perform GSEA
# IMPORTANT : Scores are very different based on scoreType param.
mulea_ranked_model <- MulEA::RankedBasedTest(
  gmt = input_gmt, testData = input_select, scores = input_select_scores, scoreType = "pos")
mulea_sub_results <- MulEA::runTest(mulea_ranked_model)
mulea_sub_reshaped_results <- MulEA::reshapeResults(
  mulea_model = mulea_ranked_model, 
  mulea_model_resuts = mulea_sub_results, 
  mulea_model_resuts_ontology_col_name='ontologyId')

warnings()

# Plot results
MulEA::plotBarplot(mulea_relaxed_resuts = mulea_sub_reshaped_results, 
                   statistics_value_cutoff=0.05)
MulEA::plotBarplot(mulea_relaxed_resuts = mulea_sub_reshaped_results, 
                   statistics_value_cutoff=0.20)
MulEA::plotHeatmap(mulea_relaxed_resuts=mulea_sub_reshaped_results, 
                   statistics_value_cutoff=0.05)
MulEA::plotHeatmap(mulea_relaxed_resuts=mulea_sub_reshaped_results, 
                   statistics_value_cutoff=0.20)
MulEA::plotGraph(mulea_relaxed_resuts=mulea_sub_reshaped_results, 
                 statistics_value_cutoff = 0.05)
MulEA::plotGraph(mulea_relaxed_resuts=mulea_sub_reshaped_results, 
                 statistics_value_cutoff = 0.20)
# DONOTRUN : Takes time, but shows interesting warnings. 
MulEA::plotGraph(mulea_relaxed_resuts=mulea_sub_reshaped_results, 
                 statistics_value_cutoff = 1.00)


# Start GSEA (Broad Institute) and read results.
edbFilePath <- 'D:/projects/science/GSEA_4.1.0/results/my_analysis.GseaPreranked.1613825545589/edb/results.edb'
MulEA::readEdbFileAsDataFrame(edbFilePath = edbFilePath)








