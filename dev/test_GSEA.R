
library(MulEA)

# Read inputs
# TODO : Change to relative package path when Eszter accepts publicity of data.
input_gmt <- MulEA::readGmtFileAsDataFrame('D:/projects/science/MulEA
                                           /inst/tests/inputs/KEGG_Paths2geneSymbols.gmt')

# Lung
input_rank_lung_fc <- read.table(file='D:/projects/science/MulEA/inst/tests/inputs/lung_FoldChange.rnk', 
                         sep='\t')
input_select_lung_fc <- input_rank_lung_fc[['V1']]
input_select_scores_lung_fc <- input_rank_lung_fc[['V2']]

input_rank_lung_netrank <- read.table(file='D:/projects/science/MulEA/inst/tests/inputs/lung_Netrank.rnk', 
                                 sep='\t')
input_select_lung_netrank <- input_rank_lung_netrank[['V1']]
input_select_scores_lung_netrank <- input_rank_lung_netrank[['V2']]

# Cervix
input_rank_cervix_fc <- read.table(file='D:/projects/science/MulEA/inst/tests/inputs/cervix_FoldChange.rnk', 
                                 sep='\t')
input_select_cervix_fc <- input_rank_cervix_fc[['V1']]
input_select_scores_cervix_fc <- input_rank_cervix_fc[['V2']]

input_rank_cervix_netrank <- read.table(file='D:/projects/science/MulEA/inst/tests/inputs/cervix_Netrank.rnk', 
                                    sep='\t')
input_select_cervix_netrank <- input_rank_cervix_netrank[['V1']]
input_select_scores_cervix_netrank <- input_rank_cervix_netrank[['V2']]

# Kidney
input_rank_kidney_fc <- read.table(file='D:/projects/science/MulEA/inst/tests/inputs/kidney_FoldChange.rnk', 
                                   sep='\t')
input_select_kidney_fc <- input_rank_kidney_fc[['V1']]
input_select_scores_kidney_fc <- input_rank_kidney_fc[['V2']]

input_rank_kidney_netrank <- read.table(file='D:/projects/science/MulEA/inst/tests/inputs/kidney_Netrank.rnk', 
                                        sep='\t')
input_select_kidney_netrank <- input_rank_kidney_netrank[['V1']]
input_select_scores_kidney_netrank <- input_rank_kidney_netrank[['V2']]


# Setup params - NOT used now. :P
number_of_steps <- 10

# Perform GSEA Lung
mulea_ranked_model_lung_fc <- MulEA::RankedBasedTest(
  gmt = input_gmt, testData = input_select_lung_fc, scores = input_select_scores_lung_fc, scoreType = "pos")
mulea_sub_results_lung_fc <- MulEA::runTest(mulea_ranked_model_lung_fc)
mulea_sub_reshaped_results_lung_fc <- MulEA::reshapeResults(
  mulea_model = mulea_ranked_model_lung_fc, 
  mulea_model_resuts = mulea_sub_results_lung_fc, 
  mulea_model_resuts_ontology_col_name='ontologyId')

warnings()

# Plot results
MulEA::plotBarplot(mulea_relaxed_resuts = mulea_sub_reshaped_results_lung_fc, 
                   statistics_value_cutoff=0.05)
MulEA::plotHeatmap(mulea_relaxed_resuts=mulea_sub_reshaped_results_lung_fc, 
                   statistics_value_cutoff=0.05)
MulEA::plotGraph(mulea_relaxed_resuts=mulea_sub_reshaped_results_lung_fc, 
                 statistics_value_cutoff = 0.05)

# Perform GSEA Lung netrank
mulea_ranked_model_lung_netrank <- MulEA::RankedBasedTest(
  gmt = input_gmt, testData = input_select_lung_netrank, scores = input_select_scores_lung_netrank, scoreType = "pos")
mulea_sub_results_lung_netrank <- MulEA::runTest(mulea_ranked_model_lung_netrank)
mulea_sub_reshaped_results_lung_netrank <- MulEA::reshapeResults(
  mulea_model = mulea_ranked_model_lung_netrank, 
  mulea_model_resuts = mulea_sub_results_lung_netrank, 
  mulea_model_resuts_ontology_col_name='ontologyId')

warnings()

# Plot results
MulEA::plotBarplot(mulea_relaxed_resuts = mulea_sub_reshaped_results_lung_netrank, 
                   statistics_value_cutoff=0.00000001)
MulEA::plotHeatmap(mulea_relaxed_resuts=mulea_sub_reshaped_results_lung_netrank, 
                   statistics_value_cutoff=0.00000001)
MulEA::plotGraph(mulea_relaxed_resuts=mulea_sub_reshaped_results_lung_netrank, 
                 statistics_value_cutoff = 0.00000001)


# Perform GSEA Cervix
mulea_ranked_model_cervix_fc <- MulEA::RankedBasedTest(
  gmt = input_gmt, testData = input_select_cervix_fc, scores = input_select_scores_cervix_fc, scoreType = "pos")
mulea_sub_results_cervix_fc <- MulEA::runTest(mulea_ranked_model_cervix_fc)
mulea_sub_reshaped_results_cervix_fc <- MulEA::reshapeResults(
  mulea_model = mulea_ranked_model_cervix_fc, 
  mulea_model_resuts = mulea_sub_results_cervix_fc, 
  mulea_model_resuts_ontology_col_name='ontologyId')

warnings()

# Plot results
MulEA::plotBarplot(mulea_relaxed_resuts = mulea_sub_reshaped_results_cervix_fc, 
                   statistics_value_cutoff=0.3)
MulEA::plotHeatmap(mulea_relaxed_resuts=mulea_sub_reshaped_results_cervix_fc, 
                   statistics_value_cutoff=0.3)
MulEA::plotGraph(mulea_relaxed_resuts=mulea_sub_reshaped_results_cervix_fc, 
                 statistics_value_cutoff = 0.3)

# Perform GSEA Cervix netrank
mulea_ranked_model_cervix_netrank <- MulEA::RankedBasedTest(
  gmt = input_gmt, testData = input_select_cervix_netrank, scores = input_select_scores_cervix_netrank, scoreType = "pos")
mulea_sub_results_cervix_netrank <- MulEA::runTest(mulea_ranked_model_cervix_netrank)
mulea_sub_reshaped_results_cervix_netrank <- MulEA::reshapeResults(
  mulea_model = mulea_ranked_model_cervix_netrank, 
  mulea_model_resuts = mulea_sub_results_cervix_netrank, 
  mulea_model_resuts_ontology_col_name='ontologyId')

warnings()

# Plot results
MulEA::plotBarplot(mulea_relaxed_resuts = mulea_sub_reshaped_results_cervix_netrank, 
                   statistics_value_cutoff=0.00000001)
MulEA::plotHeatmap(mulea_relaxed_resuts=mulea_sub_reshaped_results_cervix_netrank, 
                   statistics_value_cutoff=0.00000001)
MulEA::plotGraph(mulea_relaxed_resuts=mulea_sub_reshaped_results_cervix_netrank, 
                 statistics_value_cutoff = 0.00000001)


# Perform GSEA for Kidney
mulea_ranked_model_kidney_fc <- MulEA::RankedBasedTest(
  gmt = input_gmt, testData = input_select_kidney_fc, scores = input_select_scores_kidney_fc, scoreType = "pos")
mulea_sub_results_kidney_fc <- MulEA::runTest(mulea_ranked_model_kidney_fc)
mulea_sub_reshaped_results_kidney_fc <- MulEA::reshapeResults(
  mulea_model = mulea_ranked_model_kidney_fc, 
  mulea_model_resuts = mulea_sub_results_kidney_fc, 
  mulea_model_resuts_ontology_col_name='ontologyId')

warnings()

# Plot results
MulEA::plotBarplot(mulea_relaxed_resuts = mulea_sub_reshaped_results_kidney_fc, 
                   statistics_value_cutoff=0.001)
MulEA::plotHeatmap(mulea_relaxed_resuts=mulea_sub_reshaped_results_kidney_fc, 
                   statistics_value_cutoff=0.001)
MulEA::plotGraph(mulea_relaxed_resuts=mulea_sub_reshaped_results_kidney_fc, 
                 statistics_value_cutoff = 0.001)

# Perform GSEA for Kidney netrank
mulea_ranked_model_kidney_netrank <- MulEA::RankedBasedTest(
  gmt = input_gmt, testData = input_select_kidney_netrank, scores = input_select_scores_kidney_netrank, scoreType = "pos")
mulea_sub_results_kidney_netrank <- MulEA::runTest(mulea_ranked_model_kidney_netrank)
mulea_sub_reshaped_results_kidney_netrank <- MulEA::reshapeResults(
  mulea_model = mulea_ranked_model_kidney_netrank, 
  mulea_model_resuts = mulea_sub_results_kidney_netrank, 
  mulea_model_resuts_ontology_col_name='ontologyId')

warnings()

# Plot results
MulEA::plotBarplot(mulea_relaxed_resuts = mulea_sub_reshaped_results_kidney_netrank, 
                   statistics_value_cutoff=0.00000001)
MulEA::plotHeatmap(mulea_relaxed_resuts=mulea_sub_reshaped_results_kidney_netrank, 
                   statistics_value_cutoff=0.00000001)
MulEA::plotGraph(mulea_relaxed_resuts=mulea_sub_reshaped_results_kidney_netrank, 
                 statistics_value_cutoff = 0.00000001)



# DONOTRUN : Takes time, but shows interesting warnings. 
MulEA::plotGraph(mulea_relaxed_resuts=mulea_sub_reshaped_results, 
                 statistics_value_cutoff = 1.00)


# Start GSEA (Broad Institute) and read results.
# gsea-cli.bat GSEAPreranked -gmx D:\projects\science\MulEA\inst\tests\inputs\KEGG_Paths2geneSymbols.gmt -collapse false -mode Max_probe -norm meandiv -out D:\projects\science\MulEA\inst\tests\outputs -nperm 1000 -rnk D:\projects\science\MulEA\inst\tests\inputs\lung_FoldChange.rnk
# Lung
edb_file_lung_fc <- 'D:/projects/science/MulEA/inst/tests/outputs/my_analysis.GseaPreranked.1614331001532/edb/results.edb'
gsea_res_lung_fc <- MulEA::readEdbFileAsDataFrame(edbFilePath = edb_file_lung_fc)

edb_file_lung_netrank <- 'D:/projects/science/MulEA/inst/tests/outputs/my_analysis.GseaPreranked.1614331028430/edb/results.edb'
gsea_res_lung_netrank <- MulEA::readEdbFileAsDataFrame(edbFilePath = edb_file_lung_netrank)

# Cervix
edb_file_cervix_fc <- 'D:/projects/science/MulEA/inst/tests/outputs/my_analysis.GseaPreranked.1614331071925/edb/results.edb'
gsea_res_cervix_fc <- MulEA::readEdbFileAsDataFrame(edbFilePath = edb_file_cervix_fc)

edb_file_cervix_netrank <- 'D:/projects/science/MulEA/inst/tests/outputs/my_analysis.GseaPreranked.1614331099384/edb/results.edb'
gsea_res_cervix_netrank <- MulEA::readEdbFileAsDataFrame(edbFilePath = edb_file_cervix_netrank)

# Kidney
edb_file_kidney_fc <- 'D:/projects/science/MulEA/inst/tests/outputs/my_analysis.GseaPreranked.1614331131010/edb/results.edb'
gsea_res_kidney_fc <- MulEA::readEdbFileAsDataFrame(edbFilePath = edb_file_kidney_fc)

edb_file_kidney_netrank <- 'D:/projects/science/MulEA/inst/tests/outputs/my_analysis.GseaPreranked.1614331158898/edb/results.edb'
gsea_res_kidney_netrank <- MulEA::readEdbFileAsDataFrame(edbFilePath = edb_file_kidney_netrank)


# TODO : Comparison performed by man. Be careful, names of variables could have typo.
# mulea_sub_results_lung_fc vs gsea_res_lung_fc
# mulea_sub_results_lung_netrank vs gsea_res_lung_netrank

# mulea_sub_results_cervix_fc vs gsea_res_cervix_fc
# mulea_sub_results_cervix_netrank vs gsea_res_cervix_netrank

# mulea_sub_results_kidney_fc vs gsea_res_kidney_fc
# mulea_sub_results_kidney_netrank vs gsea_res_kidney_netrank


