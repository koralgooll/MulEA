library(tidyverse)

# Read and filter inputs (small).
mulea_path <- "D:/projects/science/MulEA_backup_2022_05_31/MulEA"
gmtFilePath <- paste(mulea_path, 
                     "/inst/tests/inputs/KEGG_example_c_elegans_Leila.gmt", sep = "")
input_gmt <- MulEA::read_gmt(gmtFilePath)
input_gmt_filtered <- MulEA::filter_ontology(gmt = input_gmt)
filteredGmtFilePath <- paste(mulea_path, 
                             "/inst/tests/inputs/KEGG_filtered.gmt", sep = "")
MulEA::write_gmt(gmt = input_gmt_filtered, file = filteredGmtFilePath)
input_gmt_small_filtered <- MulEA::read_gmt(file = filteredGmtFilePath)

# Read and filter inputs (big).
gmtFilePath <- paste(mulea_path, 
                     "/inst/tests/inputs/Pfam_Uniprot_Marton_Homo_sapiens.gmt", sep = "")
input_gmt <- MulEA::read_gmt(gmtFilePath)
input_gmt_filtered <- MulEA::filter_ontology(gmt = input_gmt)
filteredGmtFilePath <- paste(mulea_path, 
                             "/inst/tests/inputs/Pfam_Uniprot_filtered.gmt", sep = "")
MulEA::write_gmt(gmt = input_gmt_filtered, file = filteredGmtFilePath)
input_gmt_big_filtered <- MulEA::read_gmt(file = filteredGmtFilePath)


sim_mult_tests_res <- MulEA:::simulateMultipleTestsWithRatioParam(
  input_gmt_filtered = input_gmt_small_filtered,
  noise_ratio_range = seq(0.00, 0.15, 0.05),
  over_repr_ratio = 0.85,
  number_of_tests = 3, nthreads = 16)

sim_mult_tests_res_sum <- MulEA:::getMultipleTestsSummaryAcrossCutOff(
  tests_res=sim_mult_tests_res)

plot_res <- sim_mult_tests_res_sum %>% ggplot(aes(x=noise_ratio, y=TPR, fill=noise_ratio)) + 
  geom_boxplot(alpha=0.5, fatten=2) +
  geom_violin(alpha=0.3) +
  facet_grid(~method) + #facet_wrap
  theme(legend.position="right") +
  scale_fill_brewer(palette="PuBuGn") + 
  xlab('Noise Ratio') + ylab('True Positive Rate') +
  ylim(c(0,1))


