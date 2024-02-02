library(tidyverse)

# Read and filter inputs (small).
mulea_path <- "/home/cezary/science/mulea/mulea"
gmtFilePath <- paste(mulea_path, 
                     "/dev/wikipathways_entrez_Mus_musculus_Marton.gmt", sep = "")
input_gmt <- mulea::read_gmt(gmtFilePath)



input_gmt_filtered <- mulea::filter_ontology(gmt = input_gmt)
filteredGmtFilePath <- paste(mulea_path, 
                             "/dev/wikipathways_filtered.gmt", sep = "")
mulea::write_gmt(gmt = input_gmt_filtered, file = filteredGmtFilePath)
input_gmt_small_filtered <- mulea::read_gmt(file = filteredGmtFilePath)

# Read and filter inputs (big).
gmtFilePath <- paste(mulea_path, 
                     "/dev/Pfam_Uniprot_Marton_Homo_sapiens.gmt", sep = "")
input_gmt <- mulea::read_gmt(gmtFilePath)
input_gmt_filtered <- mulea::filter_ontology(gmt = input_gmt)
filteredGmtFilePath <- paste(mulea_path, 
                             "/dev/Pfam_Uniprot_filtered.gmt", sep = "")
mulea::write_gmt(gmt = input_gmt_filtered, file = filteredGmtFilePath)
input_gmt_big_filtered <- mulea::read_gmt(file = filteredGmtFilePath)

# Tests calculations
sim_mult_tests_res <- mulea:::simulateMultipleTestsWithRatioParam(
  input_gmt_filtered = input_gmt_small_filtered,
  noise_ratio_range = seq(0.00, 0.15, 0.05),
  over_repr_ratio = 0.85,
  number_of_tests = 3, 
  nthreads = 16)
sim_mult_tests_res <- mulea:::simulateMultipleTestsWithRatioParam(
  input_gmt_filtered = input_gmt_small_filtered,
  noise_ratio_range = seq(0.00, 0.15, 0.05),
  over_repr_ratio = 0.85,
  number_of_tests = 3, 
  nthreads = 32)
sim_mult_tests_res <- mulea:::simulateMultipleTestsWithRatioParam(
  input_gmt_filtered = input_gmt_small_filtered,
  noise_ratio_range = seq(0.00, 0.15, 0.05),
  over_repr_ratio = 0.85,
  number_of_tests = 3, 
  nthreads = 12)
sim_mult_tests_res <- mulea:::simulateMultipleTestsWithRatioParam(
  input_gmt_filtered = input_gmt_small_filtered,
  noise_ratio_range = seq(0.00, 0.15, 0.05),
  over_repr_ratio = 0.85,
  number_of_tests = 3, 
  nthreads = 8)
sim_mult_tests_res <- mulea:::simulateMultipleTestsWithRatioParam(
  input_gmt_filtered = input_gmt_small_filtered,
  noise_ratio_range = seq(0.00, 0.15, 0.05),
  over_repr_ratio = 0.85,
  number_of_tests = 3, 
  nthreads = 4)
sim_mult_tests_res <- mulea:::simulateMultipleTestsWithRatioParam(
  input_gmt_filtered = input_gmt_small_filtered,
  noise_ratio_range = seq(0.00, 0.15, 0.05),
  over_repr_ratio = 0.85,
  number_of_tests = 3, 
  nthreads = 2)
sim_mult_tests_res_bench <- mulea:::simulateMultipleTestsWithRatioParam(
  input_gmt_filtered = input_gmt_small_filtered,
  noise_ratio_range = seq(0.00, 0.15, 0.05),
  over_repr_ratio = 0.85,
  number_of_tests = 3, 
  nthreads = 1)

# Real calculations
# Set some seed!
sim_mult_tests_res <- mulea:::simulateMultipleTestsWithRatioParam(
  input_gmt_filtered = input_gmt_small_filtered,
  noise_ratio_range = seq(0.00, 0.20, 0.05),
  over_repr_ratio = 0.85,
  number_of_tests = 1000, 
  nthreads = 1)

# sim_mult_tests_res <- mulea:::simulateMultipleTestsWithRatioParam(
#   input_gmt_filtered = input_gmt_big_filtered,
#   noise_ratio_range = seq(0.00, 0.25, 0.05),
#   over_repr_ratio = 0.85,
#   number_of_tests = 200, 
#   nthreads = 4)

# File desc meaning sim_res_{ontology_size}_{over_repr_ratio}_{number_of_tests}.rds
# sim_mult_tests_res <- input_gmt_small_filtered
# readr::write_rds(
#   x = sim_mult_tests_res,
#   file = paste(mulea_path, "/dev/sim_res_small_wiki_085_1000.rds", sep = ""))

sim_mult_tests_res <- readr::read_rds(
  file = paste(mulea_path, "/dev/sim_res_small_wiki_085_1000.rds", sep = ""))

sim_mult_tests_res <- readr::read_rds(
  file = paste(mulea_path, "/dev/sim_res_big_085_0-025_200.rds", sep = ""))


sim_mult_tests_res_sum <- mulea:::getMultipleTestsSummaryAcrossCutOff(
  tests_res=sim_mult_tests_res, cut_off_range = seq(0, 1, 0.001))

# readr::write_rds(
#   x = sim_mult_tests_res_sum,
#   file = paste(mulea_path, "/dev/sim_res_small_wiki_085_1000_sum.rds", sep = ""))

sim_mult_tests_res_sum <- readr::read_rds(
  file = paste(mulea_path, "/dev/sim_res_small_wiki_085_1000_sum.rds", sep = ""))

sim_mult_tests_res_sum <- readr::read_rds(
  file = paste(mulea_path, "/dev/sim_res_big_085_0-025_200_sum.rds", sep = ""))

sim_mult_tests_res_to_roc <- mulea:::getSummaryToRoc(tests_res = sim_mult_tests_res)

sim_mult_tests_res_to_roc <- readr::read_rds(
  file = paste(mulea_path, "/dev/sim_res_big_085_200_roc.rds", sep = ""))


readr::write_rds(
  x = sim_mult_tests_res_to_roc,
  file = paste(mulea_path, "/dev/sim_res_small_wiki_085_1000_roc.rds", sep = ""))

plot_res <- sim_mult_tests_res_sum %>% ggplot(aes(x=noise_ratio, y=TPR, fill=noise_ratio)) + 
  geom_boxplot(alpha=0.5, fatten=2) +
  geom_violin(alpha=0.3) +
  facet_grid(~method) + #facet_wrap
  theme(legend.position="right") +
  scale_fill_brewer(palette="PuBuGn") + 
  xlab('Noise Ratio') + ylab('True Positive Rate') +
  ylim(c(0,1))





sim_mult_tests_res <- mulea:::simulateMultipleTestsWithRatioParam(
  input_gmt_filtered = input_gmt_small_filtered,
  noise_ratio_range = seq(0.00, 0.15, 0.05),
  over_repr_ratio = 0.85,
  number_of_tests = 3, 
  nthreads = 16)
