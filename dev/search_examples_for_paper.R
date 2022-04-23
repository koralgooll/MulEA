library(tidyverse)

sim_f_small <- sim_mult_tests_res_sum_small %>% 
  filter(cut_off == 0.05 & noise_ratio >=0.3) %>%
  group_by(input_select) %>%
  group_split()

sim_f_big <- sim_mult_tests_res_sum_big %>% 
  filter(cut_off == 0.05 & noise_ratio >=0.3) %>%
  group_by(input_select) %>%
  group_split()

# Small
nice_sim_f_small <- list()
all_bh_TPR_small <- c()
for (usecase in sim_f_small) {
  if (usecase$TPR[2] > usecase$TPR[3] &&
      usecase$FPR[1] > usecase$FPR[2] &&
      usecase$TPR[2] == 1 &&
      usecase$FPR[2] == 0) {
    nice_sim_f_small <- c(nice_sim_f_small, list(usecase))
    all_bh_TPR_small <<- c(all_bh_TPR_small, usecase$TPR[3])
  }
}

# Big
nice_sim_f_big <- list()
all_bh_TPR_big <- c()
for (usecase in sim_f_big) {
  if (usecase$TPR[2] > usecase$TPR[3] &&
      usecase$FPR[1] > usecase$FPR[2] &&
      usecase$TPR[2] > 0.2 && 
      usecase$FPR[2] < 0.002) {
    nice_sim_f_big <- c(nice_sim_f_big, list(usecase))
    all_bh_TPR_big <<- c(all_bh_TPR_big, usecase$TPR[3])
  }
}

quantile(all_bh_TPR_big, seq(0,1,0.1))
max(all_bh_TPR_big)
min(all_bh_TPR_big)

# For small BH TPR min is 0.5
nice_small_sim_f <- list()
all_p_FPR_small <- c()
for (usecase in nice_sim_f_small) {
  if (usecase$TPR[3] == min(all_bh_TPR_small)) {
    nice_small_sim_f <- c(nice_small_sim_f, list(usecase))
    all_p_FPR_small <<- c(all_p_FPR_small, usecase$FPR[1])
  }
}


# For big BH TPR min is 0.5
nice_big_sim_f <- list()
all_p_FPR_big <- c()
for (usecase in nice_sim_f_big) {
  if (usecase$TPR[3] < 0.02) {
    nice_big_sim_f <- c(nice_big_sim_f, list(usecase))
    all_p_FPR_big <<- c(all_p_FPR_big, usecase$FPR[1])
  }
}

max(all_p_FPR_big)

# cut_off = 0.05
# 0.3 ratio
small_ontology_1 <- sim_mult_tests_res_s[[
  nice_small_sim_f[[1]]$test_no[1]]]$test_data[
    ,c("ontologyId", "ontologyName", "listOfValues")]

saveDataFrameAsGmtFile(small_ontology_1, gmtFilePath = "dev\\small_ontology_1.gmt")
write(x=nice_small_sim_f[[1]]$input_select[[1]], "dev\\small_sample_1.txt")
write(nice_small_sim_f[[1]][["over_repr_terms"]][[1]], "dev\\small_over_repr_terms_1.txt")

small_ontology_1 <- readGmtFileAsDataFrame(gmtFilePath = "dev\\small_ontology_1.gmt")
small_sample_1 <- readLines("dev\\small_sample_1.txt")
small_over_repr_terms_1 <- readLines("dev\\small_over_repr_terms_1.txt")

# 0.4 ratio
small_ontology_2 <- sim_mult_tests_res_s[[
  nice_small_sim_f[[8]]$test_no[1]]]$test_data[
    ,c("ontologyId", "ontologyName", "listOfValues")]

saveDataFrameAsGmtFile(small_ontology_2, gmtFilePath = "dev\\small_ontology_2.gmt")
write(x=nice_small_sim_f[[8]]$input_select[[1]], "dev\\small_sample_2.txt")
write(nice_small_sim_f[[8]][["over_repr_terms"]][[1]], "dev\\small_over_repr_terms_2.txt")


# cut_off = 0.05
# 0.3 ratio
big_ontology_1 <- sim_mult_tests_res_b[[
  nice_big_sim_f[[10]]$test_no[1]]]$test_data[
    ,c("ontologyId", "ontologyName", "listOfValues")]

saveDataFrameAsGmtFile(big_ontology_1, gmtFilePath = "dev\\big_ontology_1.gmt")
write(x=nice_big_sim_f[[10]]$input_select[[1]], "dev\\big_sample_1.txt")
write(nice_big_sim_f[[10]][["over_repr_terms"]][[1]], "dev\\big_over_repr_terms_1.txt")

# 0.3 ratio
big_ontology_2 <- sim_mult_tests_res_b[[
  nice_big_sim_f[[6]]$test_no[1]]]$test_data[
    ,c("ontologyId", "ontologyName", "listOfValues")]

saveDataFrameAsGmtFile(big_ontology_2, gmtFilePath = "dev\\big_ontology_2.gmt")
write(x=nice_big_sim_f[[6]]$input_select[[1]], "dev\\big_sample_2.txt")
write(nice_big_sim_f[[6]][["over_repr_terms"]][[1]], "dev\\big_over_repr_terms_2.txt")



# Extra
readr::read_rds("dev\\sim_mult_tests_small_sample_1.rds")

readr::read_rds("dev\\sim_mult_tests_big_over_repr_terms_1.rds")
readr::read_rds("dev\\sim_mult_tests_big_over_repr_terms_2.rds")


