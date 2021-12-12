library(tidyverse)


sim_f <- sim_mult_tests_res_sum %>% 
  filter(cut_off == 0.05 & noise_ratio >=0.3) %>%
  group_by(input_select) %>%
  group_split() 

# Small
nice_sim_f <- list()
all_bh_TPR <- c()
for (usecase in sim_f) {
  if (usecase$TPR[2] > usecase$TPR[3] &&
      usecase$FPR[1] > usecase$FPR[2] &&
      usecase$TPR[2] == 1 &&
      usecase$FPR[2] == 0) {
    nice_sim_f <- c(nice_sim_f, list(usecase))
    all_bh_TPR <<- c(all_bh_TPR, usecase$TPR[3])
  }
}

# Big
nice_sim_f <- list()
all_bh_TPR <- c()
for (usecase in sim_f) {
  if (usecase$TPR[2] > usecase$TPR[3] &&
      usecase$FPR[1] > usecase$FPR[2] &&
      usecase$TPR[2] > 0.2 && 
      usecase$FPR[2] < 0.002) {
    nice_sim_f <- c(nice_sim_f, list(usecase))
    all_bh_TPR <<- c(all_bh_TPR, usecase$TPR[3])
  }
}

quantile(all_bh_TPR, seq(0,1,0.1))
max(all_bh_TPR)
min(all_bh_TPR)

# For small BH TPR min is 0.5
nice_small_sim_f <- list()
all_p_FPR <- c()
for (usecase in nice_sim_f) {
  if (usecase$TPR[3] == min(all_bh_TPR)) {
    nice_small_sim_f <- c(nice_small_sim_f, list(usecase))
    all_p_FPR <<- c(all_p_FPR, usecase$FPR[1])
  }
}


# For big BH TPR min is 0.5
nice_small_sim_f <- list()
all_p_FPR <- c()
for (usecase in nice_sim_f) {
  if (usecase$TPR[3] < 0.02) {
    nice_small_sim_f <- c(nice_small_sim_f, list(usecase))
    all_p_FPR <<- c(all_p_FPR, usecase$FPR[1])
  }
}

max(all_p_FPR)

# cut_off = 0.05
# 0.3 ratio
write_rds(nice_small_sim_f[[1]]$input_select[[1]], "dev\\sim_mult_tests_small_sample_1.rds")
# 0.4 ratio
write_rds(nice_small_sim_f[[8]]$input_select[[1]], "dev\\sim_mult_tests_small_sample_2.rds")
readr::read_rds("dev\\sim_mult_tests_small_sample_1.rds")


# cut_off = 0.05
# 0.3 ratio
write_rds(nice_small_sim_f[[10]]$input_select[[1]], "dev\\sim_mult_tests_big_sample_1.rds")
# 0.3 ratio
write_rds(nice_small_sim_f[[6]]$input_select[[1]], "dev\\sim_mult_tests_big_sample_2.rds")
readr::read_rds("dev\\sim_mult_tests_small_sample_1.rds")




