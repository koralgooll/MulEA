
# Perform noise vs method.
library(MulEA)
library(tidyverse)

# Read tests for different ratios
mult_tests_01 <- read_rds("tests_res-5-05-01-1000.rds")
mult_tests_02 <- read_rds("tests_res-5-05-02-1000.rds")
mult_tests_03 <- read_rds("tests_res-5-05-03-1000.rds")
mult_tests_04 <- read_rds("tests_res-5-05-04-1000.rds")
mult_tests_05 <- read_rds("tests_res-5-05-05-1000.rds")

# Summarize by pValue
mult_tests_01_sum_p <- MulEA:::getMultipleTestsSummary(
  tests_res = mult_tests_01, comparison_col_name = 'pValue', 
  labels = list('noise_ratio'=0.1, 'method'='p'))
mult_tests_02_sum_p <- MulEA:::getMultipleTestsSummary(
  tests_res = mult_tests_02, comparison_col_name = 'pValue', 
  labels = list('noise_ratio'=0.2, 'method'='p'))
mult_tests_03_sum_p <- MulEA:::getMultipleTestsSummary(
  tests_res = mult_tests_03, comparison_col_name = 'pValue', 
  labels = list('noise_ratio'=0.3, 'method'='p'))
mult_tests_04_sum_p <- MulEA:::getMultipleTestsSummary(
  tests_res = mult_tests_04, comparison_col_name = 'pValue', 
  labels = list('noise_ratio'=0.4, 'method'='p'))
mult_tests_05_sum_p <- MulEA:::getMultipleTestsSummary(
  tests_res = mult_tests_05, comparison_col_name = 'pValue', 
  labels = list('noise_ratio'=0.5, 'method'='p'))

comp_mult_tests_p <- rbind(mult_tests_01_sum_p, mult_tests_02_sum_p, 
                           mult_tests_03_sum_p, mult_tests_04_sum_p, 
                           mult_tests_05_sum_p)


mult_tests_01_sum_pt <- MulEA:::getMultipleTestsSummary(
  tests_res = mult_tests_01, comparison_col_name = 'adjustedPValueEmpirical', 
  labels = list('noise_ratio'=0.1, 'method'='pt'))
mult_tests_02_sum_pt <- MulEA:::getMultipleTestsSummary(
  tests_res = mult_tests_02, comparison_col_name = 'adjustedPValueEmpirical', 
  labels = list('noise_ratio'=0.2, 'method'='pt'))
mult_tests_03_sum_pt <- MulEA:::getMultipleTestsSummary(
  tests_res = mult_tests_03, comparison_col_name = 'adjustedPValueEmpirical', 
  labels = list('noise_ratio'=0.3, 'method'='pt'))
mult_tests_04_sum_pt <- MulEA:::getMultipleTestsSummary(
  tests_res = mult_tests_04, comparison_col_name = 'adjustedPValueEmpirical', 
  labels = list('noise_ratio'=0.4, 'method'='pt'))
mult_tests_05_sum_pt <- MulEA:::getMultipleTestsSummary(
  tests_res = mult_tests_05, comparison_col_name = 'adjustedPValueEmpirical', 
  labels = list('noise_ratio'=0.5, 'method'='pt'))

comp_mult_tests_pt <- rbind(mult_tests_01_sum_pt, mult_tests_02_sum_pt, 
                           mult_tests_03_sum_pt, mult_tests_04_sum_pt, 
                           mult_tests_05_sum_pt)


mult_tests_01_sum_bh <- MulEA:::getMultipleTestsSummary(
  tests_res = mult_tests_01, comparison_col_name = 'adjustedPValue', 
  labels = list('noise_ratio'=0.1, 'method'='bh'))
mult_tests_02_sum_bh <- MulEA:::getMultipleTestsSummary(
  tests_res = mult_tests_02, comparison_col_name = 'adjustedPValue', 
  labels = list('noise_ratio'=0.2, 'method'='bh'))
mult_tests_03_sum_bh <- MulEA:::getMultipleTestsSummary(
  tests_res = mult_tests_03, comparison_col_name = 'adjustedPValue', 
  labels = list('noise_ratio'=0.3, 'method'='bh'))
mult_tests_04_sum_bh <- MulEA:::getMultipleTestsSummary(
  tests_res = mult_tests_04, comparison_col_name = 'adjustedPValue', 
  labels = list('noise_ratio'=0.4, 'method'='bh'))
mult_tests_05_sum_bh <- MulEA:::getMultipleTestsSummary(
  tests_res = mult_tests_05, comparison_col_name = 'adjustedPValue', 
  labels = list('noise_ratio'=0.5, 'method'='bh'))

comp_mult_tests_bh <- rbind(mult_tests_01_sum_bh, mult_tests_02_sum_bh, 
                            mult_tests_03_sum_bh, mult_tests_04_sum_bh, 
                            mult_tests_05_sum_bh)

comp_mult_tests <- rbind(comp_mult_tests_p, comp_mult_tests_pt, comp_mult_tests_bh)


# TPR wrapped by method.
comp_mult_tests %>% ggplot(aes(x=noise_ratio, y=TPR, fill=noise_ratio)) + 
  geom_boxplot(alpha=0.5, fatten=2) +
  geom_violin(alpha=0.3) +
  geom_jitter(aes(colour=noise_ratio), size=0.4, height = 0.1) + 
  facet_grid(~method) + #facet_wrap
  theme(legend.position="right") +
  scale_fill_brewer(palette="PuBuGn") + 
  xlab('Noise Ratio') + ylab('True Positive Rate') +
  ylim(c(0,1))

# TPR wrapped by ratio.
comp_mult_tests %>% ggplot(aes(x=method, y=TPR, fill=method)) + 
  geom_boxplot(alpha=0.5, fatten=2) +
  geom_violin(alpha=0.3) +
  geom_jitter(aes(colour=noise_ratio), size=0.4, height = 0.1) + 
  facet_grid(~noise_ratio) +
  theme(legend.position="right") +
  scale_fill_brewer(palette="PuBuGn") + 
  xlab('Noise Ratio') + ylab('True Positive Rate') +
  ylim(c(0,1))

# FPR wrapped by method.
comp_mult_tests %>% ggplot(aes(x=noise_ratio, y=FPR, fill=noise_ratio)) + 
  geom_boxplot(alpha=0.5, fatten=2) +
  geom_violin(alpha=0.3) +
  geom_jitter(aes(colour=noise_ratio), size=0.4, height = 0.1) + 
  facet_grid(~method) + #facet_wrap
  theme(legend.position="right") +
  scale_fill_brewer(palette="PuBuGn") + 
  xlab('Noise Ratio') + ylab('False Positive Rate') +
  ylim(c(0,1))

# FPR wrapped by ratio.
comp_mult_tests %>% ggplot(aes(x=method, y=FPR, fill=method)) + 
  geom_boxplot(alpha=0.5, fatten=2) +
  geom_violin(alpha=0.3) +
  geom_jitter(aes(colour=noise_ratio), size=0.4, height = 0.1) + 
  facet_grid(~noise_ratio) +
  theme(legend.position="right") +
  scale_fill_brewer(palette="PuBuGn") + 
  xlab('Noise Ratio') + ylab('False Positive Rate') +
  ylim(c(0,1))

# ROC by ratio.
comp_mult_tests %>% ggplot(aes(FPR, TPR)) + 
  geom_point() +
  geom_smooth(method='lm', formula= y~x+x^2) +
  coord_cartesian(xlim=c(0,1), ylim=c(0,1))

comp_mult_tests %>% ggplot(aes(FPR, TPR)) + 
  geom_point() +
  geom_smooth(method='loess', formula= y~x+x^2) +
  coord_cartesian(xlim=c(0,1), ylim=c(0,1))



simulateMultipleTests <- function(
  input_gmt_filtered, number_of_tests = 10, 
  noise_ratio = 0.2,
  number_of_over_representation_groups = ceiling(nrow(input_gmt_filtered)*0.1), 
  number_of_under_representation_groups = 0) {
  
  print("Mulea calculation time:")
  tic()
  tests_res <- vector("list", number_of_tests)
  for (i in 1:number_of_tests) {
    print(i)
    
    input_gmt_decorated <- MulEA:::decorateGmtByUnderOvenAndNoise(
      input_gmt = input_gmt_filtered,
      number_of_over_representation_groups = number_of_over_representation_groups,
      number_of_under_representation_groups = number_of_under_representation_groups)
    
    samples <- MulEA:::generateInputSamples(
      input_gmt_decorated, 
      noise_ratio=noise_ratio, 
      over_repr_ratio=over_repr_ratio, 
      number_of_samples=number_of_samples)
    
    if (length(samples) != 1) {
      warning("sample is not size 1")
    }
    
    input_select <- unlist(samples)
    
    mulea_ora_model <- MulEA::ORA(
      gmt = input_gmt_filtered, testData = input_select, adjustMethod = "PT",
      numberOfPermutations = number_of_steps, nthreads = 16)
    
    mulea_ora_results <- MulEA::runTest(mulea_ora_model)
    tests_res[[i]]$mulea_res <- mulea_ora_results
    tests_res[[i]]$test_data <- input_gmt_decorated
  }
  toc()
}

