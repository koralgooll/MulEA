
# Perform noise vs method.
library(MulEA)
library(tidyverse)

# Read tests for different ratios
mult_tests_01 <- read_rds("dev\\tests_res-5-05-01-1000.rds")
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
  geom_smooth(method='lm', formula= y~x) +
  facet_wrap(~method) +
  coord_cartesian(xlim=c(0,1), ylim=c(0,1))

comp_mult_tests %>% ggplot(aes(FPR, TPR)) + 
  geom_point() +
  geom_smooth(method='loess', formula= y~x+x^2) +
  facet_wrap(~method) +
  coord_cartesian(xlim=c(0,1), ylim=c(0,1))



# ROC on average 3.5 noise_ratio.
mult_tests_03_04 <- c(mult_tests_03, mult_tests_04)
mult_tests_03_04 <- c(mult_tests_01, mult_tests_02, mult_tests_03, 
                      mult_tests_04, mult_tests_05)


cut_off_range <- seq(0, 1, 0.005)

mult_tests_03_04_sum <- NULL
for (cut_off in cut_off_range) {
  mult_tests_03_04_sum_p <- MulEA:::getMultipleTestsSummary(
    tests_res = mult_tests_03_04, comparison_col_name = 'pValue', 
    labels = list('method'='p', 'cut_off'=cut_off), cut_off = cut_off)
  
  mult_tests_03_04_sum_bh <- MulEA:::getMultipleTestsSummary(
    tests_res = mult_tests_03_04, comparison_col_name = 'adjustedPValue', 
    labels = list('method'='bh', 'cut_off'=cut_off), cut_off = cut_off)
  
  mult_tests_03_04_sum_pt <- MulEA:::getMultipleTestsSummary(
    tests_res = mult_tests_03_04, comparison_col_name = 'adjustedPValueEmpirical', 
    labels = list('method'='pt', 'cut_off'=cut_off), cut_off = cut_off)
  
  mult_tests_03_04_sum <- rbind(mult_tests_03_04_sum, mult_tests_03_04_sum_p, 
                                mult_tests_03_04_sum_pt, mult_tests_03_04_sum_bh)
}


# ROC curve plot.
mult_tests_03_04_sum_grouped <- mult_tests_03_04_sum %>% 
  group_by(method, cut_off) %>% 
  summarize(TPR_mean = mean(as.numeric(TPR)), FPR_mean = mean(as.numeric(FPR)))

mult_tests_03_04_sum_grouped %>% ggplot(aes(FPR_mean, TPR_mean, colour=method)) + 
  geom_point() +
  geom_step() +
  facet_wrap(~method) +
  geom_abline(color="black") +
  scale_fill_brewer(palette="PuBu")

mult_tests_03_04_sum_grouped %>% ggplot(aes(FPR_mean, TPR_mean, colour=method)) + 
  geom_point() +
  geom_step() +
  geom_abline(color="black") +
  scale_fill_brewer(palette="PuBu")


# TODO IMPORTANT : Analyze it!!
# IDEA : More enrichment points will produce better estimates of densities.
# IDEA : Use all ratios, not only 0.3 and 0.4.
mult_tests_03_04_sum %>% ggplot(aes(FPR, colour=method)) + 
  stat_ecdf() + 
  geom_abline(color="black")

mult_tests_03_04_sum %>% ggplot(aes(FPR, fill = method, colour = method)) + 
  geom_density(adjust = 4, alpha = 0.1)



mult_tests_03_04_sum %>% ggplot(aes(TPR, colour=method)) + 
  stat_ecdf() + 
  geom_abline(color="black")

mult_tests_03_04_sum %>% ggplot(aes(TPR, fill = method, colour = method)) + 
  geom_density(adjust = 1, alpha = 0.1)

mult_tests_03_04_sum %>% ggplot(aes(TPR, fill = method, colour = method)) + 
  geom_density(adjust = 9, alpha = 0.1)






