# Perform noise vs method.
library(MulEA)
library(tidyverse)



# TPR wrapped by method.
# Manipulate by jitter height for better graph.
# jitter height for small is 0.25
comp_mult_tests %>% ggplot(aes(x=noise_ratio, y=TPR, fill=noise_ratio)) + 
  geom_boxplot(alpha=0.5, fatten=2) +
  geom_violin(alpha=0.3) +
  facet_grid(~method) + #facet_wrap
  theme(legend.position="right") +
  scale_fill_brewer(palette="PuBuGn") + 
  xlab('Noise Ratio') + ylab('True Positive Rate') +
  ylim(c(0,1))

# TPR wrapped by ratio.
# jitter height for small is 0.25
comp_mult_tests %>% ggplot(aes(x=method, y=TPR, fill=method)) + 
  geom_boxplot(alpha=0.5, fatten=2) +
  geom_violin(alpha=0.3) +
  facet_grid(~noise_ratio) +
  theme(legend.position="right") +
  scale_fill_brewer(palette="PuBuGn") + 
  xlab('Noise Ratio') + ylab('True Positive Rate') +
  ylim(c(0,1))


# FPR wrapped by method.
comp_mult_tests %>% ggplot(aes(x=noise_ratio, y=FPR, fill=noise_ratio)) + 
  # geom_jitter(aes(colour=noise_ratio), size=0.4, height = 0.1) + 
  geom_boxplot(alpha=0.5, fatten=2) +
  geom_violin(alpha=0.3) +
  facet_grid(~method) + #facet_wrap
  theme(legend.position="right") +
  scale_fill_brewer(palette="PuBuGn") + 
  xlab('Noise Ratio') + ylab('False Positive Rate')
# +
#   ylim(c(0,0.30))

# FPR wrapped by ratio.
comp_mult_tests %>% ggplot(aes(x=method, y=FPR, fill=method)) + 
  # geom_jitter(aes(colour=noise_ratio), size=0.4, height = 0.1) + 
  geom_boxplot(alpha=0.5, fatten=2) +
  geom_violin(alpha=0.3) +
  facet_grid(~noise_ratio) +
  theme(legend.position="right") +
  scale_fill_brewer(palette="PuBuGn") + 
  xlab('Noise Ratio') + ylab('False Positive Rate')
# +
#   ylim(c(0,0.15))


# FDR wrapped by method.
comp_mult_tests %>% ggplot(aes(x=noise_ratio, y=FDR, fill=noise_ratio)) + 
  # geom_jitter(aes(colour=noise_ratio), size=0.4, height = 0.1) + 
  geom_boxplot(alpha=0.5, fatten=2) +
  geom_violin(alpha=0.3) +
  facet_grid(~method) + #facet_wrap
  theme(legend.position="right") +
  scale_fill_brewer(palette="PuBuGn") + 
  xlab('Noise Ratio') + ylab('False Discovery Rate')

# FDR wrapped by ratio.
comp_mult_tests %>% ggplot(aes(x=method, y=FDR, fill=method)) + 
  # geom_jitter(aes(colour=noise_ratio), size=0.4, height = 0.1) + 
  geom_boxplot(alpha=0.5, fatten=2) +
  geom_violin(alpha=0.3) +
  facet_grid(~noise_ratio) +
  theme(legend.position="right") +
  scale_fill_brewer(palette="PuBuGn") + 
  xlab('Noise Ratio') + ylab('False Discovery Rate')


# TODO : IMPORTANT!
# Not mean approach to ROC:
mult_tests_03_04_sum %>% ggplot(aes(FPR, TPR, colour=method)) + 
  geom_point() +
  # geom_step() +
  facet_wrap(~method) +
  geom_abline(color="black") +
  scale_fill_brewer(palette="PuBu")


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

mult_tests_03_04_sum_grouped %>% ggplot(aes(TPR_mean, FPR_mean, colour=method)) + 
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

mult_tests_03_04_sum_grouped %>% ggplot(aes(TPR_mean, FPR_mean, colour=method)) + 
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

mult_tests_03_04_sum %>% ggplot(aes(FPR, fill = method, colour = method)) + 
  geom_density(adjust = 10, alpha = 0.1)



mult_tests_03_04_sum %>% ggplot(aes(TPR, colour=method)) + 
  stat_ecdf() + 
  geom_abline(color="black")

mult_tests_03_04_sum %>% ggplot(aes(TPR, fill = method, colour = method)) + 
  geom_density(adjust = 1, alpha = 0.1)

mult_tests_03_04_sum %>% ggplot(aes(TPR, fill = method, colour = method)) + 
  geom_density(adjust = 9, alpha = 0.1)


mult_tests_03_04_sum %>% ggplot(aes(TPR+FPR, fill = method, colour = method)) + 
  geom_density(adjust = 8, alpha = 0.1)

mult_tests_03_04_sum %>% ggplot(aes(FPR/TPR, fill = method, colour = method)) + 
  geom_density(adjust = 8, alpha = 0.1)




