library(tidyverse)
library(MulEA)

sim_mult_tests_res <- readr::read_rds("dev\\new_tests_res\\sim_mult_tests_res_small_075_100.rds")
set_name <- "small_75"
sim_mult_tests_res <- readr::read_rds("dev\\new_tests_res\\sim_mult_tests_res_small_085_100.rds")
set_name <- "small_85"
sim_mult_tests_res <- readr::read_rds("dev\\new_tests_res\\sim_mult_tests_res_small_095_100.rds")
set_name <- "small_95"
sim_mult_tests_res <- readr::read_rds("dev\\new_tests_res\\sim_mult_tests_res_big_075_30.rds")
set_name <- "big_75"
sim_mult_tests_res <- readr::read_rds("dev\\new_tests_res\\sim_mult_tests_res_big_085_30.rds")
set_name <- "big_85"
sim_mult_tests_res <- readr::read_rds("dev\\new_tests_res\\sim_mult_tests_res_big_095_30.rds")
set_name <- "big_95"

sim_mult_tests_res_sum <- MulEA:::getMultipleTestsSummaryAcrossCutOff(tests_res=sim_mult_tests_res)
sim_mult_tests_res_to_roc <- MulEA:::getSummaryToRoc(tests_res = sim_mult_tests_res)

# Load saved data:
mulea_path <- "/home/cezary/science/MulEA/MulEA"
sim_mult_tests_res_sum <- readr::read_rds(paste(mulea_path, "/dev/sim_res_small_085_1000_sum.rds", sep = ""))
set_name <- "small_85"
sim_mult_tests_res_sum <- readr::read_rds(paste(mulea_path, "/dev/sim_res_big_085_200_sum.rds", sep = ""))
set_name <- "big_85"

# TPR:
# OK - TPR wrapped by method.
plot_res <- sim_mult_tests_res_sum %>% ggplot(aes(x=noise_ratio, y=TPR, fill=noise_ratio)) + 
  geom_boxplot(alpha=0.5, fatten=2) +
  geom_violin(alpha=0.3) +
  facet_grid(~forcats::fct_relevel(method, 'p', 'bh', 'pt'), labeller = as_labeller(
    c(`bh` = "Benjamini-Hochberg", `p` = "Uncorrected p-value", `pt` = "eFDR"))) + #facet_wrap
  theme(legend.position="right") +
  scale_fill_brewer(palette="Set1", guide="none") + 
  xlab('Noise Ratio') + ylab('True Positive Rate') +
  ylim(c(0,1))
ggsave(plot_res, filename = paste0(set_name, "_TPR_method", ".jpeg"), device = "jpeg", 
       path = paste(mulea_path, "/dev/new_plots", sep = ""))


# OK - TPR wrapped by ratio.
plot_res <- sim_mult_tests_res_sum %>% ggplot(
  aes(x=forcats::fct_relevel(method, 'p', 'bh', 'pt'), y=TPR, 
      fill=forcats::fct_relevel(method, 'p', 'bh', 'pt'))) + 
  geom_boxplot(alpha=0.5, fatten=2) +
  geom_violin(alpha=0.3) +
  facet_grid(~noise_ratio) +
  theme(legend.position="bottom", 
        axis.ticks.x=element_blank(), axis.text.x=element_blank()) +
  guides(fill=guide_legend(title='Method')) +
  scale_fill_discrete(palette = scales::hue_pal(l = 70), 
                      labels = c("Uncorrected p-value", "Benjamini-Hochberg", "eFDR")) +
  scale_x_discrete(position = "top") +
  xlab('Noise Ratio') + ylab('True Positive Rate') +
  ylim(c(0,1))
ggsave(plot_res, filename = paste0(set_name, "_TPR_ratio", ".jpeg"), device = "jpeg", 
       path = paste(mulea_path, "/dev/new_plots", sep = ""))


# FPR:
# OK - FPR wrapped by method.
plot_res <- sim_mult_tests_res_sum %>% ggplot(aes(x=noise_ratio, y=FPR, fill=noise_ratio)) + 
  geom_boxplot(alpha=0.5, fatten=2) +
  geom_violin(alpha=0.3) +
  facet_grid(~forcats::fct_relevel(method, 'p', 'bh', 'pt'), labeller = as_labeller(
    c(`bh` = "Benjamini-Hochberg", `p` = "Uncorrected p-value", `pt` = "eFDR"))) + #facet_wrap
  theme(legend.position="right") +
  scale_fill_brewer(palette="Set1", guide="none") + 
  xlab('Noise Ratio') + ylab('False Positive Rate') +
  ylim(c(0,1))
ggsave(plot_res, filename = paste0(set_name, "_FPR_method", ".jpeg"), device = "jpeg", 
       path = paste(mulea_path, "/dev/new_plots", sep = ""))


# OK - FPR wrapped by ratio.
plot_res <- sim_mult_tests_res_sum %>% ggplot(
  aes(x=forcats::fct_relevel(method, 'p', 'bh', 'pt'), y=FPR, 
      fill=forcats::fct_relevel(method, 'p', 'bh', 'pt'))) +
  geom_boxplot(alpha=0.5, fatten=2) +
  geom_violin(alpha=0.3) +
  facet_grid(~noise_ratio) +
  theme(legend.position="bottom", 
        axis.ticks.x=element_blank(), axis.text.x=element_blank()) +
  guides(fill=guide_legend(title='Method')) +
  scale_fill_discrete(palette = scales::hue_pal(l = 70), 
                      labels = c("Uncorrected p-value", "Benjamini-Hochberg", "eFDR")) +
  scale_x_discrete(position = "top") +
  xlab('Noise Ratio') + ylab('False Positive Rate') +
  ylim(c(0,1))
ggsave(plot_res, filename = paste0(set_name, "_FPR_ratio", ".jpeg"), device = "jpeg", 
       path = paste(mulea_path, "/dev/new_plots", sep = ""))


# FDR:
# OK - FDR wrapped by method.
plot_res <- sim_mult_tests_res_sum %>% ggplot(aes(x=noise_ratio, y=FDR, fill=noise_ratio)) +
  geom_boxplot(alpha=0.5, fatten=2) +
  geom_violin(alpha=0.3) +
  facet_grid(~forcats::fct_relevel(method, 'p', 'bh', 'pt'), labeller = as_labeller(
    c(`bh` = "Benjamini-Hochberg", `p` = "Uncorrected p-value", `pt` = "eFDR"))) + #facet_wrap
  theme(legend.position="right") +
  scale_fill_brewer(palette="Set1", guide="none") + 
  xlab('Noise Ratio') + ylab('False Discovery Rate') +
  ylim(c(0,1))
ggsave(plot_res, filename = paste0(set_name, "_FDR_method", ".jpeg"), device = "jpeg", 
       path = paste(mulea_path, "/dev/new_plots", sep = ""))

# FDR wrapped by ratio.
plot_res <- sim_mult_tests_res_sum %>% ggplot(
  aes(x=forcats::fct_relevel(method, 'p', 'bh', 'pt'), y=FDR, 
      fill=forcats::fct_relevel(method, 'p', 'bh', 'pt'))) +
  geom_boxplot(alpha=0.5, fatten=2) +
  geom_violin(alpha=0.3) +
  facet_grid(~noise_ratio) +
  theme(legend.position="bottom", 
        axis.ticks.x=element_blank(), axis.text.x=element_blank()) +
  guides(fill=guide_legend(title='Method')) +
  scale_fill_discrete(palette = scales::hue_pal(l = 70), 
                      labels = c("Uncorrected p-value", "Benjamini-Hochberg", "eFDR")) +
  scale_x_discrete(position = "top") +
  xlab('Noise Ratio') + ylab('False Discovery Rate') +
  ylim(c(0,1))
ggsave(plot_res, filename = paste0(set_name, "_FDR_ratio", ".jpeg"), device = "jpeg", 
       path = paste(mulea_path, "/dev/new_plots", sep = ""))


# ROC curve plot
plot_res <- sim_mult_tests_res_to_roc %>% ggplot(
  aes(FPR, TPR, 
      colour=forcats::fct_relevel(method, 'p_value', 'adjusted_p_value', 'eFDR'))) + 
  # geom_point() +
  geom_step() +
  theme(legend.position="bottom") +
  # facet_wrap(~method) +
  guides(colour=guide_legend(title='Method')) +
  scale_colour_discrete(palette = scales::hue_pal(l = 70), 
                      labels = c("Uncorrected p-value", "Benjamini-Hochberg", "eFDR")) +
  geom_abline(color="black")
ggsave(plot_res, filename = paste0(set_name, "_ROC", ".jpeg"), device = "jpeg", 
       path = paste(mulea_path, "/dev/new_plots", sep = ""))


# ROC density curve plot
plot_res <- sim_mult_tests_res_to_roc %>% ggplot(
  aes(FPR, TPR, 
      colour=forcats::fct_relevel(method, 'p_value', 'adjusted_p_value', 'eFDR'))) + 
  geom_point() +
  # geom_step() +
  theme(legend.position="bottom") +
  # facet_wrap(~method) +
  facet_wrap(~forcats::fct_relevel(method, 'p_value', 'adjusted_p_value', 'eFDR'), 
             labeller = as_labeller(
    c(`adjusted_p_value` = "Benjamini-Hochberg", 
      `p_value` = "Uncorrected p-value", `eFDR` = "eFDR"))) +
  guides(colour=guide_legend(title='Method')) +
  scale_colour_discrete(palette = scales::hue_pal(l = 70), 
                        labels = c("Uncorrected p-value", "Benjamini-Hochberg", "eFDR")) +
  geom_abline(color="black")
ggsave(plot_res, filename = paste0(set_name, "_ROC_dens", ".jpeg"), device = "jpeg", 
       path = paste(mulea_path, "/dev/new_plots", sep = ""))

######################################################################################
#### END                                                                          ####
######################################################################################












# IN PROGRESS:
plot_res <- sim_mult_tests_res_sum %>% ggplot(aes(FPR, TPR, colour=method)) + 
  geom_point() +
  # geom_step() +
  facet_wrap(~method) +
  geom_abline(color="black") +
  scale_fill_brewer(palette="PuBu")


sim_mult_tests_res_sum %>% ggplot(aes(FPR, TPR, colour=method)) + 
  geom_point() +
  geom_smooth(method='lm', formula = y~poly(x, 2)) +
  # stat_smooth() +
  # stat_density() +
  facet_wrap(~method) +
  geom_abline(color="black") +
  scale_fill_brewer(palette="PuBu")


# DEBUG:
hist(as.numeric(sim_mult_tests_res_sum$cut_off), breaks = 20)
tests_res <- sim_mult_tests_res_sum
cut_off_range = seq(0, 1, 0.1)
cut_off <- 0.3

hist(as.numeric(sim_mult_tests_res_sum[sim_mult_tests_res_sum$method == 'p', ]$FPR), breaks = 20)
hist(as.numeric(sim_mult_tests_res_sum[sim_mult_tests_res_sum$method == 'p', ]$TPR), breaks = 20)



install.packages("ISLR")
install.packages("pROC")
install.packages("ROCR")
install.packages("PRROC")


#load Default dataset from ISLR book
data <- ISLR::Default

#divide dataset into training and test set
set.seed(1)
sample <- sample(c(TRUE, FALSE), nrow(data), replace=TRUE, prob=c(0.7,0.3))
train <- data[sample, ]
test <- data[!sample, ]

#fit logistic regression model to training set
model <- glm(default~student+balance+income, family="binomial", data=train)

#use model to make predictions on test set
predicted <- predict(model, test, type="response")

#load necessary packages
library(ggplot2)
library(pROC)

#define object to plot
rocobj <- roc(test$default, predicted)
sim_mult_tests_res_sum$over_repr_terms
#create ROC plot
ggroc(rocobj)


