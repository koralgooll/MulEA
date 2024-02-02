library(tidyverse)
library(mulea)


# If calculation of summaries and roc data is necessary.
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

sim_mult_tests_res_sum <- mulea:::getMultipleTestsSummaryAcrossCutOff(tests_res=sim_mult_tests_res)
sim_mult_tests_res_to_roc <- mulea:::getSummaryToRoc(tests_res = sim_mult_tests_res)


# If start from serialized data.
# Load saved data.
mulea_path <- "/home/cezary/science/mulea/mulea"

# Small data set (~100 ontologies).
sim_mult_tests_res_sum <- readr::read_rds(paste(mulea_path, "/dev/sim_res_small_wiki_085_1000_sum.rds", sep = ""))
sim_mult_tests_res_to_roc <- readr::read_rds(paste(mulea_path, "/dev/sim_res_small_wiki_085_1000_roc.rds", sep = ""))
set_name <- "small_85"

# Big data set (~3000 ontologies).
sim_mult_tests_res_sum <- readr::read_rds(paste(mulea_path, "/dev/sim_res_big_085_0-025_200_sum.rds", sep = ""))
sim_mult_tests_res_to_roc <- readr::read_rds(paste(mulea_path, "/dev/sim_res_big_085_0-025_200_roc.rds", sep = ""))
set_name <- "big_85"

# Take slice to plots. Plot with 15 000 000 point is slow and hard to read, do not really show properties.
sim_mult_tests_res_sum_slice <- sim_mult_tests_res_sum %>% slice_sample(n=10000)
slice_name <- "_n_10k"
sim_mult_tests_res_sum_slice <- sim_mult_tests_res_sum %>% filter(cut_off == 0.05)
slice_name <- "_cut_005"
sim_mult_tests_res_sum_slice <- sim_mult_tests_res_sum %>% filter(cut_off == 0.03)
slice_name <- "_cut_003"

# TPR:
# OK - TPR wrapped by method.
plot_res <- sim_mult_tests_res_sum_slice %>% ggplot(aes(x=noise_ratio, y=TPR, fill=noise_ratio)) + 
  geom_boxplot(alpha=0.5, fatten=2) +
  geom_violin(alpha=0.3) +
  facet_grid(~forcats::fct_relevel(method, 'p', 'bh', 'pt'), labeller = as_labeller(
    c(`bh` = "Benjamini-Hochberg", `p` = "Uncorrected p-value", `pt` = "eFDR"))) + #facet_wrap
  theme_bw() + 
  theme(legend.position="right") +
  scale_fill_brewer(palette="Set1", guide="none") + 
  xlab('Noise Ratio') + ylab('True Positive Rate') +
  ylim(c(0,1))
plot_desc <- "_TPR_method"
ggsave(plot_res, filename = paste0(set_name, plot_desc, slice_name, ".jpeg"), device = "jpeg", 
       path = paste(mulea_path, "/dev/plots_2023_08_02", sep = ""))
ggsave(plot_res, filename = paste0(set_name, plot_desc, slice_name, ".svg"), device = "svg", 
       path = paste(mulea_path, "/dev/plots_2023_08_02", sep = ""))
ggsave(plot_res, filename = paste0(set_name, plot_desc, slice_name, ".pdf"), device = "pdf", 
       path = paste(mulea_path, "/dev/plots_2023_08_02", sep = ""))


# OK - TPR wrapped by ratio.
plot_res <- sim_mult_tests_res_sum_slice %>% ggplot(
  aes(x=forcats::fct_relevel(method, 'p', 'bh', 'pt'), y=TPR, 
      fill=forcats::fct_relevel(method, 'p', 'bh', 'pt'))) + 
  geom_boxplot(alpha=0.5, fatten=2) +
  geom_violin(alpha=0.3) +
  facet_grid(~noise_ratio) +
  theme_bw() + 
  theme(legend.position="bottom", 
        axis.ticks.x=element_blank(), axis.text.x=element_blank()) +
  guides(fill=guide_legend(title='Method')) +
  scale_fill_discrete(palette = scales::hue_pal(l = 70), 
                      labels = c("Uncorrected p-value", "Benjamini-Hochberg", "eFDR")) +
  scale_x_discrete(position = "top") +
  xlab('Noise Ratio') + ylab('True Positive Rate') +
  ylim(c(0,1))
plot_desc <- "_TPR_ratio"
ggsave(plot_res, filename = paste0(set_name, plot_desc, slice_name, ".jpeg"), device = "jpeg", 
       path = paste(mulea_path, "/dev/plots_2023_08_02", sep = ""))
ggsave(plot_res, filename = paste0(set_name, plot_desc, slice_name, ".svg"), device = "svg", 
       path = paste(mulea_path, "/dev/plots_2023_08_02", sep = ""))
ggsave(plot_res, filename = paste0(set_name, plot_desc, slice_name, ".pdf"), device = "pdf", 
       path = paste(mulea_path, "/dev/plots_2023_08_02", sep = ""))

# FPR:
# OK - FPR wrapped by method.
plot_res <- sim_mult_tests_res_sum_slice %>% ggplot(aes(x=noise_ratio, y=FPR, fill=noise_ratio)) + 
  geom_boxplot(alpha=0.5, fatten=2) +
  geom_violin(alpha=0.3) +
  facet_grid(~forcats::fct_relevel(method, 'p', 'bh', 'pt'), labeller = as_labeller(
    c(`bh` = "Benjamini-Hochberg", `p` = "Uncorrected p-value", `pt` = "eFDR"))) + #facet_wrap
  theme_bw() + 
  theme(legend.position="right") +
  scale_fill_brewer(palette="Set1", guide="none") + 
  xlab('Noise Ratio') + ylab('False Positive Rate') +
  ylim(c(0,1))
plot_desc <- "_FPR_method"
ggsave(plot_res, filename = paste0(set_name, plot_desc, slice_name, ".jpeg"), device = "jpeg", 
       path = paste(mulea_path, "/dev/plots_2023_08_02", sep = ""))
ggsave(plot_res, filename = paste0(set_name, plot_desc, slice_name, ".svg"), device = "svg", 
       path = paste(mulea_path, "/dev/plots_2023_08_02", sep = ""))
ggsave(plot_res, filename = paste0(set_name, plot_desc, slice_name, ".pdf"), device = "pdf", 
       path = paste(mulea_path, "/dev/plots_2023_08_02", sep = ""))

# OK - FPR wrapped by ratio.
plot_res <- sim_mult_tests_res_sum_slice %>% ggplot(
  aes(x=forcats::fct_relevel(method, 'p', 'bh', 'pt'), y=FPR, 
      fill=forcats::fct_relevel(method, 'p', 'bh', 'pt'))) +
  geom_boxplot(alpha=0.5, fatten=2) +
  geom_violin(alpha=0.3) +
  facet_grid(~noise_ratio) +
  theme_bw() + 
  theme(legend.position="bottom", 
        axis.ticks.x=element_blank(), axis.text.x=element_blank()) +
  guides(fill=guide_legend(title='Method')) +
  scale_fill_discrete(palette = scales::hue_pal(l = 70), 
                      labels = c("Uncorrected p-value", "Benjamini-Hochberg", "eFDR")) +
  scale_x_discrete(position = "top") +
  xlab('Noise Ratio') + ylab('False Positive Rate') +
  ylim(c(0,1))
plot_desc <- "_FPR_ratio"
ggsave(plot_res, filename = paste0(set_name, plot_desc, slice_name, ".jpeg"), device = "jpeg", 
       path = paste(mulea_path, "/dev/plots_2023_08_02", sep = ""))
ggsave(plot_res, filename = paste0(set_name, plot_desc, slice_name, ".svg"), device = "svg", 
       path = paste(mulea_path, "/dev/plots_2023_08_02", sep = ""))
ggsave(plot_res, filename = paste0(set_name, plot_desc, slice_name, ".pdf"), device = "pdf", 
       path = paste(mulea_path, "/dev/plots_2023_08_02", sep = ""))

# FDR:
# OK - FDR wrapped by method.
plot_res <- sim_mult_tests_res_sum_slice %>% ggplot(aes(x=noise_ratio, y=FDR, fill=noise_ratio)) +
  geom_boxplot(alpha=0.5, fatten=2) +
  geom_violin(alpha=0.3) +
  facet_grid(~forcats::fct_relevel(method, 'p', 'bh', 'pt'), labeller = as_labeller(
    c(`bh` = "Benjamini-Hochberg", `p` = "Uncorrected p-value", `pt` = "eFDR"))) + #facet_wrap
  theme_bw() + 
  theme(legend.position="right") +
  scale_fill_brewer(palette="Set1", guide="none") + 
  xlab('Noise Ratio') + ylab('False Discovery Rate') +
  ylim(c(0,1))
plot_desc <- "_FDR_method"
ggsave(plot_res, filename = paste0(set_name, plot_desc, slice_name, ".jpeg"), device = "jpeg", 
       path = paste(mulea_path, "/dev/plots_2023_08_02", sep = ""))
ggsave(plot_res, filename = paste0(set_name, plot_desc, slice_name, ".svg"), device = "svg", 
       path = paste(mulea_path, "/dev/plots_2023_08_02", sep = ""))
ggsave(plot_res, filename = paste0(set_name, plot_desc, slice_name, ".pdf"), device = "pdf", 
       path = paste(mulea_path, "/dev/plots_2023_08_02", sep = ""))

# FDR wrapped by ratio.
plot_res <- sim_mult_tests_res_sum_slice %>% ggplot(
  aes(x=forcats::fct_relevel(method, 'p', 'bh', 'pt'), y=FDR, 
      fill=forcats::fct_relevel(method, 'p', 'bh', 'pt'))) +
  geom_boxplot(alpha=0.5, fatten=2) +
  geom_violin(alpha=0.3) +
  facet_grid(~noise_ratio) +
  theme_bw() + 
  theme(legend.position="bottom", 
        axis.ticks.x=element_blank(), axis.text.x=element_blank()) +
  guides(fill=guide_legend(title='Method')) +
  scale_fill_discrete(palette = scales::hue_pal(l = 70), 
                      labels = c("Uncorrected p-value", "Benjamini-Hochberg", "eFDR")) +
  scale_x_discrete(position = "top") +
  xlab('Noise Ratio') + ylab('False Discovery Rate') +
  ylim(c(0,1))
plot_desc <- "_FDR_ratio"
ggsave(plot_res, filename = paste0(set_name, plot_desc, slice_name, ".jpeg"), device = "jpeg", 
       path = paste(mulea_path, "/dev/plots_2023_08_02", sep = ""))
ggsave(plot_res, filename = paste0(set_name, plot_desc, slice_name, ".svg"), device = "svg", 
       path = paste(mulea_path, "/dev/plots_2023_08_02", sep = ""))
ggsave(plot_res, filename = paste0(set_name, plot_desc, slice_name, ".pdf"), device = "pdf", 
       path = paste(mulea_path, "/dev/plots_2023_08_02", sep = ""))

# Plots on ROC summary data.
# ROC curve plot.
plot_res <- sim_mult_tests_res_to_roc %>% 
  dplyr::filter(!is.na(noise_ratio) & (noise_ratio <= 0.2)) %>% ggplot(
  aes(FPR, TPR, 
      colour=forcats::fct_relevel(method, 'p_value', 'adjusted_p_value', 'eFDR'))) + 
  # geom_point() +
  geom_step(direction = 'vh') +
  facet_grid(~noise_ratio) +
  theme_bw() + 
  theme(legend.position="bottom") +
  # facet_wrap(~method) +
  guides(colour=guide_legend(title='Method')) +
  scale_colour_discrete(palette = scales::hue_pal(l = 70), 
                      labels = c("Uncorrected p-value", "Benjamini-Hochberg", "eFDR")) +
  geom_abline(color="black")
plot_desc <- "_ROC_ratio"
ggsave(plot_res, filename = paste0(set_name, plot_desc, slice_name, ".jpeg"), device = "jpeg", 
       path = paste(mulea_path, "/dev/plots_2023_08_02", sep = ""))
ggsave(plot_res, filename = paste0(set_name, plot_desc, slice_name, ".svg"), device = "svg", 
       path = paste(mulea_path, "/dev/plots_2023_08_02", sep = ""))
ggsave(plot_res, filename = paste0(set_name, plot_desc, slice_name, ".pdf"), device = "pdf", 
       path = paste(mulea_path, "/dev/plots_2023_08_02", sep = ""))

# ROC curve plot
plot_res <- sim_mult_tests_res_to_roc %>% dplyr::filter(is.na(noise_ratio)) %>% ggplot(
  aes(FPR, TPR, 
      colour=forcats::fct_relevel(method, 'p_value', 'adjusted_p_value', 'eFDR'))) + 
  # geom_point() +
  geom_step(direction = 'vh') +
  theme_bw() + 
  theme(legend.position="bottom") +
  # facet_wrap(~method) +
  guides(colour=guide_legend(title='Method')) +
  scale_colour_discrete(palette = scales::hue_pal(l = 70), 
                        labels = c("Uncorrected p-value", "Benjamini-Hochberg", "eFDR")) +
  geom_abline(color="black")
plot_desc <- "_ROC"
ggsave(plot_res, filename = paste0(set_name, plot_desc, slice_name, ".jpeg"), device = "jpeg", 
       path = paste(mulea_path, "/dev/plots_2023_08_02", sep = ""))
ggsave(plot_res, filename = paste0(set_name, plot_desc, slice_name, ".svg"), device = "svg", 
       path = paste(mulea_path, "/dev/plots_2023_08_02", sep = ""))
ggsave(plot_res, filename = paste0(set_name, plot_desc, slice_name, ".pdf"), device = "pdf", 
       path = paste(mulea_path, "/dev/plots_2023_08_02", sep = ""))


# ROC density curve plot
plot_res <- sim_mult_tests_res_to_roc %>% dplyr::filter(is.na(noise_ratio)) %>% ggplot(
  aes(FPR, TPR, 
      colour=forcats::fct_relevel(method, 'p_value', 'adjusted_p_value', 'eFDR'))) + 
  geom_point() +
  # geom_step() +
  theme_bw() + 
  theme(legend.position="bottom") +
  # facet_wrap(~method) +
  facet_wrap(~forcats::fct_relevel(method, 'p_value', 'adjusted_p_value', 'eFDR'), 
             labeller = as_labeller(
    c(`adjusted_p_value` = "Benjamini-Hochberg", 
      `p_value` = "Uncorrected p-value", `eFDR` = "eFDR"))) +
  guides(colour=guide_legend(title='Method')) +
  scale_colour_discrete(palette = scales::hue_pal(l = 70), 
                        labels = c("Uncorrected p-value", "Benjamini-Hochberg", "eFDR")) +
  geom_abline(color="black") +
  stat_density_2d(geom = "point", aes(size = after_stat(density)), n = 20, contour = FALSE)
plot_desc <- "_ROC_density"
ggsave(plot_res, filename = paste0(set_name, plot_desc, slice_name, ".jpeg"), device = "jpeg", 
       path = paste(mulea_path, "/dev/plots_2023_08_02", sep = ""))
ggsave(plot_res, filename = paste0(set_name, plot_desc, slice_name, ".svg"), device = "svg", 
       path = paste(mulea_path, "/dev/plots_2023_08_02", sep = ""))
ggsave(plot_res, filename = paste0(set_name, plot_desc, slice_name, ".pdf"), device = "pdf", 
       path = paste(mulea_path, "/dev/plots_2023_08_02", sep = ""))

# F1 score plot - all
plot_res <- sim_mult_tests_res_to_roc %>% dplyr::filter(is.na(noise_ratio)) %>% ggplot(
  aes(cut_off, f1_score, 
      colour=forcats::fct_relevel(method, 'p_value', 'adjusted_p_value', 'eFDR'))) + 
  # geom_point() +
  geom_step(direction = 'vh') +
  theme_bw() + 
  theme(legend.position="bottom") +
  # facet_wrap(~method) +
  guides(colour=guide_legend(title='Method')) +
  scale_colour_discrete(palette = scales::hue_pal(l = 70), 
                        labels = c("Uncorrected p-value", "Benjamini-Hochberg", "eFDR"))
  # geom_abline(color="black")
plot_desc <- "_F1_score"
ggsave(plot_res, filename = paste0(set_name, plot_desc, slice_name, ".jpeg"), device = "jpeg", 
       path = paste(mulea_path, "/dev/plots_2023_08_02", sep = ""))
ggsave(plot_res, filename = paste0(set_name, plot_desc, slice_name, ".svg"), device = "svg", 
       path = paste(mulea_path, "/dev/plots_2023_08_02", sep = ""))
ggsave(plot_res, filename = paste0(set_name, plot_desc, slice_name, ".pdf"), device = "pdf", 
       path = paste(mulea_path, "/dev/plots_2023_08_02", sep = ""))

# F1 score plot - per noise_ratio
plot_res <- sim_mult_tests_res_to_roc %>% dplyr::filter(!is.na(noise_ratio)) %>% ggplot(
  aes(cut_off, f1_score, 
      colour=forcats::fct_relevel(method, 'p_value', 'adjusted_p_value', 'eFDR'))) + 
  # geom_point() +
  geom_step(direction = 'vh') +
  facet_grid(~noise_ratio) +
  theme_bw() + 
  theme(legend.position="bottom") +
  # facet_wrap(~method) +
  guides(colour=guide_legend(title='Method')) +
  scale_colour_discrete(palette = scales::hue_pal(l = 70), 
                        labels = c("Uncorrected p-value", "Benjamini-Hochberg", "eFDR"))
plot_desc <- "_F1_score_ratio"
ggsave(plot_res, filename = paste0(set_name, plot_desc, slice_name, ".jpeg"), device = "jpeg", 
       path = paste(mulea_path, "/dev/plots_2023_08_02", sep = ""))
ggsave(plot_res, filename = paste0(set_name, plot_desc, slice_name, ".svg"), device = "svg", 
       path = paste(mulea_path, "/dev/plots_2023_08_02", sep = ""))
ggsave(plot_res, filename = paste0(set_name, plot_desc, slice_name, ".pdf"), device = "pdf", 
       path = paste(mulea_path, "/dev/plots_2023_08_02", sep = ""))

######################################################################################
#### END                                                                          ####
######################################################################################

calculate_auc <- function(res_to_roc) {
  res_to_roc
}


# EXTRA:






# IN PROGRESS:



plot_res <- sim_mult_tests_res_sum_slice %>% ggplot(aes(FPR, TPR, colour=method)) + 
  geom_point() +
  # geom_step() +
  facet_wrap(~method) +
  geom_abline(color="black") +
  scale_fill_brewer(palette="PuBu")


sim_mult_tests_res_sum_slice %>% ggplot(aes(FPR, TPR, colour=method)) + 
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


