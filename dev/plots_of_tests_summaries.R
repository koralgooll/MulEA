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
sim_mult_tests_res_to_roc <- MulEA::getSummaryToRoc(tests_res = sim_mult_tests_res)

# TPR:
# TPR wrapped by method.
plot_res <- sim_mult_tests_res_sum %>% ggplot(aes(x=noise_ratio, y=TPR, fill=noise_ratio)) + 
  geom_boxplot(alpha=0.5, fatten=2) +
  geom_violin(alpha=0.3) +
  facet_grid(~method) + #facet_wrap
  theme(legend.position="right") +
  scale_fill_brewer(palette="PuBuGn") + 
  xlab('Noise Ratio') + ylab('True Positive Rate') +
  ylim(c(0,1))
ggsave(plot_res, filename = paste0(set_name, "_TPR_method", ".jpeg"), device = "jpeg", 
         path = "D:\\projects\\science\\MulEA\\dev\\new_plots")

# TPR wrapped by ratio.
plot_res <- sim_mult_tests_res_sum %>% ggplot(aes(x=method, y=TPR, fill=method)) + 
  geom_boxplot(alpha=0.5, fatten=2) +
  geom_violin(alpha=0.3) +
  facet_grid(~noise_ratio) +
  theme(legend.position="right") +
  scale_fill_brewer(palette="PuBuGn") + 
  xlab('Noise Ratio') + ylab('True Positive Rate') +
  ylim(c(0,1))
ggsave(plot_res, filename = paste0(set_name, "_TPR_ratio", ".jpeg"), device = "jpeg", 
       path = "D:\\projects\\science\\MulEA\\dev\\new_plots")


# FPR:
# FPR wrapped by method.
plot_res <- sim_mult_tests_res_sum %>% ggplot(aes(x=noise_ratio, y=FPR, fill=noise_ratio)) + 
  geom_boxplot(alpha=0.5, fatten=2) +
  geom_violin(alpha=0.3) +
  facet_grid(~method) + #facet_wrap
  theme(legend.position="right") +
  scale_fill_brewer(palette="PuBuGn") + 
  xlab('Noise Ratio') + ylab('False Positive Rate') +
  ylim(c(0,1))
ggsave(plot_res, filename = paste0(set_name, "_FPR_method", ".jpeg"), device = "jpeg", 
       path = "D:\\projects\\science\\MulEA\\dev\\new_plots")

# FPR wrapped by ratio.
plot_res <- sim_mult_tests_res_sum %>% ggplot(aes(x=method, y=FPR, fill=method)) +
  geom_boxplot(alpha=0.5, fatten=2) +
  geom_violin(alpha=0.3) +
  facet_grid(~noise_ratio) +
  theme(legend.position="right") +
  scale_fill_brewer(palette="PuBuGn") + 
  xlab('Noise Ratio') + ylab('False Positive Rate')
ggsave(plot_res, filename = paste0(set_name, "_FPR_ratio", ".jpeg"), device = "jpeg", 
       path = "D:\\projects\\science\\MulEA\\dev\\new_plots")


# FDR:
# FDR wrapped by method.
plot_res <- sim_mult_tests_res_sum %>% ggplot(aes(x=noise_ratio, y=FDR, fill=noise_ratio)) +
  geom_boxplot(alpha=0.5, fatten=2) +
  geom_violin(alpha=0.3) +
  facet_grid(~method) + #facet_wrap
  theme(legend.position="right") +
  scale_fill_brewer(palette="PuBuGn") + 
  xlab('Noise Ratio') + ylab('False Discovery Rate')
ggsave(plot_res, filename = paste0(set_name, "_FDR_method", ".jpeg"), device = "jpeg", 
       path = "D:\\projects\\science\\MulEA\\dev\\new_plots")

# FDR wrapped by ratio.
plot_res <- sim_mult_tests_res_sum %>% ggplot(aes(x=method, y=FDR, fill=method)) +
  geom_boxplot(alpha=0.5, fatten=2) +
  geom_violin(alpha=0.3) +
  facet_grid(~noise_ratio) +
  theme(legend.position="right") +
  scale_fill_brewer(palette="PuBuGn") + 
  xlab('Noise Ratio') + ylab('False Discovery Rate')
ggsave(plot_res, filename = paste0(set_name, "_FDR_ratio", ".jpeg"), device = "jpeg", 
       path = "D:\\projects\\science\\MulEA\\dev\\new_plots")


# ROC curve plot
plot_res <- sim_mult_tests_res_to_roc %>% ggplot(aes(FPR, TPR, colour=method)) + 
  # geom_point() +
  geom_step() +
  # facet_wrap(~method) +
  geom_abline(color="black") +
  scale_fill_brewer(palette="PuBu")
ggsave(plot_res, filename = paste0(set_name, "_ROC", ".jpeg"), device = "jpeg", 
       path = "D:\\projects\\science\\MulEA\\dev\\new_plots")


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


