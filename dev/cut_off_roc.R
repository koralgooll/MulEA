
data_to_roc <- MulEA::get(tests_res = sim_mult_tests_res)


c("noise", "over", "over") == "over"


#define object to plot
roc_pvalue <- pROC::roc(data_to_roc$sample_label, data_to_roc$pValue, 
                  smoothed = TRUE)
roc_pvalue_adj <- pROC::roc(data_to_roc$sample_label, data_to_roc$adjustedPValue,
                      smoothed = TRUE)
roc_perm_pvalue <- pROC::roc(data_to_roc$sample_label, data_to_roc$adjustedPValueEmpirical, 
                       smoothed = TRUE)


library(pROC)
rocobj1 <- roc(df$actualoutcome1, data$prediction1)
rocobj2 <- roc(df$actualoutcome1, data$prediction2)
ggroc(list("pvalue" = roc_pvalue, "pvalue_adjust" = roc_pvalue_adj, 
           "pr_empirical" = roc_perm_pvalue))


library(PRROC)

PRROC_obj <- PRROC::roc.curve(scores.class0 = data_to_roc$pValue, 
                       weights.class0 = data_to_roc$sample_label,
                       curve=TRUE)
plot(PRROC_obj)


library(plotROC)
rocplot <- ggplot(data_to_roc, aes(m = pValue, d = sample_label)) + 
  geom_roc(n.cuts=20,labels=FALSE)
rocplot + style_roc(theme = theme_grey) + geom_rocci(fill="pink") 


basicplot <- ggplot(data_to_roc, aes(d = sample_label, m = pValue)) + 
  geom_roc(n.cuts = 50, labels = FALSE) + 
  style_roc(theme = theme_grey, xlab = "FPR") + 
  geom_rocci(ci.at = quantile(data_to_roc$pValue, c(.1, .4, .5, .6, .9)))
basicplot
basicplot2 <- ggplot(data_to_roc, aes(d = sample_label, m = adjustedPValue)) + geom_roc()
basicplot2


longtest <- melt_roc(data_to_roc, "sample_label", 
                     c("pValue", "adjustedPValue", "adjustedPValueEmpirical"))
ggplot(longtest, aes(d = D, m = M, color = name)) + geom_roc() + style_roc()





sumary_res_tmp <- data.frame(
  'test_no' = i, 
  'TP' = I(list(TP)), 'TP_size' = TP_size, 
  'FP' = I(list(FP)), 'FP_size' = FP_size,
  'FN' = I(list(FN)), 'FN_size' = FN_size,
  'TN' = I(list(TN)), 'TN_size' = TN_size,
  'over_repr_terms' = I(list(over_repr_terms)))


roc_stats <- tibble(
  TP_val = numeric(),
  TN_val = numeric(),
  FP_val = numeric(),
  FN_val = numeric(),
  TPR = numeric(),
  FPR = numeric(),
  sum_test = numeric(),
  cut_off = numeric(),
  method = character()
)

for (method_name in c('pValue', 'adjustedPValue', 'adjustedPValueEmpirical')) {
  for (cut_off in seq(0, 1, 0.001)) {
    sim_mult_tests_res_to_roc_summary <- sim_mult_tests_res_to_roc %>% 
      mutate(., PP=!!as.name(method_name)<=cut_off) %>%
      mutate(., TP=(PP == TRUE & sample_label=='over'), 
             TN=(PP == FALSE & sample_label!='over'),
             FP=(PP == TRUE & sample_label!='over'),
             FN=(PP == FALSE & sample_label=='over'))
    
    sim_sum <- sim_mult_tests_res_to_roc_summary %>% summarise(
      TP_val = sum(TP), TN_val = sum(TN), FP_val = sum(FP), FN_val = sum(FN))
    
    sim_sum_roc <- sim_sum %>% mutate(
      TPR = TP_val/(TP_val+FN_val),
      FPR = FP_val/(FP_val+TN_val),
      sum_test = TP_val+TN_val+FP_val+FN_val,
      cut_off = cut_off,
      method=method_name)
    
    roc_stats <- roc_stats %>% add_row(sim_sum_roc)
  }
}









cut_off <- 1
sim_mult_tests_res_to_roc_summary <- sim_mult_tests_res_to_roc %>% 
  mutate(., PP=pValue<=cut_off) %>%
  mutate(., TP=(PP == TRUE & sample_label=='over'), 
            TN=(PP == FALSE & sample_label!='over'),
            FP=(PP == TRUE & sample_label!='over'),
            FN=(PP == FALSE & sample_label=='over'))

sim_sum <- sim_mult_tests_res_to_roc_summary %>% summarise(
  TP_val = sum(TP), TN_val = sum(TN), FP_val = sum(FP), FN_val = sum(FN))

sim_sum_roc <- sim_sum %>% mutate(
  TPR = TP_val/(TP_val+FN_val),
  FPR = FP_val/(FP_val+TN_val),
  sum_test = TP_val+TN_val+FP_val+FN_val)






