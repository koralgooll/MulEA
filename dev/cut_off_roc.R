
# Best noise:

getMultipleTestsSummary <- function(
  tests_res, 
  cut_off = 0.05, 
  comparison_col_name = 'pValue') {
  
  # Summarize results.
  print("Mulea sumary time:")
  tic()
  
  sumary_res <- data.frame(matrix(ncol = 9, nrow = 0))
  colnames(sumary_res) <- c('test_no', 
                            'TP', 'TP_size', 
                            'FP', 'FP_size',
                            'FN', 'FN_size',
                            'TN', 'TN_size'
  )
  for (i in 1:number_of_tests) {
    # Actual condition
    # Total population = P + N
    total_population <- tests_res[[i]]$test_data$ontologyId
    total_population_size <- length(total_population)
    # Positive (P)
    P <- tests_res[[i]]$test_data[tests_res[[i]]$test_data$sample_label == 'over',]$ontologyId
    P_size <- length(P)
    # Negative (N)
    N <- setdiff(total_population, P)
    N_size <- length(N)
    if (P_size + N_size != total_population_size) {
      warning("Not OK size of Actual in contingency table")
    }
    
    # Predicted condition
    # Predicted Positive (PP)
    PP <- tests_res[[i]]$mulea_res[tests_res[[i]]$mulea_res[, comparison_col_name] <= cut_off, ]$ontologyId
    PP_size <- length(PP)
    # Predicted Negative (PN)
    PN <- tests_res[[i]]$mulea_res[tests_res[[i]]$mulea_res[, comparison_col_name] > cut_off, ]$ontologyId
    PN_size <- length(PN)
    if (PP_size + PN_size != total_population_size) {
      warning("Not OK size of Predicted in contingency table")
    }
    
    # True positive (TP) : hit
    TP <- intersect(P, PP)
    TP_size <- length(TP)
    # False positive (FP) : type I error, false alarm, overestimation
    FP <- setdiff(PP, P)
    FP_size <- length(FP)
    # False negative (FN) : type II error, miss, underestimation
    FN <- setdiff(P, PP)
    FN_size <- length(FN)
    # True negative (TN) : correct rejection
    TN <- setdiff(total_population, union(P, PP))
    TN_size <- length(TN)
    
    if (TP_size + FP_size + FN_size + TN_size != total_population_size) {
      warning("Not OK size of total  contingency table")
    }
    
    sumary_res[i, ] <- data.frame(
      'test_no' = i, 
      'TP' = I(list(TP)), 'TP_size' = TP_size, 
      'FP' = I(list(FP)), 'FP_size' = FP_size,
      'FN' = I(list(FN)), 'FN_size' = FN_size,
      'TN' = I(list(TN)), 'TN_size' = TN_size)
  }
  
  sumary_res <- tibble(sumary_res) %>% 
    mutate(FPR=FP_size/(FP_size+TN_size)) %>% 
    mutate(TPR=TP_size/(TP_size+FN_size))
  
  toc()
  
  return(sumary_res)
}

# Read tests for different ratios
mult_tests_01 <-read_rds("tests_res-5-05-01-1000.rds")
mult_tests_02 <-read_rds("tests_res-5-05-02-1000.rds")
mult_tests_03 <-read_rds("tests_res-5-05-03-1000.rds")
mult_tests_04 <-read_rds("tests_res-5-05-04-1000.rds")
mult_tests_05 <-read_rds("tests_res-5-05-05-1000.rds")

mult_tests_01_sum <- getMultipleTestsSummary(tests_res = mult_tests_01)
mult_tests_02_sum <- getMultipleTestsSummary(tests_res = mult_tests_02)
mult_tests_03_sum <- getMultipleTestsSummary(tests_res = mult_tests_03)
mult_tests_04_sum <- getMultipleTestsSummary(tests_res = mult_tests_04)
mult_tests_05_sum <- getMultipleTestsSummary(tests_res = mult_tests_05)

comp_mult_tests <- data.frame(
  x = c(mult_tests_01_sum$TPR, mult_tests_02_sum$TPR, 
        mult_tests_03_sum$TPR, mult_tests_04_sum$TPR, mult_tests_05_sum$TPR), 
  g = gl(5, nrow(mult_tests_01_sum), label=c('0.1','0.2','0.3', '0.4', '0.5')))

# Bottom Left
comp_mult_tests %>% ggplot(aes(x=g, y=x, fill=g)) + 
  geom_boxplot(alpha=0.3) +
  theme(legend.position="right") +
  scale_fill_brewer(palette="BuPu") + 
  xlab('Noise Ratio') + ylab('True Positive Ratio')



mult_tests_01_sum <- getMultipleTestsSummary(tests_res = mult_tests_01)
mult_tests_02_sum <- getMultipleTestsSummary(tests_res = mult_tests_02)
mult_tests_03_sum <- getMultipleTestsSummary(tests_res = mult_tests_03)
mult_tests_04_sum <- getMultipleTestsSummary(tests_res = mult_tests_04)
mult_tests_05_sum <- getMultipleTestsSummary(tests_res = mult_tests_05)

comparison_col_name = 'adjustedPValueEmpirical'
comparison_col_name = 'adjustedPValue'

names(mulea_ora_results)


# REMOVE
comp_mult_tests %>%
  ggplot( aes(x=g, y=x, fill=x)) +
  geom_boxplot() +
  scale_fill_viridis(alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("A boxplot with jitter") +
  xlab("")


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

