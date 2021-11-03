library(MulEA)

# Read and filter inputs.
gmtFilePath <- paste(find.package("MulEA"), 
                     "/tests/inputs/KEGG_example_c_elegans_Leila.gmt", sep = "")
input_gmt <- MulEA::readGmtFileAsDataFrame(gmtFilePath)
input_gmt_filtered <- MulEA::filterOntology(input_gmt = input_gmt)
filteredGmtFilePath <- paste(find.package("MulEA"), 
                             "/tests/outputs/KEGG_filtered.gmt", sep = "")
MulEA::saveDataFrameAsGmtFile(modelDF = input_gmt_filtered, gmtFilePath = filteredGmtFilePath)
input_gmt_filtered <- MulEA::readGmtFileAsDataFrame(gmtFilePath = filteredGmtFilePath)

# Set seed.
set.seed(seed = 1234)

# Construct search space.
construct_comparison_space <- function(nubmer_of_over_repr_terms_range=5:10, 
                                       nubmer_of_under_repr_terms_range=5:10,
                                       sample_ratio_range=seq(0.1, 0.3, by=0.1),
                                       over_repr_ratio_range=seq(0.1, 0.5, by=0.2), 
                                       under_repr_ratio_range=seq(0.1, 0.2, by=0.1)) {
  comparison_space <- expand.grid(
    no_over_repr_terms=nubmer_of_over_repr_terms_range, 
    no_under_repr_terms=nubmer_of_under_repr_terms_range,
    sample_ratio=sample_ratio_range, 
    over_repr_ratio=over_repr_ratio_range,
    under_repr_ratio=under_repr_ratio_range,
    score=0)
  
  # TODO : Validation of space
  comparison_space <- comparison_space[comparison_space$sample_ratio >= comparison_space$under_repr_ratio,]
  comparison_space <- comparison_space[comparison_space$sample_ratio+comparison_space$over_repr_ratio <= 1,]
  
  return(comparison_space)
}

comparison_space <- construct_comparison_space()

# Call method on search space.
compare_by_treshold <- function(comparison_space, input_gmt_filtered, treshold=0.05) {
  # DEB : i <- 2
  comparison_res <- data.frame(matrix(ncol = 12, nrow = 0))
  colnames(comparison_res) <- c('pt_PCER', 'bh_PCER', 'pv_PCER', 
                                'pt_FDR', 'bh_FDR', 'pv_FDR', 
                                'pt_TPR', 'bh_TPR', 'pv_TPR', 
                                'pt_FPR', 'bh_FPR', 'pv_FPR')
  for (i in 1:length(comparison_space$no_over_repr_terms)) {
    comparison_case <- comparison_space[i,]
    print(i)
    
    input_generated <- MulEA::generateInputData(
      input_gmt = input_gmt_filtered, 
      sample_ratio=comparison_case$sample_ratio,
      over_repr_ratio = comparison_case$over_repr_ratio, 
      under_repr_ratio = comparison_case$under_repr_ratio,
      number_of_over_representation_groups = comparison_case$no_over_repr_terms,
      number_of_under_representation_groups = comparison_case$no_under_repr_terms)
    
    # Perform ORA test.
    input_select <- input_generated$input_select
    number_of_steps <- 10000
    
    # TODO : benchmark time.
    mulea_ora_model <- MulEA::ORA(
      gmt = input_gmt_filtered, testData = input_select, adjustMethod = "PT",
      numberOfPermutations = number_of_steps, nthreads = 16)
    
    mulea_ora_results <- NULL
    mulea_ora_results <- MulEA::runTest(mulea_ora_model)
    
    # Compare results.
    over_repr_in <- input_generated$gmt_for_generator[
      input_generated$gmt_for_generator$sample_label == "over",]$ontologyId
    under_repr_in <- input_generated$gmt_for_generator[
      input_generated$gmt_for_generator$sample_label == "under",]$ontologyId
    
    pt_res <- mulea_ora_results[
      mulea_ora_results$adjustedPValueEmpirical < treshold,]$ontologyId
    bh_res <- mulea_ora_results[
      mulea_ora_results$adjustedPValue < treshold,]$ontologyId
    pv_res <- mulea_ora_results[
      mulea_ora_results$pValue < treshold,]$ontologyId
    
    
    # Per comparison error rate (PCER): the expected value of the number
    # of Type I errors over the number of hypotheses,
    # PCER = E(V)/m
    m <- nrow(mulea_ora_results)
    
    pt_false_positive <- pt_res %in% over_repr_in
    pt_false_positive_no <- length(pt_false_positive) - sum(pt_false_positive)
    pt_PCER <- pt_false_positive_no/m
    
    bh_false_positive <- bh_res %in% over_repr_in
    bh_false_positive_no <- length(bh_false_positive) - sum(bh_false_positive)
    bh_PCER <- bh_false_positive_no/m
    
    pv_false_positive <- pv_res %in% over_repr_in
    pv_false_positive_no <- length(pv_false_positive) - sum(pv_false_positive)
    pv_PCER <- pv_false_positive_no/m
    
    # False discovery rate (FDR) is the expected proportion of Type I errors
    # among the rejected hypotheses
    # FDR = E(V/R | R>0)P(R>0)
    pt_false_positive <- pt_res %in% over_repr_in
    pt_false_positive_no <- length(pt_false_positive) - sum(pt_false_positive)
    pt_FDR <- pt_false_positive_no/length(pt_false_positive)
    
    bh_false_positive <- bh_res %in% over_repr_in
    bh_false_positive_no <- length(bh_false_positive) - sum(bh_false_positive)
    bh_FDR <- bh_false_positive_no/length(bh_false_positive)
    
    pv_false_positive <- pv_res %in% over_repr_in
    pv_false_positive_no <- length(pv_false_positive) - sum(pv_false_positive)
    pv_FDR <- pv_false_positive_no/length(pv_false_positive)
    
    # True positive rate (TPR) 
    # TPR = TP/(TP+FN)
    pt_ture_positive <- over_repr_in %in% pt_res
    pt_ture_positive_no <- sum(pt_ture_positive)
    pt_TPR <- pt_ture_positive_no/length(over_repr_in)
    
    bh_ture_positive <- over_repr_in %in% bh_res
    bh_ture_positive_no <- sum(bh_ture_positive)
    bh_TPR <- bh_ture_positive_no/length(over_repr_in)
    
    pv_ture_positive <- over_repr_in %in% pv_res
    pv_ture_positive_no <- sum(pv_ture_positive)
    pv_TPR <- pv_ture_positive_no/length(over_repr_in)
    
    # False positive rate (FPR)
    # FPR = FP/(FP+TN)
    
    pt_false_positive <- pt_res %in% over_repr_in
    pt_false_positive_no <- length(pt_false_positive) - sum(pt_false_positive)
    pt_FPR <- pt_false_positive_no/(nrow(mulea_ora_results)-length(over_repr_in))
    
    bh_false_positive <- bh_res %in% over_repr_in
    bh_false_positive_no <- length(bh_false_positive) - sum(bh_false_positive)
    bh_FPR <- bh_false_positive_no/(nrow(mulea_ora_results)-length(over_repr_in))
    
    pv_false_positive <- pv_res %in% over_repr_in
    pv_false_positive_no <- length(pv_false_positive) - sum(pv_false_positive)
    pv_FPR <- pv_false_positive_no/(nrow(mulea_ora_results)-length(over_repr_in))
    
    # over_repr_pt_res <- over_repr_in %in% pt_res
    # over_repr_bh_res <- over_repr_in %in% bh_res
    # over_repr_pv_res <- over_repr_in %in% pv_res
    # pt_res_score <- sum(over_repr_pt_res) / length(over_repr_in)
    # bh_res_score <- sum(over_repr_bh_res) / length(over_repr_in)
    # pv_res_score <- sum(over_repr_pv_res) / length(over_repr_in)
    # pt_res_score_ratio <- pt_res_score / pt_res
    # bh_res_score_ratio <- bh_res_score / bh_res
    # pv_res_score_ratio <- pv_res_score_ratio / pv_res
    comparison_res[i, ] = data.frame(pt_PCER, bh_PCER, pv_PCER, 
                                     pt_FDR, bh_FDR, pv_FDR, 
                                     pt_TPR, bh_TPR, pv_TPR,
                                     pt_FPR, bh_FPR, pv_FPR)
    
  }
  return(comparison_res)
}

nrow(comparison_space)

comp_spac_small_enr <- comparison_space[1:10,]
compare_res_small <- compare_by_treshold(comparison_space = comp_spac_small_enr, 
                                   input_gmt_filtered = input_gmt_filtered)
boxplot(compare_res_small)


comp_spac_small_big <- comparison_space[506:540,]
compare_res_big <- compare_by_treshold(comparison_space = comp_spac_small_big, 
                                   input_gmt_filtered = input_gmt_filtered)
boxplot(compare_res_big)
plot(x = compare_res_big$bh_FPR, y= compare_res_big$bh_TPR, main="FPR to TPR",
     xlab="FPR ", ylab="TPR ", xlim = c(0,1), ylim = c(0,1))
points(x = compare_res_big$pt_FPR,
       y = compare_res_big$pt_TPR,
       col = 'red')
points(x = compare_res_big$pv_FPR,
       y = compare_res_big$pv_TPR,
       col = 'blue')

# DELETE - for debug only
input_gmt=input_gmt_filtered
sample_ratio=0.3
over_repr_ratio=0.5 
under_repr_ratio=0.2
number_of_over_representation_groups=10 
number_of_under_representation_groups=10 
turn_on_log=FALSE
treshold=0.05






