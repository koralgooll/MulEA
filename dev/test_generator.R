
gmtFilePath <- paste(
  find.package("MulEA"), 
  "/tests/inputs/Bands_Drosophila_melanogaster_Flybase_Marton.gmt", sep = "")
gmtFilePath <- paste(
  find.package("MulEA"), 
  "/tests/inputs/KEGG_example_c_elegans_Leila.gmt", sep = "")


input_gmt <- MulEA::readGmtFileAsDataFrame(gmtFilePath)

filter_ontology <- function(input_gmt, min=NULL, max=NULL) {
  if (is.null(min)) {
    
  }
  filtered_input_gmt <- plyr::ddply(.data = input_gmt, .variables = c("ontologyId"), .fun = function(df_row) {
    if (length(df_row$listOfValues[[1]]) > min) {
      df_row
    } else {
      df_row[-1,]
    }
  })
  filtered_input_gmt <- plyr::ddply(.data = filtered_input_gmt, .variables = c("ontologyId"), .fun = function(df_row) {
    if (length(df_row$listOfValues[[1]]) < max) {
      df_row
    } else {
      df_row[-1,]
    }
  })
  filtered_input_gmt
}


generate_input_data <- function(input_gmt, sample_ratio=0.5,
                                group_under_over_representation_ratio=0.9,
                                number_of_over_representation_groups=1, 
                                number_of_under_representation_groups=1) {
  input_select <- plyr::llply(.data = input_gmt$listOfValues, .fun = function(term) {
    size <- length(term)
    size_of_sample <- floor(size * sample_ratio)
    input_select_indicator <- sample(1:size, size_of_sample, replace=FALSE)
    input_select_term <- term[input_select_indicator]
  })
  
  input_select <- unique(unlist(input_select))
  input_select
  
  go_size <- length(input_gmt$listOfValues)
  size_of_over_under_repr <- number_of_over_representation_groups + number_of_under_representation_groups
  go_change_repr <- sample(1:go_size, size_of_over_under_repr, replace=FALSE)
  
  go_change_repr_over <- go_change_repr[1:number_of_over_representation_groups]
  go_change_repr_under <- go_change_repr[-1:-number_of_over_representation_groups]
  
  print('***')
  
  
  over_repr_log <- sapply(as.character(go_change_repr_over), function(x) NULL)
  for (term_id in go_change_repr_over) {
    term <- input_gmt$listOfValues[[term_id]]
    size <- length(term)
    size_of_sample <- floor(size * group_under_over_representation_ratio)
    input_select_indicator <- sample(1:size, size_of_sample, replace=FALSE)
    input_select_term <- term[input_select_indicator]
    input_select <- c(input_select, input_select_term)
    over_repr_log[[as.character(term_id)]] <- paste(term_id, ' [', size, '] ', ' -> ', ' [', length(input_select_term), ']', sep = '')
  }
  

  under_repr_log <- sapply(as.character(go_change_repr_under), function(x) NULL)
  for (term_id in go_change_repr_under) {
    term <- input_gmt$listOfValues[[term_id]]
    size <- length(term)
    size_of_sample <- floor(size * group_under_over_representation_ratio)
    input_select_indicator <- sample(1:size, size_of_sample, replace=FALSE)
    input_select_term <- term[input_select_indicator]
    input_select <- setdiff(input_select, input_select_term)
    under_repr_log[[as.character(term_id)]] <- paste(term_id, ' [', size, '] ', ' -> ', ' [', length(input_select_term), ']', sep = '')
  }
  
  print("Over representation terms:")
  for (term_id in names(over_repr_log)) {
    real_term_in_sample <- intersect(input_gmt$listOfValues[[as.integer(term_id)]], 
                                     input_select)
    print(paste(over_repr_log[[term_id]], ' {', length(real_term_in_sample), '}', sep = ''))
  }
  
  print("Under representation terms:")
  for (term_id in names(under_repr_log)) {
    real_term_in_sample <- intersect(input_gmt$listOfValues[[as.integer(term_id)]], 
                                     input_select)
    print(paste(under_repr_log[[term_id]], ' {', length(real_term_in_sample), '}', sep = ''))
  }

  input_select
}


# Start of test of method.
library(MulEA)

number_of_steps <- 10

hist_data <- c()


terms_sizes <- plyr::laply(.data = input_gmt$listOfValues, .fun = function(term) {
  length(term)
})
term_size_dist_q <- quantile(terms_sizes, probs = seq(0, 1, 0.1), type = 2, na.rm = FALSE)

min_go_term_size = term_size_dist_q['20%']
max_go_term_size = term_size_dist_q['80%']

no_over_repr_terms=3 
no_under_repr_terms=2

for (i in 1:10) {
  input_gmt_filtered <- filter_ontology(input_gmt = input_gmt, 
                                        min=min_go_term_size, 
                                        max=max_go_term_size)
  
  input_select <- generate_input_data(
    input_gmt = input_gmt_filtered, 
    number_of_over_representation_groups = no_over_repr_terms,
    number_of_under_representation_groups = no_under_repr_terms)
  hist_data <- c(hist_data, length(input_select))
  
  mulea_ora_model <- MulEA::ORA(
    gmt = input_gmt_filtered, testData = input_select, adjustMethod = "PT",
    numberOfPermutations = number_of_steps)
  mulea_ora_results <- MulEA::runTest(mulea_ora_model)
}

hist(hist_data)


# Other stats 
terms_sizes_log <- log(terms_sizes)

hist(terms_sizes, breaks=100)
hist(terms_sizes_log, breaks=100)

quantile(terms_sizes, probs = seq(0, 1, 0.1), type = 2, na.rm = FALSE)
quantile(terms_sizes_log, probs = seq(0, 1, 0.1), type = 2, na.rm = FALSE)
exp(1.386294)
exp(4.394144)

mean(terms_sizes)
mean(terms_sizes_log)
exp(3.056871)
median(terms_sizes)
median(terms_sizes_log)
sd(terms_sizes_log)


# Setup params


# Perform ORA
mulea_ora_model <- MulEA::ORA(
  gmt = input_gmt, testData = input_select, adjustMethod = "PT",
  numberOfPermutations = number_of_steps)
mulea_ora_results <- MulEA::runTest(mulea_ora_model)


mulea_ora_reshaped_results <- MulEA::reshapeResults(
  mulea_model=mulea_ora_model, 
  mulea_model_resuts=mulea_ora_results, 
  category_stat_column_name='adjustedPValueEmpirical')

warnings()

# Plot results
MulEA::plotBarplot(mulea_relaxed_resuts = mulea_ora_reshaped_results, 
                   statistics_value_colname = "adjustedPValueEmpirical",
                   statistics_value_cutoff=0.05)
MulEA::plotHeatmap(mulea_relaxed_resuts=mulea_ora_reshaped_results, 
                   statistics_value_colname = "adjustedPValueEmpirical",
                   statistics_value_cutoff=0.05)
MulEA::plotGraph(mulea_relaxed_resuts=mulea_ora_reshaped_results,
                 statistics_value_colname = "adjustedPValueEmpirical",
                 statistics_value_cutoff = 0.05)
