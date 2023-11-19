#' Read GMT File
#' 
#' Read model data from a Gene Matrix Transposed (GMT) file and parse into a
#' data frame.
#'
#' @param file character; a path which points to a GMT file.
#' @return Returns data frame with the model from a specific location.
#' @export
#'

read_gmt <- function(file) {
  fileConnection <- file(file)
  tryCatchRes <- tryCatch(
    lines <- readLines(fileConnection),
    warning = function(w) {
    }
  )
  close(fileConnection)
  
  lines <- lines[!grepl('^#+', lines, fixed = FALSE)]
  lines <- lines["" != lines]
  
  gmtAsDF <-
    plyr::adply(
      .data = lines,
      .margins = 1,
      .fun = function(line) {
        fields <- strsplit(line, split = "\t")[[1]]
        category <- fields[1]
        if (startsWith(fields[2], "\"") &&
            endsWith(fields[2], "\"")) {
          description <- fields[2]
        } else {
          description <- paste("\"", fields[2], "\"", sep = "")
        }
        listOfValues <- fields[3:length(fields)]
        data.frame(
          'ontologyId' = category,
          'ontologyName' = description,
          'listOfValues' = I(list(listOfValues)),
          stringsAsFactors = FALSE
        )
      }
    )
  gmtAsDF[c("ontologyId", "ontologyName", "listOfValues")]
}


#' Write GMT file
#' 
#' Write the model data to a Gene Matrix Transposed (GMT) file.
#'
#' @param gmt A data frame which contains the data, imported from a GMT file.
#' @param file character; a path naming a file.
#' @return Returns the model as a .gmt file at a specific location.
#' @export
#' @examples
#' modelDfFromFile <- read_gmt(
#'   file = system.file(package="MulEA", "extdata", "model.gmt"))
#' write_gmt(gmt = modelDfFromFile,
#'           file = paste(system.file(package="MulEA", "extdata"),
#'           "fromDb.gmt", sep = "/"))
write_gmt <- function(gmt, file) {
  vectorOfModel <-
    plyr::daply(
      .data = gmt,
      .variables = c("ontologyId"),
      .fun = function(dataFrameRow) {
        collapsedListOfValues <-
          paste(dataFrameRow[, 3][[1]], collapse = "\t")
        paste(dataFrameRow[1], dataFrameRow[2], collapsedListOfValues, sep = "\t")
      }
    )
  fileConnection <- file(file)
  writeLines(vectorOfModel,
             con = fileConnection,
             sep = "\n",
             useBytes = FALSE)
  close(fileConnection)
}

#' @description
#' Filter Ontology
#' 
#' Filter ontology to only contain terms between given min. and max. sizes.
#'
#' @param gmt A data frame which contains the data, imported from a GMT file.
#' @param min_nr_of_elements minimum size of term. Default 20 percent from
#' quantile on term size distribution.
#' @param max_nr_of_elements maximum size of term. Default 80 percent from
#' quantile on term size distribution.
#' @return Return data frame with model from specific location.
#' @export
filter_ontology <- function(gmt,
                           min_nr_of_elements = 5,
                           max_nr_of_elements = 350) {
  # TODO : Add quantile parameters as separate! Nothing do is default.
  if (is.null(min_nr_of_elements)) {
    terms_sizes <-
      plyr::laply(
        .data = gmt$listOfValues,
        .fun = function(term) {
          length(term)
        }
      )
    term_size_dist_q <-
      stats::quantile(
        terms_sizes,
        probs = seq(0, 1, 0.1),
        type = 2,
        na.rm = FALSE
      )
    
    min_nr_of_elements = term_size_dist_q['20%']
  }
  
  if (is.null(max_nr_of_elements)) {
    terms_sizes <-
      plyr::laply(
        .data = gmt$listOfValues,
        .fun = function(term) {
          length(term)
        }
      )
    term_size_dist_q <-
      stats::quantile(
        terms_sizes,
        probs = seq(0, 1, 0.1),
        type = 2,
        na.rm = FALSE
      )
    max_nr_of_elements = term_size_dist_q['80%']
  }
  
  filtered_input_gmt <-
    plyr::ddply(
      .data = gmt,
      .variables = c("ontologyId"),
      .fun = function(df_row) {
        if (length(df_row$listOfValues[[1]]) > min_nr_of_elements) {
          df_row
        } else {
          df_row[-1, ]
        }
      }
    )
  filtered_input_gmt <-
    plyr::ddply(
      .data = filtered_input_gmt,
      .variables = c("ontologyId"),
      .fun = function(df_row) {
        if (length(df_row$listOfValues[[1]]) < max_nr_of_elements) {
          df_row
        } else {
          df_row[-1, ]
        }
      }
    )
  filtered_input_gmt
}


# PUBLIC API
#' @description
#' \code{decorateGmtByUnderOvenAndNoise}
#'
#' \code{decorateGmtByUnderOvenAndNoise} decorates GO with labels (over, under, noise) per term.
#'
#' @param input_gmt input dataframe, read from gmt file.
#' @param number_of_over_representation_groups set the number of groups
#' which will be chosen to over representation.
#' @param number_of_under_representation_groups set the number of groups
#' which will be chosen to under representation.
#' @return Return data frame with model from specific location.
#' @noRd
decorateGmtByUnderOvenAndNoise <- function(input_gmt,
                                           number_of_over_representation_groups = 1,
                                           number_of_under_representation_groups = 0) {
  # Initialize all by noise labels.
  sample_label <- rep('noise', length(input_gmt$ontologyId))
  gmt_for_generator <-
    data.frame(input_gmt, "sample_label" = sample_label)
  
  # Choose and label terms for over and under representation.
  go_size <- length(gmt_for_generator$listOfValues)
  size_of_over_under_repr <-
    number_of_over_representation_groups + number_of_under_representation_groups
  go_change_repr <-
    sample(1:go_size, size_of_over_under_repr, replace = FALSE)
  
  over_under_label <-
    c(
      rep('over', number_of_over_representation_groups),
      rep('under', number_of_under_representation_groups)
    )
  
  terms_to_manipulation <- data.frame('term_id' = go_change_repr,
                                      'over_under_label' = over_under_label)
  
  for (i in 1:length(terms_to_manipulation$term_id)) {
    term_row <- terms_to_manipulation[i, ]
    gmt_for_generator[term_row$term_id, ]$sample_label <-
      term_row$over_under_label
  }
  
  return(gmt_for_generator)
}


#' Convert list to GMT data frame
#'
#' Convert ontology representation from list to gmt data frame.
#'
#' @param gmt_list List with element names as ontologyId and genes in each
#' element.
#' @return Return data frame with model.
#' @export
list_to_gmt <- function(gmt_list) {
  listAsGmtDataFrame <-
    plyr::ldply(
      .data = gmt_list,
      .id = c('ontologyId'),
      .fun = function(element) {
        print(element)
        ontology_name <- stringi::stri_rand_strings(length = 5, n = 1)
        data.frame(
          'ontologyName' = ontology_name,
          'listOfValues' = I(list(element)),
          stringsAsFactors = FALSE
        )
      }
    )
  return(listAsGmtDataFrame)
}


# PUBLIC API
#' @description
#' \code{generateInputSamples} Generates artificial GO with specific terms (under or over represented).
#'
#' @param input_gmt input dataframe, read from gmt file.
#' @param noise_ratio level of noise in the sample, value from 0 to 1.
#' @param group_under_over_representation_ratio ratio of over represented group.
#' @param number_of_over_representation_groups number of groups chosen to over representation.
#' @param number_of_under_representation_groups number of groups chosen to under representation.
#' @return Return data frame with model from specific location.
#' @noRd
generateInputSamples <-
  function(input_gmt_decorated,
           noise_ratio = 0.2,
           over_repr_ratio = 0.5,
           under_repr_ratio = 0.05,
           rand_from_unique = TRUE,
           number_of_samples = 1) {
    all_genes_in_ontology <- NULL
    all_genes_in_enrichment <- NULL
    if (rand_from_unique) {
      all_genes_in_ontology <-
        unique(unlist(input_gmt_decorated$listOfValues))
      all_genes_in_enrichment <- unique(unlist(input_gmt_decorated[input_gmt_decorated$sample_label == 'over', ]$listOfValues))
    } else {
      all_genes_in_ontology <- unlist(input_gmt_decorated$listOfValues)
      all_genes_in_enrichment <- unlist(input_gmt_decorated[input_gmt_decorated$sample_label == 'over', ]$listOfValues)
    }
    
    size_of_ontology <- length(all_genes_in_ontology)
    size_of_noise <- ceiling(size_of_ontology * noise_ratio)
    
    size_of_enrichment <-
      ceiling(length(all_genes_in_enrichment) * over_repr_ratio)
    
    samples <- vector("list", number_of_samples)
    for (i in 1:length(samples)) {
      sample_noise <- all_genes_in_ontology[sample(1:length(all_genes_in_ontology), size_of_noise, replace =
                                                     FALSE)]
      sample_enrichment <- all_genes_in_enrichment[sample(1:length(all_genes_in_enrichment),
                                                          size_of_enrichment,
                                                          replace = FALSE)]
      #samples[[i]] <- unique(c(sample_noise, sample_enrichment))
      samples[[i]] <- list('sample_noise' = sample_noise, 'sample_enrichment' = sample_enrichment)
    }

    return(samples)
  }

# PUBLIC API
#' @description
#' \code{getMultipleTestsSummary}
#'
#' \code{getMultipleTestsSummary} generate artificial GO with specific terms under or over represented.
#'
#' @param tests_res list of multiple tests results.
#' @param comparison_col_name column name which indicated data to compare on.
#' @param labels label datatable by additional columns with values.
#' @param cut_off threshold for value selected by comparison_col_name
#' @return Return data frame with FDR. TPRs per test.
#' @noRd
#' @importFrom magrittr %>%
#' @importFrom data.table :=
getMultipleTestsSummary <- function(tests_res,
                                    comparison_col_name,
                                    labels = list(),
                                    cut_off = 0.05) {
  # Summarize results.
  print("Mulea sumary time:")
  tictoc::tic()
  
  metadata_len <- length(tests_res[[1]]$metadata)
  
  
  sumary_res <-
    data.frame(matrix(ncol = 10 + metadata_len, nrow = 0))
  colnames(sumary_res) <- c(
    'test_no',
    'TP',
    'TP_size',
    'FP',
    'FP_size',
    'FN',
    'FN_size',
    'TN',
    'TN_size',
    'over_repr_terms',
    names(tests_res[[1]]$metadata)
  )
  number_of_tests <- length(tests_res)
  for (i in 1:number_of_tests) {
    # Actual condition
    # Total population = P + N
    total_population <- tests_res[[i]]$test_data$ontologyId
    total_population_size <- length(total_population)
    # Positive (P)
    P <- tests_res[[i]]$test_data[tests_res[[i]]$test_data$sample_label == 'over', ]$ontologyId
    P_size <- length(P)
    # Negative (N)
    N <- tests_res[[i]]$test_data[tests_res[[i]]$test_data$sample_label != 'over', ]$ontologyId
    N_size <- length(N)
    if (P_size + N_size != total_population_size) {
      warning("Not OK size of Actual in contingency table")
    }
    
    # Predicted condition
    # Predicted Positive (PP)
    PP <- tests_res[[i]]$mulea_res[tests_res[[i]]$mulea_res[, comparison_col_name] <= cut_off,]$ontology_id
    PP_size <- length(PP)
    # Predicted Negative (PN)
    PN <- tests_res[[i]]$mulea_res[tests_res[[i]]$mulea_res[, comparison_col_name] > cut_off,]$ontology_id
    PN_size <- length(PN)
    if (PP_size + PN_size != total_population_size) {
      warning("Not OK size of Predicted in contingency table")
    }
    
    # True positive (TP) : hit
    TP <- intersect(P, PP)
    TP_size <- length(TP)
    # False positive (FP) : type I error, false alarm, overestimation
    FP <- intersect(N, PP)
    FP_size <- length(FP)
    # False negative (FN) : type II error, miss, underestimation
    FN <- intersect(P, PN)
    FN_size <- length(FN)
    # True negative (TN) : correct rejection
    TN <- intersect(N, PN)
    TN_size <- length(TN)
    
    if (TP_size + FP_size + FN_size + TN_size != total_population_size) {
      warning("Not OK size of total  contingency table")
    }
    
    over_repr_terms <- tests_res[[i]]$test_data[tests_res[[i]]$test_data$sample_label == 'over', ]$ontologyId
    
    sumary_res_tmp <- data.frame(
      'test_no' = i,
      'TP' = I(list(TP)),
      'TP_size' = TP_size,
      'FP' = I(list(FP)),
      'FP_size' = FP_size,
      'FN' = I(list(FN)),
      'FN_size' = FN_size,
      'TN' = I(list(TN)),
      'TN_size' = TN_size,
      'over_repr_terms' = I(list(over_repr_terms))
    )
    
    for (metadata_entry in names(tests_res[[i]]$metadata)) {
      if ('input_select' == metadata_entry) {
        sumary_res_tmp <- cbind(sumary_res_tmp,
                                'input_select' = I(tests_res[[i]]$metadata[metadata_entry]))
      } else {
        sumary_res_tmp <- cbind(sumary_res_tmp,
                                metadata_entry = as.character(tests_res[[i]]$metadata[metadata_entry]))
      }
    }
    
    sumary_res[i,] <- sumary_res_tmp
  }
  
  sumary_res <- tibble::tibble(sumary_res) %>%
    dplyr::mutate(FPR = FP_size / (FP_size + TN_size)) %>%
    dplyr::mutate(TPR = TP_size / (TP_size + FN_size)) %>%
    dplyr::mutate(FDR = FP_size / (TP_size + FP_size)) %>%
    dplyr::mutate(NPV = TN_size / (FN_size + TN_size))
  
  for (label_id in seq_along(labels)) {
    # IMPORTANT : Labels are as characters in datatable
    label_name <- as.character(names(labels)[[label_id]])
    label_value <- as.character(labels[[label_id]])
    sumary_res <-
      sumary_res %>% dplyr::mutate(!!label_name := label_value)
  }
  
  tictoc::toc()
  
  return(sumary_res)
}

# PUBLIC API
#' @description
#' \code{getSummaryToRoc}
#'
#' \code{getSummaryToRoc} generate artificial GO with specific terms under or over represented.
#'
#' @param tests_res list of multiple tests results.
#' @return Return data frame which is the base to count ROC.
#' @noRd
#' @importFrom magrittr %>%
#' @importFrom plyr .
#' @importFrom rlang .data
getSummaryToRoc <- function(tests_res,
                            cut_off_resolution = 0.01,
                            methods_names = c("p_value", "adjusted_p_value", "eFDR")) {
  print("Mulea ROC data calculation time:")
  tictoc::tic()
  
  number_of_tests <- length(tests_res)
  data_to_roc <- data.frame(
    "sample_label" = c(),
    "pValue" = c(),
    "adjustedPValue" = c(),
    "adjustedPValueEmpirical" = c(),
    "noise_ratio" = c()
  )
  for (i in 1:number_of_tests) {
    # tests_res[[i]]$mulea_res[, c("p_value", "adjusted_p_value", "eFDR")]
    data_to_roc <-
      rbind(
        data_to_roc,
        data.frame("sample_label" = 
                     tests_res[[i]]$test_data[, c("sample_label")],
                   tests_res[[i]]$mulea_res[, c("p_value", 
                                                "adjusted_p_value",
                                                "eFDR")],
                   "noise_ratio" = 
                     rep(tests_res[[i]]$metadata$noise_ratio, 
                         times=(nrow(tests_res[[i]]$test_data)))
                   )
      )
  }
  
  roc_stats <- tibble::tibble(
    TP_val = numeric(),
    TN_val = numeric(),
    FP_val = numeric(),
    FN_val = numeric(),
    TPR = numeric(),
    FPR = numeric(),
    PPV = numeric(),
    f1_score = numeric(),
    sum_test = numeric(),
    cut_off = numeric(),
    method = character(),
    noise_ratio = numeric()
  )
  
  for (method_name in methods_names) {
    for (cut_off in seq(0, 1, cut_off_resolution)) {
      sim_mult_tests_res_to_roc_summary <- data_to_roc %>%
        dplyr::mutate(PP = !!as.name(method_name) <= cut_off) %>%
        dplyr::mutate(TP = (.data$PP == TRUE & .data$sample_label == 'over'),
                      TN = (.data$PP == FALSE & .data$sample_label != 'over'),
                      FP = (.data$PP == TRUE & .data$sample_label != 'over'),
                      FN = (.data$PP == FALSE & .data$sample_label == 'over')
        )

      sim_sum <- 
        sim_mult_tests_res_to_roc_summary %>% dplyr::summarise(
           TP_val = sum(.data$TP),
           TN_val = sum(.data$TN),
           FP_val = sum(.data$FP),
           FN_val = sum(.data$FN)
         )
      
      sim_sum_new <- sim_mult_tests_res_to_roc_summary %>% 
        dplyr::group_by(noise_ratio) %>% 
        dplyr::summarise(
          TP_val = sum(.data$TP),
          TN_val = sum(.data$TN),
          FP_val = sum(.data$FP),
          FN_val = sum(.data$FN)
      )
      
      sim_sum_roc <- sim_sum %>% dplyr::mutate(
        TPR = .data$TP_val / (.data$TP_val + .data$FN_val),
        FPR = .data$FP_val / (.data$FP_val + .data$TN_val),
        PPV = .data$TP_val / (.data$TP_val + .data$FP_val),
        f1_score = .data$TP_val / (.data$TP_val + 0.5 * (.data$FP_val + .data$FN_val)),
        sum_test = .data$TP_val + .data$TN_val + .data$FP_val + .data$FN_val,
        cut_off = cut_off,
        method = method_name,
        noise_ratio = NA_real_
      )
      
      sim_sum_roc_new <- sim_sum_new %>% dplyr::mutate(
        TPR = .data$TP_val / (.data$TP_val + .data$FN_val),
        FPR = .data$FP_val / (.data$FP_val + .data$TN_val),
        PPV = .data$TP_val / (.data$TP_val + .data$FP_val),
        f1_score = .data$TP_val / (.data$TP_val + 0.5 * (.data$FP_val + .data$FN_val)),
        sum_test = .data$TP_val + .data$TN_val + .data$FP_val + .data$FN_val,
        cut_off = cut_off,
        method = method_name
      )
      
      roc_stats <- roc_stats %>% 
        tibble::add_row(sim_sum_roc) %>% 
        tibble::add_row(sim_sum_roc_new)
    }
  }
  
  tictoc::toc()
  
  return(roc_stats)
}

# PUBLIC API
#' @description
#' \code{getMultipleTestsSummaryAcrossCutOff}
#'
#' \code{getMultipleTestsSummaryAcrossCutOff} doing summary across cutoff range.
#'
#' @param tests_res list of multiple tests results.
#' @param cut_off_range threshold for value selected by comparison_col_name
#' @return Return data frame with FDR. TPRs per test.
#' @noRd
getMultipleTestsSummaryAcrossCutOff <- function(tests_res,
                                                cut_off_range = seq(0, 1, 0.1),
                                                prod_run = TRUE) {
  tests_res_sum <- NULL
  for (cut_off in cut_off_range) {
    print(cut_off)
    tests_res_sum_p <- getMultipleTestsSummary(
      tests_res = tests_res,
      comparison_col_name = 'p_value',
      labels = list('method' = 'p', 'cut_off' = cut_off),
      cut_off = cut_off
    )
    
    tests_res_sum_bh <- getMultipleTestsSummary(
      tests_res = tests_res,
      comparison_col_name = 'adjusted_p_value',
      labels = list('method' = 'bh', 'cut_off' = cut_off),
      cut_off = cut_off
    )
    
    tests_res_sum_pt <- getMultipleTestsSummary(
      tests_res = tests_res,
      comparison_col_name = 'eFDR',
      labels = list('method' = 'pt', 'cut_off' = cut_off),
      cut_off = cut_off
    )
    
    if (prod_run) {
      columns_to_plot <- c('test_no', 'FPR', 'TPR', 'FDR', 'NPV', 
                           'method', 'cut_off', 'noise_ratio', 
                           'number_of_tests')
      tests_res_sum <- rbind(tests_res_sum,
                             tests_res_sum_p[, columns_to_plot],
                             tests_res_sum_pt[, columns_to_plot],
                             tests_res_sum_bh[, columns_to_plot])
    } else {
      tests_res_sum <- rbind(tests_res_sum,
                             tests_res_sum_p,
                             tests_res_sum_pt,
                             tests_res_sum_bh)
    }
  }
  return(tests_res_sum)
}


# PUBLIC API
#' @description
#' \code{simulateMultipleTests}
#'
#' \code{simulateMultipleTests} generate artificial GO with specific terms under or over represented.
#'
#' @param input_gmt_filtered gmt data frame with ontology for tests.
#' @param number_of_tests number of tests to perform.
#' @param noise_ratio ratio of noise in data from [0,1] interval.
#' @param number_of_over_representation_groups number of terms to over represent.
#' @param number_of_under_representation_groups number of terms to under represent.
#' @return Return data frame with FDR. TPRs per test.
#' @noRd
simulateMultipleTests <- function(input_gmt_filtered,
                                  number_of_tests = 10,
                                  noise_ratio = 0.35,
                                  over_repr_ratio  = 0.5,
                                  number_of_over_representation_groups = ceiling(nrow(input_gmt_filtered) *
                                                                                   0.1),
                                  number_of_under_representation_groups = 0,
                                  number_of_steps = 5000,
                                  nthreads = 16) {
  print("Mulea calculation time:")
  tictoc::tic()
  number_of_samples <- 1
  tests_res <- vector("list", number_of_tests)
  for (i in 1:number_of_tests) {
    print(i)
    
    input_gmt_decorated <-
      decorateGmtByUnderOvenAndNoise(
        input_gmt = input_gmt_filtered,
        number_of_over_representation_groups = number_of_over_representation_groups,
        number_of_under_representation_groups = number_of_under_representation_groups
      )
    
    samples <- generateInputSamples(
      input_gmt_decorated,
      noise_ratio = noise_ratio,
      over_repr_ratio = over_repr_ratio,
      number_of_samples = number_of_samples
    )
    
    if (length(samples) != 1) {
      warning("sample is not size 1")
    }
    
    input_select <- unique(unlist(samples))
    
    mulea_ora_model <- MulEA::ora(
      gmt = input_gmt_filtered,
      element_names = input_select,
      p_value_adjustment_method = "eFDR",
      number_of_permutations = number_of_steps,
      number_of_cpu_threads = nthreads
    )
    
    mulea_ora_model_bh <- MulEA::ora(
      gmt = input_gmt_filtered,
      element_names = input_select,
      p_value_adjustment_method = "BH",
      number_of_permutations = number_of_steps,
      number_of_cpu_threads = nthreads
    )
    
    mulea_ora_results <- MulEA::run_test(mulea_ora_model)
    mulea_ora_results_bh <- MulEA::run_test(mulea_ora_model_bh)
    mulea_ora_results$adjusted_p_value <- mulea_ora_results_bh$adjusted_p_value
    
    tests_res[[i]]$mulea_res <- mulea_ora_results
    tests_res[[i]]$test_data <- input_gmt_decorated
    tests_res[[i]]$metadata <- list(
      'noise_ratio' = noise_ratio,
      'number_of_tests' = number_of_tests,
      'over_repr_ratio' = over_repr_ratio,
      'number_of_over_representation_groups' = number_of_over_representation_groups,
      'input_select' = input_select
    )
  }
  tictoc::toc()
  return(tests_res)
}


# PUBLIC API
#' @description
#' \code{simulateMultipleTestsWithRatioParam}
#'
#' \code{simulateMultipleTestsWithRatioParam} generate artificial GO with specific terms under or over represented.
#'
#' @param input_gmt_filtered gmt data frame with ontology for tests.
#' @param number_of_tests number of tests to perform.
#' @param noise_ratio_range range of ratios of noise in data from [0,1] interval.
#' @param number_of_over_representation_groups number of terms to over represent.
#' @param number_of_steps number of steps.
#' @return Return data frame with FDR. TPRs per test.
#' @noRd
simulateMultipleTestsWithRatioParam <- function(input_gmt_filtered,
                                                noise_ratio_range = seq(0.1, 0.5, 0.1),
                                                number_of_tests = 100,
                                                over_repr_ratio = 0.5,
                                                number_of_over_representation_groups = ceiling(nrow(input_gmt_filtered) *
                                                                                                 0.2),
                                                number_of_steps = 10000,
                                                nthreads = 16) {
  tictoc::tic()
  sim_mult_tests <- list()
  for (noise_ratio in noise_ratio_range) {
    print("noise_ratio")
    print(noise_ratio)
    sim_mult_tests <- c(
      sim_mult_tests,
      simulateMultipleTests(
        input_gmt_filtered = input_gmt_filtered,
        number_of_tests = number_of_tests,
        noise_ratio = noise_ratio,
        over_repr_ratio = over_repr_ratio,
        number_of_over_representation_groups = number_of_over_representation_groups,
        number_of_under_representation_groups = 0,
        number_of_steps = number_of_steps,
        nthreads = nthreads
      )
    )
  }
  print('MulEA : ratio search calculation time:')
  tictoc::toc()
  return(sim_mult_tests)
}
