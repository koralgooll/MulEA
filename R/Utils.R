
# PUBLIC API
#' @description
#' \code{readGmtFileAsDataFrame}
#'
#' \code{readGmtFileAsDataFrame} read model in data frame form from gmt file.
#'
#' @param gmtFilePath path with name of file, where the file is localized or where to save model. Example: "/R/MulEA/extdata/model.gmt"
#'
#'
#' @title Input/Output Functions
#' @name  InputOutputFunctions
#' @export
#'
#' @return Return data frame with model from specific location.
readGmtFileAsDataFrame <- function(gmtFilePath) {
    fileConnection <- file(gmtFilePath)
    tryCatchRes <- tryCatch(
        lines <- readLines(fileConnection),
        warning=function(w) {}
    )
    close(fileConnection)
    
    lines <- lines[!grepl('^#+', lines, fixed = FALSE)]
    lines <- lines["" != lines]
    
    gmtAsDF <- plyr::adply(.data = lines, .margins = 1, .fun = function(line) {
        fields <- strsplit(line, split = "\t")[[1]]
        category <- fields[1]
        if (startsWith(fields[2], "\"") && endsWith(fields[2], "\"")) {
            description <- fields[2]
        } else {
            description <- paste("\"", fields[2], "\"", sep = "")
        }
        listOfValues <- fields[3:length(fields)]
        data.frame('ontologyId' = category, 'ontologyName' = description, 'listOfValues' = I(list(listOfValues)), stringsAsFactors = FALSE)
    })
    gmtAsDF[c("ontologyId", "ontologyName", "listOfValues")]
}

#TODO : Is that hepler needed?
readGmtFileAsPlaneDF <- function(gmtFilePath) {
    maxColLength <- max(count.fields(gmtFilePath, sep = '\t', quote = "\""))
    model <- read.table(file = gmtFilePath, header = FALSE, fill = TRUE,
                        stringsAsFactors = FALSE, sep = "\t", strip.white = TRUE,
                        col.names = paste0("V",seq_len(maxColLength)), quote = "\"")
    model
}

# PUBLIC API
#' @description
#' \code{saveDataFrameAsGmtFile}
#'
#' \code{saveDataFrameAsGmtFile} saves copy of the model from data frame in gmt file.
#'
#' @param modelDF data frame with model.
#'
#' @rdname InputOutputFunctions
#' @export
#'
#' @return Return gmt file under specific location which include model in gmt format.
#' @examples
#' modelDfFromFile <- MulEA::readGmtFileAsDataFrame(gmtFilePath = system.file(package="MulEA", "extdata", "model.gmt"))
#' MulEA::saveDataFrameAsGmtFile(modelDF = modelDfFromFile, gmtFilePath = paste(system.file(package="MulEA", "extdata"), "fromDb.gmt", sep = "/"))
saveDataFrameAsGmtFile <- function(modelDF, gmtFilePath) {
    vectorOfModel <- plyr::daply(.data = modelDF, .variables = c("ontologyId"), .fun = function(dataFrameRow){
        collapsedListOfValues <- paste(dataFrameRow[,3][[1]], collapse = "\t")
        paste(dataFrameRow[1], dataFrameRow[2], collapsedListOfValues, sep = "\t")
    })
    fileConnection <- file(gmtFilePath)
    writeLines(vectorOfModel, con = fileConnection, sep = "\n", useBytes = FALSE)
    close(fileConnection)
}

# PUBLIC API
#' @description
#' \code{readEdbFileAsDataFrame}
#'
#' \code{readEdbFileAsDataFrame} read GSEA results in data frame form from .edb file.
#'
#' @param edbFilePath path with name of file, where the file is localized or where to save model.
#'
#'
#' @title Input/Output Functions
#' @name  InputOutputFunctions
#' @export
#'
#' @return Return data frame with model from specific location.
readEdbFileAsDataFrame <- function(edbFilePath) {
    
    xml_parsed <- XML::xmlTreeParse(edbFilePath, useInternalNodes = T)
    geneset_list <- XML::xpathApply(xml_parsed, path = '/EDB/DTG/@GENESET')
    es_list <- XML::xpathApply(xml_parsed, path = '/EDB/DTG/@ES')
    nes_list <- XML::xpathApply(xml_parsed, path = '/EDB/DTG/@NES')
    fdr_list <- XML::xpathApply(xml_parsed, path = '/EDB/DTG/@FDR')
    np_list <- XML::xpathApply(xml_parsed, path = '/EDB/DTG/@NP')
    
    gsea_bi_res_df <- data.frame('ontology_id'=c(), 'es'=c(), 'nes'=c(), 'fdr'=c(), 'np'=c())
    
    
    for (i in 1:length(geneset_list)) {
        onlology_id <- strsplit(geneset_list[[i]][['GENESET']], '#')[[1]][2]
        es <- es_list[[i]][['ES']]
        nes <- nes_list[[i]][['NES']]
        fdr <- fdr_list[[i]][['FDR']]
        np <- np_list[[i]][['NP']]
        gsea_bi_res_df_i <- data.frame('ontology_id'=onlology_id, 'es'=es, 'nes'=nes, 'fdr'=fdr, 'np'=np)
        gsea_bi_res_df <- rbind(gsea_bi_res_df, gsea_bi_res_df_i)
    }
    
    gsea_bi_res_df
}


# PUBLIC API
# TODO : Add quantile parameters as separate! Nothing do is default.
#' @description
#' \code{filterOntology}
#'
#' \code{filterOntology} cut ontology to specific terms sizes.
#'
#' @param input_gmt input dataframe, read from gmt file.
#' @param min minimum size of term. Default 20% from quantile on term size distribution. 
#' @param max maximum size of term. Default 80% from quantile on term size distribution.
#'
#'
#' @title Input/Output Functions
#' @name  InputOutputFunctions
#' @export
#'
#' @return Return data frame with model from specific location.
filterOntology <- function(input_gmt, min=NULL, max=NULL) {
    if (is.null(min)) {
        terms_sizes <- plyr::laply(.data = input_gmt$listOfValues, .fun = function(term) {
            length(term)
        })
        term_size_dist_q <- quantile(terms_sizes, probs = seq(0, 1, 0.1), type = 2, na.rm = FALSE)
        
        min = term_size_dist_q['20%']
    }
    
    if (is.null(max)) {
        terms_sizes <- plyr::laply(.data = input_gmt$listOfValues, .fun = function(term) {
            length(term)
        })
        term_size_dist_q <- quantile(terms_sizes, probs = seq(0, 1, 0.1), type = 2, na.rm = FALSE)
        max = term_size_dist_q['80%']
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


# PUBLIC API
#' @description
#' \code{generateInputData}
#'
#' \code{generateInputData} generate artificial GO with specific terms under or over represented.
#'
#' @param input_gmt input dataframe, read from gmt file.
#' @param sample_ratio  
#' @param group_under_over_representation_ratio
#' @param number_of_over_representation_groups
#' @param number_of_under_representation_groups
#'
#'
#' @title Input/Output Functions
#' @name  InputOutputFunctions
#' @export
#'
#' @return Return data frame with model from specific location.
generateInputData <- function(input_gmt, sample_ratio=0.5,
                              #TODO : Separate under_over parameter int two.
                                group_under_over_representation_ratio=0.1,
                                number_of_over_representation_groups=1, 
                                number_of_under_representation_groups=1, turn_on_log=FALSE) {
    
    # Initialize all by noise size. 
    sample_label <- rep('noise', length(input_gmt$ontologyId))
    sample_values <- rep(c(''), length(input_gmt$ontologyId))
    gmt_for_generator <- data.frame(input_gmt, "sample_label"=sample_label, "listOfSampleValues"=sample_values)
    
    # Choose terms for over and under representation.
    go_size <- length(gmt_for_generator$listOfValues)
    size_of_over_under_repr <- number_of_over_representation_groups + number_of_under_representation_groups
    go_change_repr <- sample(1:go_size, size_of_over_under_repr, replace=FALSE)
    
    over_under_label <- c(rep('over', number_of_over_representation_groups), rep('under', number_of_under_representation_groups))
    over_under_label <- sample(over_under_label)
    terms_to_manipulation <- data.frame('term_id' = go_change_repr, 'over_under_label' = over_under_label)
    
    for (i in 1:length(terms_to_manipulation$term_id)) {
        term_row <- terms_to_manipulation[i,]
        gmt_for_generator[term_row$term_id,]$sample_label <- term_row$over_under_label
    }
    
    # Permute all rows in gmt_for_generator.
    permutation_mask <- sample(1:go_size)
    gmt_for_generator <- gmt_for_generator[permutation_mask,]
    
    
    for (i in 1:length(gmt_for_generator$ontologyId)) {
        term_row <- gmt_for_generator[i,]
        term_size <- length(term_row$listOfValues[[1]])
        if (term_row$sample_label == "over") {
            over_repr_ratio <- sample_ratio + group_under_over_representation_ratio
            # Think to randomly floor or ceiling.
            size_of_sample <- floor(term_size * over_repr_ratio)
        }
        if (term_row$sample_label == "under") {
            under_repr_ratio <- sample_ratio - group_under_over_representation_ratio
            size_of_sample <- floor(term_size * under_repr_ratio)
        }
        if (term_row$sample_label == "noise") {
            size_of_sample <- floor(term_size * sample_ratio)
        }
        
        input_select_indicator <- sample(1:term_size, size_of_sample, replace=FALSE)
        input_select_term <- term_row$listOfValues[[1]][input_select_indicator]
        
        gmt_for_generator[i,]$listOfSampleValues <- I(list(input_select_term))
    }
    
    input_select <- unique(unlist(gmt_for_generator$listOfSampleValues))
    
    to_return <- list()
    to_return$input_select <- input_select
    to_return$gmt_for_generator <- gmt_for_generator
    
    return(to_return)
}
