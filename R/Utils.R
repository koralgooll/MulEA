
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