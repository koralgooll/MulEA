
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
