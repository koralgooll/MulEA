checkIfPoolIncludeSample <-
  function(model, sampleVector, poolVector = NULL) {
    # Chcking if experiment data all in model data.
    if (0 != length(poolVector)) {
      if (0 != sum(!(sampleVector %in% poolVector))) {
        warning("element_names are outside of pool.",
                " ",
                paste(setdiff(sampleVector, unique(poolVector)), collapse = ", "))
        return(setdiff(sampleVector, setdiff(sampleVector, unique(poolVector))))
      }
    } else {
      if (0 != sum(!(sampleVector %in% unique(unlist(
        model$listOfValues
      ))))) {
        warning("element_names are outside of gmt.",
                " ",
                paste(setdiff(sampleVector, unique(
                  unlist(model$listOfValues)
                )),
                collapse = ", "))
        return(setdiff(sampleVector, setdiff(sampleVector, unique(
          unlist(model$listOfValues)
        ))))
      }
    }
    return(sampleVector)
  }


cutGmtToPool <- function(gmt, pool) {
  cutDF <- plyr::ddply(
    .data = gmt,
    .variables = c("ontologyId"),
    .fun = function(dfRow) {
      dfRow$listOfValues[[1]] <- intersect(dfRow$listOfValues[[1]], pool)
      dfRow
    }
  )
  cutDF
}
