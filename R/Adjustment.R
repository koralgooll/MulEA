
Adjustment <- setClass("Adjustment",
                         slots = list(
                           setBasedEnrichmentTest = "function"
                         ))

setMethod("initialize", "Adjustment",
          function(.Object,
                   setBasedEnrichmentTest = NULL,
                   ...) {

            .Object@setBasedEnrichmentTest <- function() {
                print("SET BASE ENRICHMENT TESR WILL WORK HERE")
            }

            .Object
          })
