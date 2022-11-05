#####################################################################################################
# Function for FDR corrected hypergeometric enrichment test
#####################################################################################################

# snow and rlecuyer packages should be installed


## Arguments of HyperGeomFDR function:
#	steps:		the rounds of simulations (a single number)
#	pool:		background genes (character vector)
#	select:		genes to investigate (character vector)
#	DB: 		the genes set used for enrichment analysis (character list)
#	nthreads:	number of threads to use (a single number)

## Description of the hypergeometric test
# phyper(q, m, n, k, lower.tail = TRUE, log.p = FALSE)
# Arguments:
#        q: vector of quantiles representing the number of white balls
#           drawn without replacement from an urn which contains both
#           black and white balls.
#        m: the number of white balls in the urn.
#        n: the number of black balls in the urn.
#        k: the number of balls drawn from the urn.
#
# x=length(intersect(select,DB_i))    	#Number of common genes between DB and select
# m=length(intersect(pool,DB_i))        #Number of common genes between DB and pool
# n=length(pool)-length(intersect(pool,DB_i))     #Number of non-pool genes among DB (setdiff)
# k=length(select)                    	#Number of genes in select
# P_val=dhyper(length(intersect(select,DB_i)), length(intersect(pool,DB_i)), length(pool)-length(intersect(pool,DB_i)), length(select))
#
# wikipedia
# N = length(pool)
# K = length(intersect(pool,DB_i))
# n = length(select)
# k = length(intersect(select,DB_i))
#
# (choose(length(intersect(pool,DB_i)),length(intersect(select,DB_i)))*
#      choose(length(pool)-length(intersect(pool,DB_i)),length(select)-length(intersect(select,DB_i))))/
#      choose(length(pool),length(select))

# m= intersect(select,DB_i)
# n= intersect(bg,DB_i)-m
# k=length (DB_i)


# P_val=(choose(length(intersect(BG,DB_i)),length(intersect(IG,DB_i)))*choose(length(BG)-length(intersect(BG,DB_i)),length(IG)-length(intersect(IG,DB_i))))/choose(length(BG),length(IG))

#' PRIVATE class : An S4 class to represent a Hypergeometric tests in Mulea.
#'
#' @slot gmt A data.frame representing GMT's reprezentation of model.
#' @slot element_names A data from expeciment to analize accross model.
#' @slot pool A background data to count test.
#' @slot nthreads Number of processor's threads used in calculations.
#' @return SetBasedEnrichmentTest object. Used as private function.
#' @examples
#' \dontrun{
#' #It is a private s4 object. Look at ora's examples.
#' }
SetBasedEnrichmentTest <- setClass(
  "SetBasedEnrichmentTest",
  slots = list(
    gmt = "data.frame",
    element_names = "character",
    pool = "character",
    number_of_permutations = "numeric",
    only_hyper_geometric_test = "logical",
    nthreads = "numeric",
    test = "function"
  )
)

setMethod("initialize", "SetBasedEnrichmentTest",
          function(.Object,
                   gmt = data.frame(),
                   element_names = character(),
                   pool = character(),
                   number_of_permutations = 10000,
                   test = NULL,
                   only_hyper_geometric_test = FALSE,
                   nthreads = 4,
                   ...) {
            .Object@gmt <- gmt
            .Object@element_names <- element_names
            .Object@pool <- pool
            .Object@number_of_permutations <- number_of_permutations
            .Object@only_hyper_geometric_test <-
              only_hyper_geometric_test
            .Object@nthreads <- nthreads
            
            .Object@test <- function(model) {
              pool <- NULL
              if (0 == length(model@pool)) {
                pool <- unique(unlist(.Object@gmt[, 'listOfValues']))
              } else {
                pool <- unique(model@pool)
              }
              
              element_names <- .Object@element_names
              if (!all(element_names %in% pool)) {
                element_names <- element_names[element_names %in% pool]
                warning(
                  "Not all elements of element_names (sample) are from pool. ",
                  "TestData vector is automatically cut off to pool vector."
                )
              }
              
              DB <- .Object@gmt[, 'listOfValues']
              names(DB) <- .Object@gmt$ontologyId
              
              testResults <-
                set.based.enrichment.test.wrapper(
                  steps = .Object@number_of_permutations,
                  pool = pool,
                  select = element_names,
                  DB = DB,
                  only_hyper_geometric_test =
                    model@only_hyper_geometric_test,
                  nthreads = model@nthreads
                )
              testResults
            }
            
            .Object
          })

#' @describeIn SetBasedEnrichmentTest runs test calculations.
#' @param model Object of s4 class represents Mulea Test.
#' @return run_test method for SetBasedEnrichmentTest object. Used as private function.
#' @examples
#' \dontrun{
#' #It is a private method. Look at run_test of ora's examples.
#' }
setMethod("run_test",
          signature(model = "SetBasedEnrichmentTest"),
          function(model) {
            model@test(model)
          })


set.based.enrichment.test.wrapper = function(steps,
                                             pool,
                                             select,
                                             DB,
                                             nthreads = 4,
                                             only_hyper_geometric_test = FALSE) {
  # print("nthreads")
  # print(nthreads)
  # print("steps")
  # print(steps)
  # print("pool")
  # print(pool)
  # print("select")
  # print(select)
  # print("DB")
  # print(DB)
  
  setEnrTestRes <-
    set.based.enrichment.test(
      steps = steps,
      pool = pool,
      select = select,
      DB = DB,
      nthread = nthreads
    )
  return(setEnrTestRes)
}
