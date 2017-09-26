
setClass("pathNetTest",
         slots = list(
           # testData = "character",
           testData = "numeric",
           scores = "numeric",
           graph = "data.frame",
           test = "function",
           p = "numeric",
           numberOfPermutations = "numeric"
         ))

setMethod("initialize", "pathNetTest",
          function(.Object,
                   # testData = character(),
                   testData = numeric(),
                   scores = numeric(),
                   graph = data.frame(),
                   test = NULL,
                   p = 1,
                   numberOfPermutations = 1000,
                   ...) {

            .Object@testData <- testData
            .Object@scores <- scores
            .Object@graph <- graph
            .Object@p <- p
            .Object@numberOfPermutations <- numberOfPermutations

            .Object@test <- function(dataObject, testObject) {

              directEvidenceData <- cbind(testObject@testData, testObject@scores)
              graphObject <- igraph::graph_from_data_frame(testObject@graph, directed = TRUE, vertices = NULL)
              adjacencyMatrix <- as.matrix(get.adjacency(graphObject))

              results <- PathNet(Enrichment_Analysis = TRUE,
                                 DirectEvidence_info = directEvidenceData,
                                 Adjacency = adjacencyMatrix,
                                 pathway = testObject@graph,
                                 Column_DirectEvidence = 2,
                                 n_perm = 2, threshold = 0.05)

              # listmodelDfFromFile <- dataObject@gmt$listOfValues
              # names(listmodelDfFromFile) <- dataObject@gmt$ontologyId
              #
              # samplesToAnalisys <- testObject@scores
              # names(samplesToAnalisys) <- testObject@testData
              #
              # fgseaRes <- fgsea(pathways = listmodelDfFromFile,
              #                   stats = samplesToAnalisys,
              #                   gseaParam = testObject@p, nperm = testObject@numberOfPermutations)
              #
              # resultDf <- merge(dataObject@gmt, fgseaRes, by.x = "ontologyId", by.y = "pathway", all = TRUE)[c("ontologyId", "ontologyName", "listOfValues", "pval")]
              # colnames(resultDf) <- c("ontologyId", "ontologyName", "listOfValues", "p.value")
              # resultDf
              results
            }

            .Object

          })
