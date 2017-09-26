setClass("muleaData",
         slots = list(
           gmt = "data.frame",
           taxID = "character",
           ontologyName = "character",
           version = "character",
           scientificName = "character",
           commonEnglishName = "character",
           description = "character",
           modelBinaryMatrix = "data.frame" #TODO : this will be removed, privat
         ))

setMethod("initialize", "muleaData",
          function(.Object,
                   gmt = data.frame(),
                   taxID = character(),
                   ontologyName = character(),
                   version = character(),
                   scientificName = character(),
                   commonEnglishName = character(),
                   description = character(),
                   # TODO : Implement.
                   modelBinaryMatrix = data.frame(),
                   ...) {

            .Object@gmt <- gmt
            .Object@taxID <- taxID
            .Object@ontologyName <- ontologyName
            .Object@version <- version
            .Object@scientificName <- scientificName
            .Object@commonEnglishName <- commonEnglishName
            .Object@description <- description

            # TODO : An algorithm to count it from clean model.
            .Object@modelBinaryMatrix <- data.frame(a = c(1,2), b = c("asa", "bsb"))

            .Object

          })
