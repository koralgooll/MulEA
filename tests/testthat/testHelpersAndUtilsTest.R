test_that("Methods : checkIfPoolIncludeSample false", {
  gmtMock <- data.frame(
    ontologyId = "GO:0000001",
    ontologyName = "Imagin gen ontology to tests.",
    listOfValues = I(list(c("a", "b", "c", "d"))),
    stringsAsFactors = FALSE
  )
  testDataMock <- c("a", "b", "e", "f")
  
  testthat::expect_warning(checkIfPoolIncludeSample(model = gmtMock, sampleVector = testDataMock))
})

test_that("Methods : checkIfPoolIncludeSample false", {
  gmtMock <- data.frame(
    ontologyId = "GO:0000001",
    ontologyName = "Imagin gen ontology to tests.",
    listOfValues = I(list(c("a", "b", "c", "d"))),
    stringsAsFactors = FALSE
  )
  testDataMock <- c("a", "b", "d")
  
  testthat::expect_equal(
    checkIfPoolIncludeSample(model = gmtMock, sampleVector = testDataMock),
    c("a", "b", "d")
  )
})


test_that("Methods : cutGmtToPool", {
  gmtMock <- data.frame(
    ontologyId = "GO:0000001",
    ontologyName = "Imagin gen ontology to tests.",
    listOfValues = I(list(c("a", "b", "c", "d"))),
    stringsAsFactors = FALSE
  )
  poolMock <- c("a", "b", "e", "f", "g", "h")
  
  testthat::expect_equal(cutGmtToPool(gmt = gmtMock, pool = poolMock)[['listOfValues']][[1]],
                         list(c("a", "b"))[[1]])
})
