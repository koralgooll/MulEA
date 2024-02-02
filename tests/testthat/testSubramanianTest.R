test_that("GSEA : object creation test.", {
  gmtMock <- data.frame(
    ontologyId = "GO:0000001",
    ontologyName = "Imagin gen ontology to tests.",
    listOfValues = I(list(c("a", "b", "c"))),
    stringsAsFactors = FALSE
  )
  testDataMock <- c("a", "b", "c")
  scoreDataMock <- c(0.1, 0.5, 1)
  
  mulea_ranked_based_test_model <- mulea::gsea(gmt = gmtMock,
                                               element_names = testDataMock,
                                               element_scores = scoreDataMock)
  
  testthat::expect_equal(mulea_ranked_based_test_model@gmt, gmtMock)
  testthat::expect_equal(mulea_ranked_based_test_model@element_names, c("a", "b", "c"))
  testthat::expect_equal(mulea_ranked_based_test_model@element_scores, c(0.1, 0.5, 1))
})

test_that("GSEA : no element_scores vector.", {
  gmtMock <- data.frame(
    ontologyId = "GO:0000001",
    ontologyName = "Imagin gen ontology to tests.",
    listOfValues = I(list(c("a", "b", "c"))),
    stringsAsFactors = FALSE
  )
  testDataMock <- c("a", "b", "d")
  
  mulea_ranked_based_test_model <- mulea::gsea(gmt = gmtMock, element_names = testDataMock)
  testthat::expect_error(muleaTestRes <-
                           mulea::run_test(mulea_ranked_based_test_model))
})

test_that("GSEA : out of ontology elements.", {
  set.seed(1)
  gmtMock <- data.frame(
    ontologyId = "GO:0000001",
    ontologyName = "Imagin gen ontology to tests.",
    listOfValues = I(list(c("a", "b", "c", "d"))),
    stringsAsFactors = FALSE
  )
  testDataMock <- c("a", "b", "e", "f")
  scoreDataMock <- c(0.1, 0.5, 0.7, 1)
  
  mulea_ranked_based_test_model <- mulea::gsea(gmt = gmtMock,
                                               element_names = testDataMock,
                                               element_scores = scoreDataMock)
  
  testthat::expect_warning(muleaTestRes <-
                             mulea::run_test(mulea_ranked_based_test_model))
  testthat::expect_equal(muleaTestRes$pValue, 151/495)
})

test_that("GSEA : DB1 + DB2.", {
  set.seed(1)
  gmtMock1 <- data.frame(
    ontologyId = "GO:0000001",
    ontologyName = "Imagin gen ontology to tests.",
    listOfValues = I(list(c("a", "b", "c", "d"))),
    stringsAsFactors = FALSE
  )
  gmtMock2 <- data.frame(
    ontologyId = "GO:0000002",
    ontologyName = "Imagin gen ontology to tests.",
    listOfValues = I(list(c("e", "f", "g", "h"))),
    stringsAsFactors = FALSE
  )
  gmtMock <- rbind(gmtMock1, gmtMock2)
  testDataMock <- c("d", "e", "f")
  scoreDataMock <- c(-0.3, 0.4, 0.5)
  
  mulea_ranked_based_test_model <- mulea::gsea(gmt = gmtMock,
                                               element_names = testDataMock,
                                               element_scores = scoreDataMock)
  
  muleaTestRes <- mulea::run_test(mulea_ranked_based_test_model)
  testthat::expect_equal(muleaTestRes$pValue, c(169 / 330, 107 / 219))
})

test_that("GSEA : DB1 + DB2 out of background_element_names.", {
  set.seed(1)
  gmtMock1 <- data.frame(
    ontologyId = "GO:0000001",
    ontologyName = "Imagin gen ontology to tests.",
    listOfValues = I(list(c("a", "b", "c", "d"))),
    stringsAsFactors = FALSE
  )
  gmtMock2 <- data.frame(
    ontologyId = "GO:0000002",
    ontologyName = "Imagin gen ontology to tests.",
    listOfValues = I(list(c("e", "f", "c", "d"))),
    stringsAsFactors = FALSE
  )
  gmtMock <- rbind(gmtMock1, gmtMock2)
  testDataMock <- c("b", "d", "e", "f")
  scoreDataMock <- c(-1, 0, 0, 1)
  
  mulea_ranked_based_test_model <- mulea::gsea(gmt = gmtMock,
                                               element_names = testDataMock,
                                               element_scores = scoreDataMock)
  
  muleaTestRes <- mulea::run_test(mulea_ranked_based_test_model)
  testthat::expect_equal(muleaTestRes$pValue, c(156 / 653, 251 / 492))
})
