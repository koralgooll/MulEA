context("ORA")

test_that("ORA : object creation test.", {
  gmtMock <- data.frame(ontologyId = "GO:0000001",
                        ontologyName = "Imagin gen ontology to tests.",
                        listOfValues = I(list(c("a", "b", "c"))),
                        stringsAsFactors = FALSE)
  testDataMock <- c("a", "b", "c")
  poolMock <- c("a", "c", "d")
  
  mulea_ora_model <- MulEA::ORA(
    gmt = gmtMock, testData = testDataMock, pool=poolMock,
    adjustMethod = "PT")
  
  testthat::expect_equal(mulea_ora_model@gmt, gmtMock)
  testthat::expect_equal(mulea_ora_model@testData, c("a", "b", "c"))
  testthat::expect_equal(mulea_ora_model@pool, c("a", "c", "d"))
})

test_that("ORA : testData out of DB model.", {
  gmtMock <- data.frame(ontologyId = "GO:0000001",
                        ontologyName = "Imagin gen ontology to tests.",
                        listOfValues = I(list(c("a", "b", "c"))),
                        stringsAsFactors = FALSE)
  testDataMock <- c("a", "b", "d")

  mulea_ora_model <- MulEA::ORA(gmt = gmtMock, testData = testDataMock, adjustMethod = "PT")
  
  testthat::expect_warning(muleaTestRes <- MulEA::runTest(mulea_ora_model))
  testthat::expect_equal(muleaTestRes$pValue, 1)
})

test_that("ORA : testData out of pool.", {
  gmtMock <- data.frame(ontologyId = "GO:0000001",
                        ontologyName = "Imagin gen ontology to tests.",
                        listOfValues = I(list(c("a", "b", "c"))),
                        stringsAsFactors = FALSE)
  testDataMock <- c("a", "b", "c")
  poolMock <- c("a", "b", "d")

  mulea_ora_model <- MulEA::ORA(
    gmt = gmtMock, testData = testDataMock, pool=poolMock,
    adjustMethod = "PT")
  
  testthat::expect_warning(muleaTestRes <- MulEA::runTest(mulea_ora_model))
  testthat::expect_equal(muleaTestRes$pValue, 1/3)
})

test_that("ORA : matrix 2,2,2,2.", {
  gmtMock <- data.frame(ontologyId = "GO:0000001",
                        ontologyName = "Imagin gen ontology to tests.",
                        listOfValues = I(list(c("a", "b", "c", "d"))),
                        stringsAsFactors = FALSE)
  testDataMock <- c("a", "b", "e", "f")
  poolMock <- c("a", "b", "c", "d", "e", "f", "g", "h")

  mulea_ora_model <- MulEA::ORA(
    gmt = gmtMock, testData = testDataMock, pool=poolMock,
    adjustMethod = "PT")
  
  muleaTestRes <- MulEA::runTest(mulea_ora_model)
  testthat::expect_equal(muleaTestRes$pValue, 53/70)
})

test_that("ORA : pool >> var + DBi, matrix 2,2,2,18.", {
  gmtMock <- data.frame(ontologyId = "GO:0000001",
                        ontologyName = "Imagin gen ontology to tests.",
                        listOfValues = I(list(c("a", "b", "c", "d"))),
                        stringsAsFactors = FALSE)
  testDataMock <- c("a", "b", "e", "f")
  poolMock <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p",
                "q", "r", "s", "t", "u", "w", "x", "y")

  mulea_ora_model <- MulEA::ORA(
    gmt = gmtMock, testData = testDataMock, pool=poolMock,
    adjustMethod = "PT")
  
  muleaTestRes <- MulEA::runTest(mulea_ora_model)
  testthat::expect_equal(muleaTestRes$pValue, 37/322)
})

test_that("ORA : DBi not include pool, matrix 2,0,2,2.", {
  gmtMock <- data.frame(ontologyId = "GO:0000001",
                        ontologyName = "Imagin gen ontology to tests.",
                        listOfValues = I(list(c("a", "b", "c", "d"))),
                        stringsAsFactors = FALSE)
  testDataMock <- c("a", "b", "e", "f")
  poolMock <- c("a", "b", "e", "f", "g", "h")

  mulea_ora_model <- MulEA::ORA(
    gmt = gmtMock, testData = testDataMock, pool=poolMock,
    adjustMethod = "PT")
  
  muleaTestRes <- MulEA::runTest(mulea_ora_model)
  testthat::expect_equal(muleaTestRes$pValue, 0.4)
})

test_that("ORA : DB1 + DB2 => pool, matrix 1,3,2,2 and 2,2,1,3.", {
  gmtMock1 <- data.frame(ontologyId = "GO:0000001",
                        ontologyName = "Imagin gen ontology to tests.",
                        listOfValues = I(list(c("a", "b", "c", "d"))),
                        stringsAsFactors = FALSE)
  gmtMock2 <- data.frame(ontologyId = "GO:0000002",
                        ontologyName = "Imagin gen ontology to tests.",
                        listOfValues = I(list(c("e", "f", "g", "h"))),
                        stringsAsFactors = FALSE)
  gmtMock <- rbind(gmtMock1, gmtMock2)
  testDataMock <- c("d", "e", "f")

  mulea_ora_model <- MulEA::ORA(gmt = gmtMock, 
    testData = testDataMock, adjustMethod = "PT")
  
  muleaTestRes <- MulEA::runTest(mulea_ora_model)
  testthat::expect_equal(muleaTestRes$pValue, c(13/14, 0.5))
})

test_that("ORA : DB1 + DB2 => pool, matrix 2,2,2,0 and 2,2,1,3.", {
  gmtMock1 <- data.frame(ontologyId = "GO:0000001",
                         ontologyName = "Imagin gen ontology to tests.",
                         listOfValues = I(list(c("a", "b", "c", "d"))),
                         stringsAsFactors = FALSE)
  gmtMock2 <- data.frame(ontologyId = "GO:0000002",
                         ontologyName = "Imagin gen ontology to tests.",
                         listOfValues = I(list(c("e", "f", "c", "d"))),
                         stringsAsFactors = FALSE)
  gmtMock <- rbind(gmtMock1, gmtMock2)
  testDataMock <- c("b", "d", "e", "f")
  poolMock <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p",
                "q", "r", "s", "t", "u", "w", "x", "y")

  mulea_ora_model <- MulEA::ORA(
    gmt = gmtMock, testData = testDataMock, pool=poolMock,
    adjustMethod = "PT")
  
  muleaTestRes <- MulEA::runTest(mulea_ora_model)
  testthat::expect_equal(muleaTestRes$pValue, c(37/322, 27/3542))
})
