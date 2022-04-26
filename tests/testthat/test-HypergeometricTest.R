
library(MulEA)
context("MuleaHypergeometricTest")

test_that("MuleaHypergeometricTest : object creation test.", {
  gmtMock <- data.frame(ontologyId = "GO:0000001",
                        ontologyName = "Imagin gen ontology to tests.",
                        listOfValues = I(list(c("a", "b", "c"))),
                        stringsAsFactors = FALSE)
  testDataMock <- c("a", "b", "c")
  poolMock <- c("a", "c", "d")

  muleaHypergeometricTest <- MuleaHypergeometricTest(gmt = gmtMock, testData = testDataMock, pool = poolMock)

  expect_equal(muleaHypergeometricTest@gmt, gmtMock)
  expect_equal(muleaHypergeometricTest@testData, c("a", "b", "c"))
  expect_equal(muleaHypergeometricTest@pool, c("a", "c", "d"))
})

test_that("MuleaHypergeometricTest : testData out of DB model.", {
  gmtMock <- data.frame(ontologyId = "GO:0000001",
                        ontologyName = "Imagin gen ontology to tests.",
                        listOfValues = I(list(c("a", "b", "c"))),
                        stringsAsFactors = FALSE)
  testDataMock <- c("a", "b", "d")

  muleaHypergeometricTest <- MuleaHypergeometricTest(gmt = gmtMock, testData = testDataMock)
  suppressWarnings(muleaHypergeometricTestRes <- runTest(muleaHypergeometricTest))

  expect_equal(muleaHypergeometricTestRes$p.value,
               fisher.test(matrix(c(2, 1, 0, 0), 2, 2), alternative = "greater")$p.value)
})

test_that("MuleaHypergeometricTest : testData out of pool.", {
  gmtMock <- data.frame(ontologyId = "GO:0000001",
                        ontologyName = "Imagin gen ontology to tests.",
                        listOfValues = I(list(c("a", "b", "c"))),
                        stringsAsFactors = FALSE)
  testDataMock <- c("a", "b", "c")
  poolMock <- c("a", "b", "d")

  muleaHypergeometricTest <- MuleaHypergeometricTest(gmt = gmtMock, testData = testDataMock, pool = poolMock)


  expect_warning(muleaHypergeometricTestRes <- runTest(muleaHypergeometricTest))
  expect_equal(muleaHypergeometricTestRes$p.value,
              fisher.test(matrix(c(2, 0, 0, 1), 2, 2), alternative = "greater")$p.value)
})

test_that("MuleaHypergeometricTest : matrix 2,2,2,2.", {
  gmtMock <- data.frame(ontologyId = "GO:0000001",
                        ontologyName = "Imagin gen ontology to tests.",
                        listOfValues = I(list(c("a", "b", "c", "d"))),
                        stringsAsFactors = FALSE)
  testDataMock <- c("a", "b", "e", "f")
  poolMock <- c("a", "b", "c", "d", "e", "f", "g", "h")

  muleaHypergeometricTest <- MuleaHypergeometricTest(gmt = gmtMock, testData = testDataMock, pool = poolMock)
  muleaHypergeometricTestRes <- runTest(muleaHypergeometricTest)

  expect_equal(muleaHypergeometricTestRes$p.value,
               fisher.test(matrix(c(2, 2, 2, 2), 2, 2), alternative = "greater")$p.value)
})

test_that("MuleaHypergeometricTest : pool >> var + DBi, matrix 2,2,2,18.", {
  gmtMock <- data.frame(ontologyId = "GO:0000001",
                        ontologyName = "Imagin gen ontology to tests.",
                        listOfValues = I(list(c("a", "b", "c", "d"))),
                        stringsAsFactors = FALSE)
  testDataMock <- c("a", "b", "e", "f")
  poolMock <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p",
                "q", "r", "s", "t", "u", "w", "x", "y")

  muleaHypergeometricTest <- MuleaHypergeometricTest(gmt = gmtMock, testData = testDataMock, pool = poolMock)
  muleaHypergeometricTestRes <- runTest(muleaHypergeometricTest)

  expect_equal(muleaHypergeometricTestRes$p.value,
               fisher.test(matrix(c(2, 2, 2, 18), 2, 2), alternative = "greater")$p.value)
})

test_that("MuleaHypergeometricTest : DBi not include pool, matrix 2,0,2,2.", {
  gmtMock <- data.frame(ontologyId = "GO:0000001",
                        ontologyName = "Imagin gen ontology to tests.",
                        listOfValues = I(list(c("a", "b", "c", "d"))),
                        stringsAsFactors = FALSE)
  testDataMock <- c("a", "b", "e", "f")
  poolMock <- c("a", "b", "e", "f", "g", "h")

  muleaHypergeometricTest <- MuleaHypergeometricTest(gmt = gmtMock, testData = testDataMock, pool = poolMock)
  muleaHypergeometricTestRes <- runTest(muleaHypergeometricTest)

  expect_equal(muleaHypergeometricTestRes$p.value,
               fisher.test(matrix(c(2, 0, 2, 2), 2, 2), alternative = "greater")$p.value)
})

test_that("MuleaHypergeometricTest : DB1 + DB2 => pool, matrix 1,3,2,2 and 2,2,1,3.", {
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

  muleaHypergeometricTest <- MuleaHypergeometricTest(gmt = gmtMock, testData = testDataMock)
  muleaHypergeometricTestRes <- runTest(muleaHypergeometricTest)

  expect_equal(muleaHypergeometricTestRes$p.value,
               c(fisher.test(matrix(c(1, 3, 2, 2), 2, 2), alternative = "greater")$p.value,
                 fisher.test(matrix(c(2, 2, 1, 3), 2, 2), alternative = "greater")$p.value))
})

test_that("MuleaHypergeometricTest : DB1 + DB2 => pool, matrix 2,2,2,0 and 2,2,1,3.", {
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

  muleaHypergeometricTest <- MuleaHypergeometricTest(gmt = gmtMock, testData = testDataMock)
  muleaHypergeometricTestRes <- runTest(muleaHypergeometricTest)

  expect_equal(muleaHypergeometricTestRes$p.value,
               c(fisher.test(matrix(c(2, 2, 2, 0), 2, 2), alternative = "greater")$p.value,
                 fisher.test(matrix(c(3, 1, 1, 1), 2, 2), alternative = "greater")$p.value))
})
