
library(MulEA)
context("KolmogorovSmirnovTest")

test_that("KolmogorovSmirnovTest : object creation test (pool is testData in this test), (score is testData to).", {
  gmtMock <- data.frame(ontologyId = "GO:0000001",
                        ontologyName = "Imagin gen ontology to tests.",
                        listOfValues = I(list(c("a", "b", "c"))),
                        stringsAsFactors = FALSE)
  testDataMock <- c("a", "b", "c")

  ksTest <- KolmogorovSmirnovTest(gmt = gmtMock, testData = testDataMock)

  expect_equal(ksTest@gmt, gmtMock)
  expect_equal(ksTest@testData, c("a", "b", "c"))
})

test_that("KolmogorovSmirnovTest : DB model partialy out of testData.", {
  gmtMock <- data.frame(ontologyId = "GO:0000001",
                        ontologyName = "Imagin gen ontology to tests.",
                        listOfValues = I(list(c("a", "b", "c"))),
                        stringsAsFactors = FALSE)
  testDataMock <- c("a", "b", "d")

  ksTest <- KolmogorovSmirnovTest(gmt = gmtMock, testData = testDataMock)

  testthat::expect_warning(rankedTestRes <- runTest(ksTest))
  testthat::expect_gte(object = rankedTestRes$p.value, expected = 0.974)
  testthat::expect_lte(object = rankedTestRes$p.value, expected = 0.978)
})

test_that("KolmogorovSmirnovTest : empty intersection of experiment data with category value.", {
  gmtMock <- data.frame(ontologyId = "GO:0000001",
                        ontologyName = "Imagin gen ontology to tests.",
                        listOfValues = I(list(c("a", "b", "c"))),
                        stringsAsFactors = FALSE)
  testDataMock <- c("d", "e", "f", "g")

  ksTest <- KolmogorovSmirnovTest(gmt = gmtMock, testData = testDataMock)
  rankedTestRes <- runTest(ksTest)

  expect_equal(rankedTestRes$p.value, NA)
})

test_that("KolmogorovSmirnovTest : DB model fully inside testData.", {
  gmtMock <- data.frame(ontologyId = "GO:0000001",
                        ontologyName = "Imagin gen ontology to tests.",
                        listOfValues = I(list(c("c", "f", "j"))),
                        stringsAsFactors = FALSE)
  testDataMock <- c("a", "b", "c", "d", "e", "f", "g", "h", "j")

  ksTest <- KolmogorovSmirnovTest(gmt = gmtMock, testData = testDataMock)
  suppressWarnings(rankedTestRes <- runTest(ksTest))

  testthat::expect_gte(object = rankedTestRes$p.value, expected = 0.898)
  testthat::expect_lte(object = rankedTestRes$p.value, expected = 0.934)

})

test_that("KolmogorovSmirnovTest : DB model the same as testData.", {
  gmtMock <- data.frame(ontologyId = "GO:0000001",
                        ontologyName = "Imagin gen ontology to tests.",
                        listOfValues = I(list(c("a", "b", "c", "g", "h", "j"))),
                        stringsAsFactors = FALSE)
  testDataMock <- c("a", "b", "c", "g", "h", "j")

  ksTest <- KolmogorovSmirnovTest(gmt = gmtMock, testData = testDataMock)
  suppressWarnings(rankedTestRes <- runTest(ksTest))

  expect_equal(rankedTestRes$p.value, 1)
})
