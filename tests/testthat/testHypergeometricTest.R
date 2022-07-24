library(MulEA)
context("ORA")

test_that("ORA : object creation test.", {
  gmtMock <- data.frame(
    ontologyId = "GO:0000001",
    ontologyName = "Imagin gen ontology to tests.",
    listOfValues = I(list(c("a", "b", "c"))),
    stringsAsFactors = FALSE
  )
  testDataMock <- c("a", "b", "c")
  poolMock <- c("a", "c", "d")
  
  mulea_ora_model <- MulEA::ORA(
    gmt = gmtMock,
    testData = testDataMock,
    pool = poolMock,
    adjustMethod = "PT", 
    nthreads = 2)
  
  testthat::expect_equal(mulea_ora_model@gmt, gmtMock)
  testthat::expect_equal(mulea_ora_model@testData, c("a", "b", "c"))
  testthat::expect_equal(mulea_ora_model@pool, c("a", "c", "d"))
})

test_that("ORA : testData out of DB model.", {
  gmtMock <- data.frame(
    ontologyId = "GO:0000001",
    ontologyName = "Imagin gen ontology to tests.",
    listOfValues = I(list(c("a", "b", "c"))),
    stringsAsFactors = FALSE
  )
  testDataMock <- c("a", "b", "d")
  mulea_ora_model <-
    MulEA::ORA(gmt = gmtMock,
               testData = testDataMock,
               adjustMethod = "PT", 
               nthreads = 2)
  
  testthat::expect_warning(muleaTestRes <-
                             MulEA::run_test(mulea_ora_model))
  testthat::expect_equal(muleaTestRes$pValue, 1)
})

test_that("ORA : testData out of pool.", {
  gmtMock <- data.frame(
    ontologyId = "GO:0000001",
    ontologyName = "Imagin gen ontology to tests.",
    listOfValues = I(list(c("a", "b", "c"))),
    stringsAsFactors = FALSE
  )
  testDataMock <- c("a", "b", "c")
  poolMock <- c("a", "b", "d")
  
  mulea_ora_model <- MulEA::ORA(
    gmt = gmtMock,
    testData = testDataMock,
    pool = poolMock,
    adjustMethod = "PT", 
    nthreads = 2)
  
  testthat::expect_warning(muleaTestRes <-
                             MulEA::run_test(mulea_ora_model))
  testthat::expect_equal(muleaTestRes$pValue, 1 / 3)
})

test_that("ORA : matrix 2,2,2,2.", {
  gmtMock <- data.frame(
    ontologyId = "GO:0000001",
    ontologyName = "Imagin gen ontology to tests.",
    listOfValues = I(list(c("a", "b", "c", "d"))),
    stringsAsFactors = FALSE
  )
  testDataMock <- c("a", "b", "e", "f")
  poolMock <- c("a", "b", "c", "d", "e", "f", "g", "h")
  
  mulea_ora_model <- MulEA::ORA(
    gmt = gmtMock,
    testData = testDataMock,
    pool = poolMock,
    adjustMethod = "PT", 
    nthreads = 2)
  
  muleaTestRes <- MulEA::run_test(mulea_ora_model)
  testthat::expect_equal(muleaTestRes$pValue, 53 / 70)
})

test_that("ORA : pool >> var + DBi, matrix 2,2,2,18.", {
  gmtMock <- data.frame(
    ontologyId = "GO:0000001",
    ontologyName = "Imagin gen ontology to tests.",
    listOfValues = I(list(c("a", "b", "c", "d"))),
    stringsAsFactors = FALSE
  )
  testDataMock <- c("a", "b", "e", "f")
  poolMock <-
    c(
      "a",
      "b",
      "c",
      "d",
      "e",
      "f",
      "g",
      "h",
      "i",
      "j",
      "k",
      "l",
      "m",
      "n",
      "o",
      "p",
      "q",
      "r",
      "s",
      "t",
      "u",
      "w",
      "x",
      "y"
    )
  
  mulea_ora_model <- MulEA::ORA(
    gmt = gmtMock,
    testData = testDataMock,
    pool = poolMock,
    adjustMethod = "PT", 
    nthreads = 2)
  
  muleaTestRes <- MulEA::run_test(mulea_ora_model)
  testthat::expect_equal(muleaTestRes$pValue, 37 / 322)
})

test_that("ORA : DBi not include pool, matrix 2,0,2,2.", {
  gmtMock <- data.frame(
    ontologyId = "GO:0000001",
    ontologyName = "Imagin gen ontology to tests.",
    listOfValues = I(list(c("a", "b", "c", "d"))),
    stringsAsFactors = FALSE
  )
  testDataMock <- c("a", "b", "e", "f")
  poolMock <- c("a", "b", "e", "f", "g", "h")
  
  mulea_ora_model <- MulEA::ORA(
    gmt = gmtMock,
    testData = testDataMock,
    pool = poolMock,
    adjustMethod = "PT", 
    nthreads = 2)
  
  muleaTestRes <- MulEA::run_test(mulea_ora_model)
  testthat::expect_equal(muleaTestRes$pValue, 0.4)
})

test_that("ORA : DB1 + DB2 => pool, matrix 1,3,2,2 and 2,2,1,3.", {
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

  mulea_ora_model <- MulEA::ORA(gmt = gmtMock,
                                testData = testDataMock,
                                adjustMethod = "PT", 
                                nthreads = 2)
  
  muleaTestRes <- MulEA::run_test(mulea_ora_model)
  testthat::expect_equal(muleaTestRes$pValue, c(13 / 14, 0.5))
})

test_that("ORA : DB1 + DB2 => pool, matrix 2,2,2,0 and 2,2,1,3.", {
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
  poolMock <-
    c(
      "a",
      "b",
      "c",
      "d",
      "e",
      "f",
      "g",
      "h",
      "i",
      "j",
      "k",
      "l",
      "m",
      "n",
      "o",
      "p",
      "q",
      "r",
      "s",
      "t",
      "u",
      "w",
      "x",
      "y"
    )
  
  mulea_ora_model <- MulEA::ORA(
    gmt = gmtMock,
    testData = testDataMock,
    pool = poolMock,
    adjustMethod = "PT", 
    nthreads = 2)
  
  muleaTestRes <- MulEA::run_test(mulea_ora_model)
  testthat::expect_equal(muleaTestRes$pValue, c(37 / 322, 27 / 3542))
})


test_that("ORA : DB1 + DB2 => pool, matrix 2,2,2,0 and 2,2,1,3.", {
  gmtMock <- MulEA::readGmtFileAsDataFrame(gmtFilePath = system.file(package="MulEA", "extdata", "model.gmt"))
  testDataMock <- c("FBgn0004407", "FBgn0010438", "FBgn0037044", "FBgn0002887", "FBgn0028434", "FBgn0030170", "FBgn0263831")
  poolMock <- unique(c(c("FBgn0033690", "FBgn0261618", "FBgn0004407", "FBgn0010438", "FBgn0032154", "FBgn0039930", "FBgn0040268", "FBgn0013674",
                         "FBgn0037008", "FBgn0003116", "FBgn0037743", "FBgn0035401", "FBgn0037044", "FBgn0051005", "FBgn0026737", "FBgn0026751",
                         "FBgn0038704", "FBgn0002887", "FBgn0028434", "FBgn0030170", "FBgn0263831", "FBgn0000579"),
                       c("FBgn0066666", "FBgn0000000", "FBgn0099999", "FBgn0011111", "FBgn0022222", "FBgn0777777", "FBgn0333333")))
  
  mulea_ora_model <- MulEA::ORA(
    gmt = gmtMock, testData = testDataMock, pool=poolMock, nthreads = 2)
  
  muleaTestRes <- MulEA::run_test(mulea_ora_model)
  testthat::expect_equal(muleaTestRes$pValue[c(1,2)], c(1, 19205/26423))
})

