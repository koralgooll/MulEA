library(MulEA)
context("set.based.enrichment.test")

create_random_db <- function() {
  DB <- list()
  for (cat_i in 1:10) {
    DB_cat_values <- c()
    for (el_i in 1:floor(runif(1, min=2, max=10))) {
      DB_cat_values <- append(DB_cat_values, paste0("el_", el_i))
    }
    cat_id <- paste0("cat_", cat_i)
    DB[[cat_id]]  <- DB_cat_values
  }
  return(DB)
}

test_that("set.based.enrichment.test.", {
  steps <- 100
  DB <- create_random_db()
  pool <- unique(unlist(DB))
  select <- c("el_4", "el_5", "el_6")
  nthread <- 2
  res <- MulEA:::set.based.enrichment.test(
    steps = steps,
    pool = pool,
    select = select,
    DB = DB,
    nthread = nthread)
  
  testthat::expect_gt(res[["FDR"]][1], 0)
})
