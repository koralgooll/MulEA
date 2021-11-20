#testing the three MulEA scripts a hundred times

#setwd('~/Documents/MulEA_testing/')
#Database

#Creating a list of DB IDs and associated genes from text file

library(tidyverse)
library(tictoc)

N1<-40 # ennyi GO kategoria lesz
N2<-30 # minden kategoriaban ennyi gen lesz

# olyan teszt adatbazist keszitek, ahol minden ketegoriaban ugyanannyi gen van, es nincsenek atfedesek
DB<-list()
for(i1 in 1:N1)
{
  DB[[sprintf("CAT_%04i_%04i",(i1-1)*N2+1,i1*N2)]] <-sprintf("gene_%04i",(i1-1)*N2+1:N2)
}

list_of_all_genes<-unlist(DB)
names(list_of_all_genes)<-NULL

# ------------------------- CPP based test -------------------------

DB_df <- MulEA:::convertListToGmtDataFrame(ontologyReprAsList = DB)
number_of_steps <- 5000


tic()
mulea_ora_model <- MulEA::ORA(
  gmt = DB_df, testData = DB_df$listOfValues[[3]], adjustMethod = "PT",
  numberOfPermutations = number_of_steps, nthreads = 16, pool = list_of_all_genes)
result_tbl <- MulEA::runTest(mulea_ora_model)
toc()


tic()
select <- unlist(DB_df$listOfValues[3:4])
mulea_ora_model <- MulEA::ORA(
  gmt = DB_df, testData = select, adjustMethod = "PT",
  numberOfPermutations = number_of_steps, nthreads = 16, pool = list_of_all_genes)
result_tbl <- MulEA::runTest(mulea_ora_model)
toc()


tic()
select <- unlist(DB_df$listOfValues[1:6])
mulea_ora_model <- MulEA::ORA(
  gmt = DB_df, testData = select, adjustMethod = "PT",
  numberOfPermutations = number_of_steps, nthreads = 16, pool = list_of_all_genes)
result_tbl <- MulEA::runTest(mulea_ora_model)
toc()


tic()
select<-list_of_all_genes
mulea_ora_model <- MulEA::ORA(
  gmt = DB_df, testData = select, adjustMethod = "PT",
  numberOfPermutations = number_of_steps, nthreads = 16, pool = list_of_all_genes)
result_tbl <- MulEA::runTest(mulea_ora_model)
toc()


tic()
select <- unlist(DB_df$listOfValues[1:30])
mulea_ora_model <- MulEA::ORA(
  gmt = DB_df, testData = select, adjustMethod = "PT",
  numberOfPermutations = number_of_steps, nthreads = 16, pool = list_of_all_genes)
result_tbl <- MulEA::runTest(mulea_ora_model)
toc()



tic()
list1<-list()
for( i in 1:800)
{
  # FIX : sort will slow down. sort shouldn't impact stats tests.
  # select <- sort(sample(list_of_all_genes, size = 100, replace = FALSE))
  select <- sample(list_of_all_genes, size = 100, replace = FALSE)
  mulea_ora_model <- MulEA::ORA(
    gmt = DB_df, testData = select, adjustMethod = "PT",
    numberOfPermutations = number_of_steps, nthreads = 16, pool = list_of_all_genes)
  result_tbl <- MulEA::runTest(mulea_ora_model)
  list1[[i]] <- result_tbl
}
toc()
# write_rds(list1, "test-results-of-random-samples-001.rds")
# list1<-read_rds( "test-results-of-random-samples-001.rds")


# Tests by plotting, no enrichment at all.
# Clear p-values
x <- sapply(list1, function(r){r$pValue[[7]]})
x <- as.numeric(unlist(sapply(list1, function(r){r$pValue})))
tibble(P_val=x) %>% ggplot(mapping = aes(x=P_val))+stat_ecdf()+geom_abline(color="green")

# Our permutation test approach
x <- sapply(list1, function(r){r$adjustedPValueEmpirical[[7]]})
# QUESTION : What is the idea behind this cutoff? 
x <- as.numeric(unlist(sapply(list1, function(r){r$adjustedPValueEmpirical[[7]]<0.15})))
x <- as.numeric(unlist(sapply(list1, function(r){r$adjustedPValueEmpirical[[7]]})))
tibble(P_val=x) %>% ggplot(mapping = aes(x=P_val))+stat_ecdf()+geom_abline(color="green")

# BH adjustment method
x <- sapply(list1, function(r){r$adjustedPValue[[7]]})
x <- as.numeric(unlist(sapply(list1, function(r){r$adjustedPValue[[7]]<0.15})))
x <- as.numeric(unlist(sapply(list1, function(r){r$adjustedPValue})))
tibble(P_val=x) %>% ggplot(mapping = aes(x=P_val))+stat_ecdf()+geom_abline(color="green")


really_enriched <- 1:4
# gene_list_in <- sort(unique(unlist(DB[really_enriched])))
gene_list_in <- unique(unlist(DB_df$listOfValues[really_enriched]))
gene_list_out <- setdiff(list_of_all_genes, gene_list_in)
list1<-list()
tic()
for( i in 1:1500)
{
  select_in <- sample(gene_list_in, size = 10, replace = FALSE)
  select_out <- sample(gene_list_out, size = 40, replace = FALSE)
  select <- c(select_in, select_out)
  mulea_ora_model <- MulEA::ORA(
    gmt = DB_df, testData = select, adjustMethod = "PT",
    numberOfPermutations = number_of_steps, nthreads = 16, pool = list_of_all_genes)
  result_tbl <- MulEA::runTest(mulea_ora_model)
  list1[[i]] <- result_tbl
}
toc()


# write_rds(list1, "test-results-with-enr-001-004.rds")
# list1<-read_rds( "test-results-of-random-samples-001.rds")

result_tbl2<-tibble(FDR=rep(as.numeric(NA),length(list1)))
for( i in seq(list1))
{
  result_tbl <- list1[[i]]
  result_tbl <- result_tbl %>%
    mutate(significant=adjustedPValueEmpirical<0.05) %>%
    mutate(really_true=1:n() %in% really_enriched)  %>% 
    mutate(c=case_when(
      significant & really_true ~ "TP",
      !significant & !really_true ~ "TN",
      significant & !really_true ~ "FP",
      !significant & really_true ~ "FN",
      TRUE~"ERROR_346"
      ))
  
  TP<-result_tbl %>% filter(    significant & really_true) %>%  nrow()
  TN<-result_tbl %>% filter(    !significant & !really_true) %>%  nrow()
  FP<-result_tbl %>% filter(    significant & !really_true) %>%  nrow()
  FN<-result_tbl %>% filter(    !significant & really_true) %>%  nrow()
  
  result_tbl2$FDR[[i]]<-ifelse((TP+FP)!=0,  FP/(TP+FP) , 0)

}

result_tbl2 %>% ggplot( mapping = aes(x=FDR)) +  geom_histogram()+geom_rug()
mean(result_tbl2$FDR, na.rm = TRUE)
