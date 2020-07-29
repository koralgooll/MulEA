# Eszter's first example.
rm(list = ls())
root_of_cpp_project <- 'D:/projects/science/enrichment-analysis'
setwd(root_of_cpp_project)


convert_model_list_to_model_df <- function(model_list) {
  model_df <- plyr::adply(.data = names(DB), .margins=1, .fun = function(list_name){
    
    data.frame('ontologyId' = list_name, 'ontologyName' = list_name, 
               'listOfValues' = I(DB[list_name]), stringsAsFactors = FALSE)
  })[c('ontologyId', 'ontologyName', 'listOfValues')]
  model_df
}


# Test Case 01
source("src/set-based-enrichment-test.R") # It contains the R source of enrichment analizis and it compiles the C++ part.
load(file="test-parameters-01.RData") # It loas a dummy databaseinto variabla 'DB', and a 'pool' and 'selet' lists; The DB contains only disjoint sets.  One of the sets is empty.
number_of_steps <- 100

your_res_01 <- set.based.enrichment.test(steps=number_of_steps, pool=pool, 
                                         select=select, DB=DB)

model_df <- convert_model_list_to_model_df(DB)
mulea_ora_01 <- MulEA::ORA(gmt = model_df, testData = select,
                           pool = pool, adjustMethod = "PT",
                           numberOfPermutations = number_of_steps)
mulea_res_01 <- MulEA::runTest(mulea_ora_01)



# Plot Processing.

library(ggplot2)
library(dplyr)

selection_vector = c(1, 7, 15, 24, 27, 51, 61, 84, 86)
mulea_res_01[selection_vector,] %>% ggplot( aes(x=DB_names, y=FDR, fill=FDR)) +
  geom_bar(stat="identity") +
  scale_fill_gradient2(mid='red', high='darkgreen', space='Lab') +
  coord_flip() +
  xlab("") +
  theme_bw()

library(reshape2)

mulea_res_01.long <- melt(mulea_res_01[selection_vector,])
mulea_res_01.long$DB_names <- as.character(mulea_res_01.long$DB_names)
mulea_res_01.long$variable <- as.character(mulea_res_01.long$variable)

mulea_res_01.long %>% filter(variable %in% c('P', 'FDR')) %>% 
  ggplot(aes(DB_names, value, alpha=variable, fill=value)) +
  geom_bar(stat="identity", position = "dodge") +
  scale_alpha_ordinal(range = c(0.5, 0.9)) +
  scale_fill_gradient2(mid='red', high='darkgreen', space='Lab') +
  coord_flip() + 
  xlab("") +
  theme_bw()



mulea_res_01.long %>% filter(variable %in% c('P', 'FDR')) %>% 
  ggplot(aes(DB_names, value, alpha=variable, fill=value)) +
  geom_bar(stat="identity", position = "dodge") +
  scale_alpha_ordinal(range = c(0.5, 0.9)) +
  scale_fill_gradient2(mid='red', high='darkgreen', space='Lab') +
  coord_flip() + 
  xlab("") +
  theme_bw()












