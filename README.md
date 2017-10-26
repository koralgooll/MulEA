
# MulEA - A Tool for Multi Enrichment Analysis

Functional interpretation of the biological data typically involves identifying key genes, molecules, reactions or pathways by finding non-random changes between two or more conditions or phenotypes. Performing enrichment analysis on set of molecules selected from  differential omics analysis is a method of choice. Among many packages that can be applied for this task, only few could be applied either to multiple species, ontology types or providing an access to multiple statistics.

MulEA is addressing this gap in addition providing improved way to calculate correction for multiple testing that assume partial dependence between ontology terms and in result limits number of correct associations falsely scored as insignificant. Besides the commonly applied tests, MulEA provides a unique permutation based, empirical false discovery rate correction of the p-values to substitute the too conservative Bonferroni and Benjamini-Hochberg procedures.

MulEA allows enrichment analysis using most popular gene and pathway ontologies (GO, KEGG, Reactome). In addition, one can test enrichment in genomic locations and in gene expression, protein domain, miRNA and transcription factors data bases, all created from publicly available resources and presented in standardized manner. Beyond genes or proteins, MulEA even allows working with basically any kind of data types, i.e. small molecules, chromosome region, enhancers, molecular interactions or any other information defined by the user.

Mulea currently supports 25 organisms from bacteria to human. Because, in addition to knowledge-bases provided alongside the package, the user may provide its own ontology files, MulEA can work with any biological species.

To analyse the data MulEA provide multiple types of statistics in one tool, which allows the user to calculate over-representations using the hypergeometric test, and enrichment analyses of ranked input by the Kolmogorov-Smirnov test.
                   
To conclude, MulEA is a comprehensive enrichment software that allows expansive analyses using diverse ontologies, statistical models and p-value correction procedures that can extend our understanding of the results of various high-throughput analyses and therefore expand our knowledge.

An R-package for fast analysis of bioligical data. The package implements three different approaches of this type of analysis. This file include blueprint of package possibilities, to see more ... 


## Installation

```{r}
```


## Run

Loading libraries

```{r}
library(MulEA)
```

Loading example input data:

```{r}
muleaPkgDir <- find.package("MulEA")
modelDfFromFile <- MulEA::readGmtFileAsDataFrame(gmtFilePath = paste(muleaPkgDir,"/example/model.gmt", sep = ""))
dataFromExperiment <- c("FBgn0004407", "FBgn0010438", "FBgn0003742", "FBgn0029709", "FBgn0030341", "FBgn0037044", "FBgn0002887", "FBgn0028434", "FBgn0030170", "FBgn0263831", "FBgn0261618", "FBgn0038704", "FBgn0000579")
dataFromExperimentScores <- c(0.09, 0.11, 0.15, 0.20, 0.21, 0.24, 0.28, 0.30, 0.45, 0.50, 0.53, 0.60, 0.61)
dataFromExperimentPool <- unique(c(c("FBgn0033690", "FBgn0261618", "FBgn0004407", "FBgn0010438", "FBgn0032154", "FBgn0039930", "FBgn0040268", "FBgn0013674",
                                   "FBgn0037008", "FBgn0003116", "FBgn0037743", "FBgn0035401", "FBgn0037044", "FBgn0051005", "FBgn0026737", "FBgn0026751",
                                   "FBgn0038704", "FBgn0002887", "FBgn0028434", "FBgn0030170", "FBgn0263831", "FBgn0000579"),
                                 c("FBgn0066666", "FBgn0000000", "FBgn0099999", "FBgn0011111", "FBgn0022222", "FBgn0777777", "FBgn0333333", "FBgn0003742",
                                   "FBgn0029709", "FBgn0030341")))
```

Running MulEA (two implemented approaches):

- Set based tests represented by `SetBasedTest` class:

```{r}
setBasedTestWithPoolAndAdjust <- SetBasedTest(gmt = modelDfFromFile, testData = dataFromExperiment, pool = dataFromExperimentPool, adjustMethod = "BH")
setBasedTestWithPoolAndAdjustRes <- MulEA::runTest(setBasedTestWithPoolAndAdjust)
```

- Ranked based tests represented by `RankedBasedTest` class:

```{r}
rankedBasedTestSubramanian <- RankedBasedTest(method = "Subramanian", gmt = modelDfFromFile, testData = dataFromExperiment, scores = dataFromExperimentScores)
rankedBasedTestSubramanianRes <- MulEA::runTest(rankedBasedTestSubramanian)
```


## Results

Example of results in data.frame form:

```
|ontologyId |ontologyName                                         |listOfValues                                                                                           |   p.value|  
|:----------|:----------------------------------------------------|:------------------------------------------------------------------------------------------------------|---------:|  
|ID:0000001 |"mitochondrion inheritance"                          |FBgn0033690, FBgn0261618                                                                               | 0.4523810|
|ID:0000002 |"mitochondrial genome maintenance"                   |FBgn0004407, FBgn0010438, FBgn0032154, FBgn0039930, FBgn0040268, FBgn0013674, FBgn0037008, FBgn0003116 | 0.0256410|
|ID:0000009 |"alpha-1,6-mannosyltransferase activity"             |FBgn0037743, FBgn0035401                                                                               |        NA|
|ID:0000010 |"trans-hexaprenyltranstransferase activity"          |FBgn0037044, FBgn0051005                                                                               | 0.8587571|
|ID:0000012 |"single strand break repair"                         |FBgn0026737, FBgn0026751, FBgn0038704                                                                  | 0.2820513|
|ID:0000014 |"single-stranded DNA endodeoxyribonuclease activity" |FBgn0002887, FBgn0028434, FBgn0030170, FBgn0263831                                                     | 0.1901596|
|ID:0000015 |"phosphopyruvate hydratase complex"                  |FBgn0000579                                                                                            | 0.1410256|
```
