
# MulEA
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
modelDfFromFile <- MulEA::readGmtFileAsDF(gmtFilePath = paste(muleaPkgDir,"/example/model.gmt", sep = ""))
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
