% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Plotting.R, R/Utils.R
\name{InputOutputFunctions}
\alias{InputOutputFunctions}
\alias{reshapeResults}
\alias{plotGraph}
\alias{plotBarplot}
\alias{plotHeatmap}
\alias{readGmtFileAsDataFrame}
\alias{saveDataFrameAsGmtFile}
\alias{readEdbFileAsDataFrame}
\alias{filterOntology}
\alias{generateInputData}
\title{Input/Output Functions}
\usage{
reshapeResults(
  mulea_model = NULL,
  mulea_model_resuts = NULL,
  mulea_model_ontology_col_name = "ontologyId",
  mulea_model_resuts_ontology_col_name = "ontologyId",
  category_stat_column_name = "adjustedPValue",
  cut_off_to_test_data = TRUE
)

plotGraph(
  mulea_relaxed_resuts,
  edge_weight_cutoff = 0,
  statistics_value_colname = "adjustedPValue",
  ontology_id_column_name = "ontologyId",
  gen_id_in_ontology_column_name = "genIdInOntology",
  statistics_value_cutoff = 0.05
)

plotBarplot(
  mulea_relaxed_resuts,
  selection_vector = NULL,
  categories_names = "ontologyId",
  statistics_value_colname = "adjustedPValue",
  statistics_value_cutoff = 0.05
)

plotHeatmap(
  mulea_relaxed_resuts,
  statistics_value_colname = "adjustedPValue",
  gen_id_in_ontology_column_name = "genIdInOntology",
  statistics_value_cutoff = 0.05
)

readGmtFileAsDataFrame(gmtFilePath)

saveDataFrameAsGmtFile(modelDF, gmtFilePath)

readEdbFileAsDataFrame(edbFilePath)

filterOntology(input_gmt, min = NULL, max = NULL)

generateInputData(
  input_gmt,
  sample_ratio = 0.5,
  group_under_over_representation_ratio = 0.9,
  number_of_over_representation_groups = 1,
  number_of_under_representation_groups = 1,
  turn_on_log = FALSE
)
}
\arguments{
\item{mulea_model}{MulEA object represents model. For example created by MulEA::ORA.}

\item{mulea_model_resuts}{Results from model, in most cases it is returned by MulEA::runTest generic method.}

\item{mulea_relaxed_resuts}{data.table in relaxed form.}

\item{selection_vector}{vector for selecting variables to plot.}

\item{gmtFilePath}{path with name of file, where the file is localized or where to save model. Example: "/R/MulEA/extdata/model.gmt"}

\item{modelDF}{data frame with model.}

\item{edbFilePath}{path with name of file, where the file is localized or where to save model.}

\item{input_gmt}{input dataframe, read from gmt file.}

\item{min}{minimum size of term. Default 20% from quantile on term size distribution.}

\item{max}{maximum size of term. Default 80% from quantile on term size distribution.}

\item{sample_ratio}{}

\item{number_of_under_representation_groups}{}
}
\value{
Return detailed and relaxed datatable where model and results are merged for plotting purposes.

Return plot.

Return plot.

Return plot.

Return data frame with model from specific location.

Return gmt file under specific location which include model in gmt format.

Return data frame with model from specific location.

Return data frame with model from specific location.

Return data frame with model from specific location.
}
\description{
\code{reshapeResults}

\code{reshapeResults} merge model and model relsuts into 
one relaxed datatable for easy resutls graphical interpretation.

\code{plotGraph}

\code{plotGraph} merge model and model relsuts into 
one relaxed datatable for easy resutls graphical interpretation.

\code{plotBarplot}

\code{plotBarplot} barplot of p-values.

\code{plotHeatmap}

\code{plotHeatmap} merge model and model relsuts into 
one relaxed datatable for easy resutls graphical interpretation.

\code{readGmtFileAsDataFrame}

\code{readGmtFileAsDataFrame} read model in data frame form from gmt file.

\code{saveDataFrameAsGmtFile}

\code{saveDataFrameAsGmtFile} saves copy of the model from data frame in gmt file.

\code{readEdbFileAsDataFrame}

\code{readEdbFileAsDataFrame} read GSEA results in data frame form from .edb file.

\code{filterOntology}

\code{filterOntology} cut ontology to specific terms sizes.

\code{generateInputData}

\code{generateInputData} generate artificial GO with specific terms under or over represented.
}
\examples{
modelDfFromFile <- MulEA::readGmtFileAsDataFrame(gmtFilePath = system.file(package="MulEA", "extdata", "model.gmt"))
MulEA::saveDataFrameAsGmtFile(modelDF = modelDfFromFile, gmtFilePath = paste(system.file(package="MulEA", "extdata"), "fromDb.gmt", sep = "/"))
}
