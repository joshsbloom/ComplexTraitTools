% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/xQTL_stats.R
\name{plotSummary}
\alias{plotSummary}
\title{Plot summary contrast and standard error, as well as -log10(p)}
\usage{
plotSummary(results, simulatedQTL = NULL, effective.n.tests = 600)
}
\arguments{
\item{results}{data.frame of experiment results}

\item{simulatedQTL}{data frame of chromosome, and position of simulated QTL (optional)}

\item{effective.n.tests}{effective number of tests genomewide (integer)}
}
\value{
ggplot object of experiment results
}
\description{
Plot summary contrast and standard error, as well as -log10(p)
}
