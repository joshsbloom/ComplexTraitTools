% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/xQTL_sim.R
\name{generateCounts}
\alias{generateCounts}
\title{Simulate reference-based 0/1 counts for ref and alt alleles}
\usage{
generateCounts(y, sx, sel.frac, depth, lower.tail = F)
}
\arguments{
\item{y}{vector of simulated phenotype, scaled so total phenotypic variance is 1}

\item{sx}{XSnpMatrix}

\item{sel.frac}{fraction of population selected}

\item{depth}{average depth per variant site (summed over both alleles)}

\item{lower.tail}{selecting lower tail (T/F, default=F)}
}
\value{
data.frame of variant ID, expected allele frequence, ref counts and alt counts
}
\description{
Simulate reference-based 0/1 counts for ref and alt alleles
}
