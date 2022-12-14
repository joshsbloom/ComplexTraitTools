% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/xQTL_sim.R
\name{getCrossVCF}
\alias{getCrossVCF}
\title{Get a subset of a VCF for two input parents from a biparental cross
Read in a reference VCF file, assuming variant calls are assuming ploidy = 1,
only extract information at biallelic variants, and ensure that only sites that vary
between the two parents are output.}
\usage{
getCrossVCF(ref.vcf, p1.name, p2.name)
}
\arguments{
\item{ref.vcf}{Path to reference vcf file}

\item{p1.name}{string corresponding to parent designated as parent 1 (column name needs to exist in vcf )}

\item{p2.name}{string corresponding to parent designated as parent 2 (column name needs to exist in vcf )}
}
\value{
a vcfR object for the biallelic variants for that cross
}
\description{
Get a subset of a VCF for two input parents from a biparental cross
Read in a reference VCF file, assuming variant calls are assuming ploidy = 1,
only extract information at biallelic variants, and ensure that only sites that vary
between the two parents are output.
}
