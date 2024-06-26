% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot-circle.R
\name{medam_circle_plot}
\alias{medam_circle_plot}
\title{Circle plot}
\usage{
medam_circle_plot(
  medamres,
  sort.by = "vip",
  max_metabo_len = 30,
  start = 0,
  .offset = 0.1
)
}
\arguments{
\item{medamres}{Result of MeDAM batch search.}

\item{sort.by}{Order the metabolites by selected variable, one of "vip",
"abs(log2fc)", "-log10(pvalue)" and "auc".}

\item{max_metabo_len}{A numeric value sets wrap length for metabolite labels.
It will wrap names longer that 30 characters by default.}

\item{start}{Offset of starting point from 12 o'clock in radians. Offset is
applied clockwise.}

\item{.offset}{Offset of each statistical measure in different circular ring.}
}
\description{
Show statistical measure of metabolites which have
differentially abundance in the two groups.
}
