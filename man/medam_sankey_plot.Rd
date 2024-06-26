% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot-sankey.R
\name{medam_sankey_plot}
\alias{medam_sankey_plot}
\title{Sankey plot}
\usage{
medam_sankey_plot(
  medamres,
  metabolite,
  branch,
  disterms,
  matchkey,
  max_metabo_len = 30,
  max_desc_len = 30
)
}
\arguments{
\item{medamres}{Result of MeDAM batch search.}

\item{metabolite}{Metabolite.}

\item{branch}{One of "target_proteins", "ssim_metabolites" and
"coab_metabolites"}

\item{disterms}{The disease related functional terms or pathways.}

\item{matchkey}{"Description" or "ID".}

\item{max_metabo_len}{A numeric value sets wrap length for labels of
metabolites. It will wrap names longer that 30 characters by default.}

\item{max_desc_len}{A numeric value sets wrap length for labels of terms.
It will wrap names longer that 30 characters by default.}
}
\description{
Show interactions between: (1) the potential proteins and the
query nodes in MeDAM branch (e.g. the target proteins, the structural similar
metabolites and the co-abundant metabolites); (2) the potential proteins and
functional terms.
}
