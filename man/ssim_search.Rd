% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ssim-search.R
\name{ssim_search}
\alias{ssim_search}
\title{Structurally similar compounds}
\usage{
ssim_search(medam, cid)
}
\arguments{
\item{medam}{a Pool object connected MeDAM.db}

\item{cid}{PubChem cid, e.g., 2244.}
}
\value{
A tibble contained query cid and ssimcid (compounds shared structural
similarity in 2D).
}
\description{
Search compounds shared similarity structure.
}
\details{
\code{ssim_search} is based on substructure key-based 2D Tanimoto
similarity search (tanimoto >= 90\%). The synonyms of compound can be
converted to "cid" (PubChem cid) using \code{\link{compound2cid}}.
}
\examples{
\dontrun{
# remotes::install_github("rstudio/pool")
library(pool)
library(dplyr)

medamdb <- dbPool(drv = RSQLite::SQLite(), dbname = "MeDAM.db")
c2cid <- compound2cid(medamdb, "ONONETIN")
cid <- pull(c2cid, cid)
cid2ssim <- ssim_search(medamdb, cid)
poolClose(medamdb)
}
}
