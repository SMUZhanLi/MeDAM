#' @title Search disease-related compounds
#' @description Retrieve disease-related compounds.
#' @param medam a Pool object connected MeDAM.db
#' @param cid PubChem cid, e.g. '2244'.
#' @param doid Disease Ontology ID (DOID), e.g. 'DOID:178'.
#' @details Disease-related compounds collected from Human Metabolome Database
#' (HMDB) and Comparative Toxicogenomics Database (CTD).
#' @return A table in which each row represents the compound-disease (cid-doid)
#' association.
#' @examples
#' \dontrun{
#' # remotes::install_github("rstudio/pool")
#' library(pool)
#' library(dplyr)
#'
#' medamdb <- dbPool(drv = RSQLite::SQLite(), dbname = "MeDAM.db")
#' cid <- c("1", "3715")
#' doid <- c("DOID:83", "DOID:3021")
#' cid2doid <- drcid_search(medamdb, cid, doid)
#' poolClose(medamdb)
#' }
#' @export
drcid_search <- function(medam, cid, doid) {
  cid2doid <- medam |>
    dbquery("cid2doid", cid %in% !!cid, doid %in% !!doid)
  return(cid2doid)
}
