#' @title Search disease-related compounds
#' @description Retrieve disease-related compounds.
#' @param medam a Pool object connected MeDAM.db
#' @param cid PubChem cid, e.g. '2244'.
#' @param doid Disease Ontology ID (DOID), e.g. 'DOID:178'.
#' @details Disease-related compounds collected from Human Metabolome Database
#' (HMDB) and Comparative Toxicogenomics Database (CTD).
drcid_search <- function(medam, cid, doid) {
  cid2doid <- medam |>
    dbquery("cid2doid", cid %in% !!cid, doid %in% !!doid)
  return(cid2doid)
}