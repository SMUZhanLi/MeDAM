#' @title Search disease-related genes
#' @description Retrieve disease-related genes.
#' @param medam a Pool object connected MeDAM.db
#' @param doid Disease Ontology ID (DOID).
#' @importFrom clusterProfiler bitr
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @details The aliases of disease can be converted to DOID using
#' \code{\link{disease2doid}}.
#' @return entrez gene id
#' @examples
#' \dontrun{
#' # remotes::install_github("rstudio/pool")
#' library(pool)
#' library(dplyr)
#'
#' medamdb <- dbPool(drv = RSQLite::SQLite(), dbname = "MeDAM.db")
#' pe_doid <- disease2doid(medamdb, "pre-eclampsia")
#' pe_doid <- pull(pe_doid, doid) ## DOID:10591
#' pe_drg <- drgene_search(medamdb, pe_doid)
#' poolClose(medamdb)
#' }
#' @export
drgene_search <- function(medam, doid) {
  drg <- medam |>
    dbquery("disease2gene", doid == !!doid) |>
    pull(entrezid) |>
    strsplit(split = ",") |>
    unlist()
  drg <- suppressMessages(
    bitr(drg, "ENTREZID", "SYMBOL", org.Hs.eg.db, drop = FALSE))
  drg <- as_tibble(drg)
  return(drg)
}

#' @title Disease-related genes ORA
#' @description Over-representation analysis (ORA) for disease-related genes.
#' Ref: fig 1.C in Huimin Zheng, et al. (2022).
#' @param nodes The tibble \code{nodes} from STRING or STITCH network.
#' @param drg Disease-related genes (entrez gene id).
#' @param universe Count of universal background protein-coding genes.
#' For example, background of homo sapiens is 21306.
#' @importFrom stats phyper
#' @return A list contained:
#' * \code{ngene} number of interaction genes (proteins) in "nodes"
#' * \code{indrg} number of intersection between genes (proteins) in "nodes" and
#' disease-related genes.
#' * \code{pvalue} p value of hypergeometric test
#' * \code{overlap} overlap genes between genes (proteins) in "nodes" and
#' disease-related genes
#' @references Huimin Zheng, et al. (2022). In silico method to maximise the
#' biological potential of understudied metabolomic biomarkers: a study in
#' pre-eclampsia.
#' @export
drgene_ora <- function(nodes, drg, universe) {
  res <- list(ngene = 0, indrg = 0, pvalue = NA, overlap = NA)
  if (!is.null(nodes)) {
    iprot <- nodes |> filter(type == "protein")
    entrezid <- pull(iprot, external)
    q <- sum(entrezid %in% drg) - 1
    m <- length(drg)
    n <- universe - m
    k <- length(entrezid)
    res$ngene <- k
    if (q > 0) {
      res$indrg <- q + 1
      res$pvalue <- phyper(q, m, n, k, lower.tail = FALSE)
      res$overlap <- iprot |>
        filter(external %in% drg) |>
        pull(name) |>
        paste0(collapse = "/")
    }
  }
  return(res)
}
