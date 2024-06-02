#' @title Structurally similar compounds
#' @description Search compounds shared similarity structure.
#' @param medam a Pool object connected MeDAM.db
#' @param cid PubChem cid, e.g., 2244.
#' @details  \code{ssim_search} is based on substructure key-based 2D Tanimoto
#' similarity search (tanimoto >= 90%). The synonyms of compound can be
#' converted to "cid" (PubChem cid) using \code{\link{compound2cid}}.
#' @return A tibble contained query cid and ssimcid (compounds shared structural
#' similarity in 2D).
#' @examples
#' \dontrun{
#' # remotes::install_github("rstudio/pool")
#' library(pool)
#' library(dplyr)
#'
#' medamdb <- dbPool(drv = RSQLite::SQLite(), dbname = "MeDAM.db")
#' c2cid <- compound2cid(medamdb, "ONONETIN")
#' cid <- pull(c2cid, cid)
#' cid2ssim <- ssim_search(medamdb, cid)
#' poolClose(medamdb)
#' }
#' @export
ssim_search <- function(medam, cid) {
  c2ssimcid <- medam |>
    dbquery("cid2ssimcid", cid %in% !!cid) |>
    separate_longer_delim(ssimcid, delim = ",") |>
    mutate(stitch = cid2stitchid(ssimcid)) |>
    filter(!is.na(stitch))
  ssimm_stitchid <- c2ssimcid |> pull(stitch)
  simm2cheminfo <- medam |>
    dbquery("stitchv5compoundinfo", stitchv5 %in% !!ssimm_stitchid) |>
    rename(stitch = 1) |>
    select(stitch, name, description)
  c2ssimcid <- c2ssimcid |>
    left_join(simm2cheminfo, by = "stitch") |>
    select(cid, ssimcid, stitch, name, description)
  return(c2ssimcid)
}