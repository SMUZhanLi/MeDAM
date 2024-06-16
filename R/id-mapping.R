#' @title compound2cid
#' @description Convert synonyms of compound to PubChem cid and STITCH cid
#' (which are are derived from PubChem).
#' @param medam A Pool object connected MeDAM.db
#' @param compound Synonyms of compound (ignored case), e.g. "aspirin" or
#' some external id (see details).
#' @importFrom dplyr pull mutate group_by arrange slice_head ungroup left_join select
#' @importFrom tidyr separate_longer_delim
#' @details some external id of "aspirin":
#' * \code{KEGG compound} D00109
#' * \code{ChEBI} CHEBI:15365
#' * \code{HMDB} HMDB0001879
#' * \code{CAS-RN} 50-78-2
#' * etc.
#' @examples
#' \dontrun{
#' # remotes::install_github("rstudio/pool")
#' library(pool)
#'
#' medamdb <- dbPool(drv = RSQLite::SQLite(), dbname = "MeDAM.db")
#' metabolites <- c("aspirin", "D00109", "CHEBI:15365", "HMDB0001879")
#' c2cid <- compound2cid(medamdb, metabolites)
#' poolClose(medamdb)
#' }
#' @return a tibble contained query compound, PubChem cid and STITCH cid.
#' @export
compound2cid <-  function(medam, compound) {
  is_cid <- grepl("^CID\\d+$", compound)
  cid <- sub("^CID", "", compound[is_cid])
  c2cid <- tibble(compound = compound[is_cid], cid = cid,
                  stitch = cid2stitchid(cid))
  non_cid <- compound[!is_cid]
  if (length(non_cid) > 0) {
    c2md5sum <- compound2md5sum(non_cid)
    md5sum <- pull(c2md5sum, md5sum)
    md5sum2cid <- medam |>
      dbquery("compound2cid", md5sum %in% !!md5sum) |>
      separate_longer_delim(cid, delim = ",") |>
      mutate(stitch = cid2stitchid(cid)) |>
      group_by(md5sum) |>
      arrange(as.integer(cid), .by_group = TRUE) |>
      slice_head(n = 1) |>
      ungroup()
    c2cid <- c2md5sum |>
      left_join(md5sum2cid, by = "md5sum") |>
      select(compound, cid, stitch) |>
      bind_rows(c2cid)
  }
  return(c2cid)
}

## Convert synonyms of compound to md5sum string.
#' @importFrom digest digest
#' @importFrom tibble tibble
#' @keywords internal
compound2md5sum <- function(compound) {
  md5sum <- compound |>
    tolower() |>
    lapply(digest, algo = "md5", serialize = FALSE) |>
    unlist()
  compound2md5sum <- tibble(compound = compound, md5sum = md5sum)
  return(compound2md5sum)
}

## Convert PubChem cid to STTICH id.
#' @keywords internal
cid2stitchid <- function(cid) {
  cnum <- as.integer(cid)
  cnum[is.na(cnum) | cnum >= 10^8] <- 0
  stitchid <- sprintf("CIDs%08d", cnum)
  i <- which(cnum == 0)
  stitchid[i] <- NA
  return(stitchid)
}

#' @title protein2stringid
#' @description Convert synonyms of protein to STRING id.
#' @param medam A Pool object connected MeDAM.db.
#' @param protein Synonyms of protein(ignored case), e.g. "TP53".
#' @importFrom dplyr rename
#' @return A tibble contained query protein and STRING id.
#' @examples
#' \dontrun{
#' # remotes::install_github("rstudio/pool")
#' library(pool)
#'
#' medamdb <- dbPool(drv = RSQLite::SQLite(), dbname = "MeDAM.db")
#' proteins <- c("PTCH1", "TP53", "BRCA1", "BRCA2")
#' protein2stringid(medamdb, proteins)
#' poolClose(medamdb)
#' }
#' @export
protein2stringid <- function(medam, protein) {
  prot2stringid <- medam |>
    dbquery("stringv12protein2stringid",
            protein %in% !!protein) |>
    rename(string = 2) |>
    mutate(tmpgroup = toupper(protein)) |>
    group_by(tmpgroup) |>
    arrange(string, .by_group = TRUE) |>
    slice_head(n = 1) |>
    ungroup() |>
    select(-tmpgroup)
  return(prot2stringid)
}

#' @title disease2doid
#' @description Convert aliases of disease to Diease Ontology id (DOID).
#' @param medam A Pool object connected MeDAM.db.
#' @param disease Aliases of disease(ignored case), e.g. "pre-eclampsia".
#' @param fixed	If TRUE match "disease" exactly, otherwise use \code{str_like()}
#' to query disease fuzzy faintly. See details.
#' @importFrom stringr str_like
#' @details Fuzzy matching follows the conventions of the SQL LIKE operator:
#' * Must match the entire string.
#' * \code{⁠_}⁠ matches a single character (like \code{.}).
#' * \code{⁠%}⁠ matches any number of characters (like \code{⁠.*}⁠).
#' * \code{⁠\%}⁠ and \code{⁠\_}⁠ match literal \code{⁠%}⁠ and \code{⁠_}⁠.
#' * The match is case insensitive by default.
#' @seealso \code{\link[stringr]{str_like}}
#' @return A tibble contained DOID, name and superclass of query disease.
#' @examples
#' \dontrun{
#' # remotes::install_github("rstudio/pool")
#' library(pool)
#'
#' medamdb <- dbPool(drv = RSQLite::SQLite(), dbname = "MeDAM.db")
#' disease2doid(medamdb, "pre-eclampsia")
#' ## fuzzy matching
#' disease2doid(medamdb, "%eclampsia", fixed = FALSE)
#' poolClose(medamdb)
#' }
#' @export
disease2doid <-  function(medam, disease, fixed = TRUE) {
  if (grepl("DOID:\\d+", disease)) {
    d2doid <- medam |>
      dbquery("disease2doid", doid == !!disease)
  } else {
    if (fixed) {
      d2doid <- medam |>
        dbquery("disease2doid", disease == !!disease)
    } else {
      disease <- paste0("%", disease, "%")
      d2doid <- medam |>
        dbquery("disease2doid", str_like(disease, !!disease))
    }
  }
  return(d2doid)
}