#' @title Query database method
#' @description Query a database using dplyr and pool.
#' @param dbconn A Pool object connecting database
#' @param dbtbl A table in the database
#' @param ... The arguments passed on to \code{\link[dplyr]{filter}} and allows
#' you to do queries using dplyr that can be transformed into a SELECT SQL
#' statement.
#' @importFrom dplyr tbl filter collect
#' @return tibble object
#' @examples
#' \dontrun{
#' library(dplyr)
#'
#' medamdb <- medamdbPool("MeDAM.db")
#' disease <- c("Pre-eclampsia", "Endometriosis")
#' medamdb |>
#'   dbquery("disease2doid", disease %in% !!disease)
#' poolClose(medamdb)
#' }
#' @export
dbquery <- function(dbconn, dbtbl, ...) {
  dbconn |>
    tbl(dbtbl) |>
    filter(...) |>
    collect()
}

#' @title MeDAM dbPool connections
#' @description Create a pool of database connections for MeDAM.db
#' @param dbname MeDAM.db
#' @importFrom RSQLite SQLite
#' @importFrom pool dbPool
#' @return SQLiteConnection objects
#' @examples
#' \dontrun{
#' medamdb <- medamdbPool(dbname = "MeDAM.db")
#' }
#' @export

medamdbPool <- function(dbname) {
  dbPool(drv = RSQLite::SQLite(), dbname = dbname)
}