#' @title Scaling abundance table
#' @description Wrapper for "X-Scaling" in R pkg "ropls".
#' @param xMN Abundance table, a matrix or data frame in which columns are
#' metabolites and rows ar samples.
#' @param scaleC Scaling methods for "abundance", one of "none", "center",
#' "pareto", and "standard" (default).
#' @importFrom stats sd
#' @export
scale_abundance <- function(xMN, scaleC) {
  xMeanVn <- apply(xMN, 2, function(colVn) mean(colVn, na.rm = TRUE))
  switch(scaleC,
         none = {
           xMeanVn <- rep(0, ncol(xMN))
           xSdVn <- rep(1, times = ncol(xMN))
         },
         center = {
           xSdVn <- rep(1, times = ncol(xMN))
         },
         pareto = {
           xSdVn <- apply(xMN, 2, function(colVn) sqrt(sd(colVn, na.rm = TRUE)))
         },
         standard = {
           xSdVn <- apply(xMN, 2, function(colVn) sd(colVn, na.rm = TRUE))
         })
  xMN <- scale(xMN, center = xMeanVn, scale = xSdVn)
  return(xMN)
}