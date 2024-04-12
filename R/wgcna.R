#' @title WGCNA modules
#' @description Identify abundance-correlated metabolites (i.e., consistency
#' abundance module) using weighted correlation network analysis (WGCNA).
#' @param abundance Abundance table, a matrix or data frame in which columns are
#' metabolites and rows ar samples.
#' @param networkType Network type, one of "signed", "unsigned" and
#' "signed hybrid".
#' @param corMethod Correlation algorithm for network construction, e.g.,
#' "spearman" or "pearson"
#' @param minModuleSize Minimum size of module or cluster, default is 5.
#' @param cutHeight Maximum dissimilarity (i.e., 1-correlation) that qualifies
#' modules for merging.
#' @param nthreads Default 10.
#' @importFrom WGCNA pickSoftThreshold adjacency mergeCloseModules
#' @importFrom dynamicTreeCut cutreeDynamic
#' @importFrom flashClust flashClust
#' @importFrom stats as.dist
#' @return A \code{tibble} combined metabolites and modules.
#' @examples
#' \dontrun{
#' data(PE.metabolomics)
#' abun <- PE.metabolomics$abundance
#' wgcna_mods <- wgcna_analysis(abun)
#' }
#' @export
wgcna_analysis <- function(abundance,
                           networkType = "signed",
                           corMethod = "spearman",
                           minModuleSize = 5,
                           cutHeight = 0.25,
                           nthreads = 10) {
  powers <- c(1:10, seq(12, 30, 2))
  sft <- pickSoftThreshold(abundance,
                           powerVector = powers,
                           networkType = networkType)
  softPower <- sft$powerEstimate
  adjMat <- adjacency(datExpr = abundance,
                      power = softPower,
                      type = networkType,
                      corFnc = "cor",
                      corOptions = list(method = corMethod))
  TOM <- TOMsimilarity_parallel(adjMat, nthreads)
  dissTOM <- 1 - TOM
  metaboTree <- flashClust(as.dist(dissTOM), method = "complete")
  moduleLabels <- cutreeDynamic(dendro = metaboTree,
                                distM = dissTOM,
                                method = "hybrid",
                                deepSplit = 2,
                                pamRespectsDendro = FALSE,
                                minClusterSize = minModuleSize)
  merge <- mergeCloseModules(
    exprData = abundance,
    colors =  moduleLabels,
    corFnc = "cor",
    corOptions = list(method = corMethod),
    cutHeight = cutHeight
  )
  res <- tibble(metabolite = colnames(abundance),
                module = as.character(merge$colors))
  return(res)
}

## get co-abundance compound
#' @keywords internal
get_cocompound <- function(wgcna, compound) {
  mod <- wgcna |>
    filter(metabolite == compound) |>
    pull(module)
  cocompound <- wgcna |>
    filter(module == mod) |>
    pull(metabolite)
}