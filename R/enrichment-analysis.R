#' @title Enrichment analysis
#' @description Functional or pathway enrichment analysis for genes/proteins
#' interactied with metabolites. Ref: fig 1.C in Huimin Zheng, et al. (2022).
#' @param nodes The tibble \code{nodes} from STRING or STITCH network.
#' @param padj P-values adjusted methods, one of "holm", "hochberg", "hommel",
#' "bonferroni", "BH", "BY", "fdr", "none".
#' @param pcutoff 0.05.
#' @param ont c("BP", "MF", "CC", "KEGG", "WikiPathways", "Reactome", "DO").
#' @importFrom clusterProfiler enricher
#' @return A \code{\link[DOSE]{enrichResult-class}} instance for GO-BP, GO-MF,
#' GO-CC, KEGG, WikiPathways, Reactome and DO (Disease Ontology) enrichment
#' analysis.
#' @references Huimin Zheng, et al. (2022). In silico method to maximise the
#' biological potential of understudied metabolomic biomarkers: a study in
#' pre-eclampsia.
#' @export
enrichment_analysis <- function(nodes,
                                padj = "BH",
                                pcutoff = 0.05,
                                ont = c("BP", "MF", "CC", "KEGG",
                                        "WikiPathways", "Reactome", "DO")) {
  if (!is.null(nodes)) {
    iprot <- nodes |> filter(type == "protein")
    entrezid <- pull(iprot, external)
    res <- lapply(ont, function(term) {
      gsid2gene <- medam_terms[[term]]$gsid2gene
      gsid2name <- medam_terms[[term]]$gsid2name
      enricher(entrezid,
               pAdjustMethod = padj,
               pvalueCutoff = pcutoff,
               qvalueCutoff = 1,
               TERM2GENE = gsid2gene,
               TERM2NAME = gsid2name)
    })
    names(res) <- ont
    res <- lapply(res, entrez2gename, iprot = iprot)
  }
  return(res)
}

#' @title Add result of enrichment analysis
#' @description Add result of enrichment analysis into the slot `$enriment`
#' in result of MeDAM batch search when select a metabolite and its branch.
#' @param medamres Result of MeDAM batch search.
#' @param metabolite Metabolite name.
#' @param branch MeDAM branch.
#' @param padj P-values adjusted methods, one of "holm", "hochberg", "hommel",
#' "bonferroni", "BH", "BY", "fdr", "none".
#' @param pcutoff 0.05.
#' @param ont c("BP", "MF", "CC", "KEGG", "WikiPathways", "Reactome", "DO").
#' @importFrom clusterProfiler enricher
#' @return Result of MeDAM batch search.
#' @export
medam_add_enrichment <- function(medamres,
                                 metabolite,
                                 branch,
                                 padj = "BH",
                                 pcutoff = 0.05,
                                 ont = c("BP", "MF", "CC", "KEGG",
                                         "WikiPathways", "Reactome", "DO")) {
  nodes <- medamres[[metabolite]][[branch]]$nodes
  enres <- enrichment_analysis(nodes, padj, pcutoff, ont)
  medamres$enriment[[metabolite]][[branch]] <- enres
  return(medamres)
}

## Convert entrez gene id to gene symbol for \code{enrichResult}.
#' @keywords internal
entrez2gename <- function(enrichres, iprot) {
  if (inherits(enrichres, "enrichResult")) {
    res <- enrichres@result
    entrezid <- strsplit(res$geneID, split = "/")
    gename <- lapply(entrezid, function(x) {
      i <- match(x, iprot$external)
      paste0(iprot$name[i], collapse = "/")
    })
    res$geneID <- unlist(gename)
    enrichres@result <- res
  }
  return(enrichres)
}