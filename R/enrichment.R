#' @title Enrichment analysis
#' @description Functional or pathway enrichment analysis for genes/proteins
#' interactied with metabolites. Ref: fig 1.C in Huimin Zheng, et al. (2022).
#' @param nodes The tibble \code{nodes} from STRING or STITCH network.
#' @param padj P-values adjusted methods, one of "holm", "hochberg", "hommel",
#' "bonferroni", "BH", "BY", "fdr", "none".
#' @param pcutoff 0.05.
#' @param ont c("BP", "MF", "CC", "KEGG", "WikiPathways", "Reactome", "DO").
#' @importFrom clusterProfiler enrichGO enricher
#' @importFrom ReactomePA enrichPathway
#' @importFrom DOSE enrichDO
#' @return A \code{\link[DOSE]{enrichResult-class}} instance for GO-BP, GO-MF,
#' GO-CC, KEGG, WikiPathways, Reactome and DO (Disease Ontology) enrichment
#' analysis.
#' @references Huimin Zheng, et al. (2022). In silico method to maximise the
#' biological potential of understudied metabolomic biomarkers: a study in
#' pre-eclampsia.
#' @export
enrichment <- function(nodes,
                       padj = "BH",
                       pcutoff = 0.05,
                       ont = c("BP", "MF", "CC", "KEGG",
                               "WikiPathways", "Reactome", "DO")) {
  res <- vector("list", length(ont))
  names(res) <- ont
  if (!is.null(nodes)) {
    iprot <- nodes |> filter(type == "protein")
    entrezid <- pull(iprot, external)
    for (x in c("BP", "MF", "CC")) {
      if (x %in% ont) {
        res[[x]] <- enrichGO(entrezid,
                        OrgDb = "org.Hs.eg.db",
                        ont = x,
                        pAdjustMethod = padj,
                        pvalueCutoff = pcutoff,
                        qvalueCutoff = 1)
        # res[[x]] <- simplify(ego)
      }
    }
    if ("KEGG" %in% ont) {
      res$KEGG <- enricher(entrezid,
                           pAdjustMethod = padj,
                           pvalueCutoff = pcutoff,
                           qvalueCutoff = 1,
                           gson = kegg_gson_hsa)
    }
    if ("WikiPathways" %in% ont) {
      res$WikiPathways <- enricher(entrezid,
                                   pAdjustMethod = padj,
                                   pvalueCutoff = pcutoff,
                                   qvalueCutoff = 1,
                                   gson = wp_gson_homo)
    }
    if ("Reactome" %in% ont) {
      res$Reactome <- enrichPathway(entrezid,
                                    pAdjustMethod = padj,
                                    pvalueCutoff = pcutoff,
                                    qvalueCutoff = 1)
    }
    if ("DO" %in% ont) {
      res$DO <- enrichDO(entrezid,
                         pAdjustMethod = padj,
                         pvalueCutoff = pcutoff,
                         qvalueCutoff = 1)
    }
    res <- lapply(res, entrez2gename, iprot = iprot)
  }
  return(res)
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