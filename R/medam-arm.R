#' @title The target proteins branch in MeDAM pipeline
#' @description The first arm of pipeline in medam. Ref: fig 1.C in
#' Huimin Zheng, et al. (2022).
#' @param medam A Pool object connected MeDAM.db.
#' @param compound Target compound, e.g. Ononetin.
#' @param protein Target proteins regulated by compound and its intermediates,
#' e.g. TRPM3 (be inhibited by Ononetin) and TRPA1(be activated).
#' @param score Threshold of significance to include a interaction, a number
#' between 0 and 1000 (default: 700).
#' @param use_ecpi \code{TRUE}, see details.
#' @param ecpi_score Threshold of significance for metabilite-protein
#' experimental interactions (from 0 to 1000, default: 800).
#' @param drg Disease-related genes(entrez gene id).
#' @param universe Count of universal background protein-coding genes.
#' For example, background of homo sapiens is 21306.
#' @param padj P-values adjusted methods, one of "holm", "hochberg", "hommel",
#' "bonferroni", "BH", "BY", "fdr", "none".
#' @param pcutoff 0.05.
#' @param ont c("BP", "MF", "CC", "KEGG", "WikiPathways", "Reactome", "DO").
#' @details In general, the \code{protein} obtained from literature and database
#'  meant a increase in work. So, we recommend to search protein regulated by
#' target compound using \code{stitch_expt_cpi()} (STITCH metabilite-protein
#' experimental interactions network) required a \code{ecpi_score} (e.g. 700).
#' @seealso \code{\link{stitch_expt_cpi}}, \code{\link{string_network}}
#' @return list contained weights of each edge, annotations of each node
#' result of disease-related genes ORA analysis and result of enrichment
#' analysis as following:
#'
#' the tibble \code{edges} with fields:
#' * \code{node1, node2} protein- or compound-nodes in the network
#' * \code{type} type of interaction edges in network, e.g., cpi
#' (compound-protein), cci (compound-compound) and ppi (protein-protein)
#' * \code{score} combined score
#' * \code{experiment} experimental score
#' * \code{neighborhood} gene neighborhood score
#' * \code{fusion} gene fusion score
#' * \code{prediction} prediction score
#' * \code{cooccurence} cooccurence score
#' * \code{coexpression} coexpression score
#' * \code{database} database score
#' * \code{textmining} textmining score
#' * \code{evidence} more convincing interaction, the rank are: 1) experiment,
#' 2) neighborhood/fusion, 3) prediction/cooccurence/coexpression,
#' 4) database/textmining
#'
#' the tibble \code{nodes} with fields:
#' * \code{id} STRING id (compound) or STTICH id (protein)
#' * \code{type} type of nodes in network, protein or metabolite
#' * \code{name} common synonyms of compound or protein
#' * \code{external} external id of compound (PubChem cid) or protein
#' (entrez gene id)
#' * \code{description} description of compound or protein
#' @references Huimin Zheng, et al. (2022). In silico method to maximise the
#' biological potential of understudied metabolomic biomarkers: a study in
#' pre-eclampsia.
#' @export
medam_target_proteins <- function(medam,
                                  compound,
                                  protein = NULL,
                                  score = 700,
                                  use_ecpi = TRUE,
                                  ecpi_score = 800,
                                  drg = NULL,
                                  universe = 21306,
                                  padj = "BH",
                                  pcutoff = 0.05,
                                  ont = c("BP", "MF", "CC", "KEGG",
                                          "WikiPathways", "Reactome", "DO")) {
  stringid <- c()
  if (use_ecpi) {
    stitchid <- compound2cid(medam, compound) |> pull(stitch)
    ecpi <- stitch_expt_cpi(medam, stitchid, ecpi_score)
    ecpi_prot <- pull(ecpi, stringv12)
    stringid <- c(stringid, ecpi_prot)
  }
  if (!is.null(protein)) {
    prior_prot <- protein2stringid(medam, protein) |>
      pull(string) |>
      unique()
    stringid <- c(stringid, prior_prot)
  }
  if (length(stringid) > 0) {
    res <- string_network(medam, stringid, score)
    res$drg <- drgene_ora(res$nodes, drg, universe)
    res$enrichment <- enrichment_analysis(res$nodes, padj, pcutoff, ont)
  } else {
    res <- NULL
  }
  return(res)
}

#' @title The co-abundance metabolites branch in MeDAM pipeline
#' @description the branch 1 of second arm of pipeline in medam. Ref: fig 1.C in
#' Huimin Zheng, et al. (2022).
#' @param medam A Pool object connected MeDAM.db
#' @param compound A single target compound or a cluster of compounds correlated
#' with abundance (i.e., consistency module).
#' @param score Threshold of significance to include a interaction, a number
#' between 0 and 1000 (default: 700).
#' @param use_wgcna \code{TRUE}, see details.
#' @param abundance Abundance table, a matrix or data frame in which columns are
#' metabolites and rows ar samples.
#' @param networkType Network type, one of "signed", "unsigned" and
#' "signed hybrid".
#' @param corMethod Correlation algorithm for network construction, e.g.,
#' "spearman" or "pearson".
#' @param minModuleSize Maximum size of module or cluster, default is 5.
#' @param cutHeight Maximum dissimilarity (i.e., 1-correlation) that qualifies
#' modules for merging.
#' @param drg Disease-related genes(entrez gene id).
#' @param universe Count of universal background protein-coding genes.
#' For example, background of homo sapiens is 21306.
#' @param padj P-values adjusted methods, one of "holm", "hochberg", "hommel",
#' "bonferroni", "BH", "BY", "fdr", "none".
#' @param pcutoff 0.05.
#' @param ont c("BP", "MF", "CC", "KEGG", "WikiPathways", "Reactome", "DO").
#' @details Use \code{wgcna_analysis()} to calculate consistency modules 
#' (co-abundance compounds) and select the module contained target compound
#' while length of \code{compound} is 1, \code{use_wgcna} is TRUE and
#' \code{abundance} is not null.
#' @seealso \code{\link{wgcna_analysis}}, \code{\link{stitch_network}}
#' @return list contained weights of each edge, annotations of each node
#' result of disease-related genes ORA analysis and result of enrichment
#' analysis as following:
#'
#' the tibble \code{edges} with fields:
#' * \code{node1, node2} protein- or compound-nodes in the network
#' * \code{type} type of interaction edges in network, e.g., cpi
#' (compound-protein), cci (compound-compound) and ppi (protein-protein)
#' * \code{score} combined score
#' * \code{experiment} experimental score
#' * \code{neighborhood} gene neighborhood score
#' * \code{fusion} gene fusion score
#' * \code{prediction} prediction score
#' * \code{cooccurence} cooccurence score
#' * \code{coexpression} coexpression score
#' * \code{database} database score
#' * \code{textmining} textmining score
#' * \code{evidence} more convincing interaction, the rank are: 1) experiment,
#' 2) neighborhood/fusion, 3) prediction/cooccurence/coexpression,
#' 4) database/textmining
#'
#' the tibble \code{nodes} with fields:
#' * \code{id} STRING id (compound) or STTICH id (protein)
#' * \code{type} type of nodes in network, protein or metabolite
#' * \code{name} common synonyms of compound or protein
#' * \code{external} external id of compound (PubChem cid) or protein
#' (entrez gene id)
#' * \code{description} description of compound or protein
#' @references Huimin Zheng, et al. (2022). In silico method to maximise the
#' biological potential of understudied metabolomic biomarkers: a study in
#' pre-eclampsia.
#' @export
medam_coab_metabolites <- function(medam,
                                   compound,
                                   score = 700,
                                   use_wgcna = TRUE,
                                   abundance = NULL,
                                   networkType = "signed",
                                   corMethod = "spearman",
                                   minModuleSize = 5,
                                   cutHeight = 0.25,
                                   drg = NULL,
                                   universe = 21306,
                                   padj = "BH",
                                   pcutoff = 0.05,
                                   ont = c("BP", "MF", "CC", "KEGG",
                                           "WikiPathways", "Reactome", "DO")) {
  if (length(compound) == 1 && use_wgcna && !is.null(abundance)) {
    wgcna <- wgcna_analysis(abundance, networkType, corMethod,
                            minModuleSize, cutHeight)
    cocompound <- get_cocompound(wgcna, compound)
  } else {
    cocompound <- compound
  }
  stitchid <- compound2cid(medam, cocompound) |> pull(stitch)
  if (length(stitchid) > 0) {
    res <- stitch_network(medam, stitchid, score = score)
    res$drg <- drgene_ora(res$nodes, drg, universe)
    res$enrichment <- enrichment_analysis(res$nodes, padj, pcutoff, ont)
  } else {
    res <- NULL
  }
  return(res)
}

#' @title The structural similar metabolites branch in MeDAM pipeline
#' @description The branch 2 of second arm of pipeline in medam. Ref: fig 1.C in
#' Huimin Zheng, et al. (2022).
#' @param medam A Pool object connected MeDAM.db
#' @param compound A single target compound or a set of compounds shared similar
#' structure.
#' @param score Threshold of significance to include a interaction, a number
#' between 0 and 1000 (default: 700).
#' @param search_ssim \code{TRUE}, see details.
#' @param drg Disease-related genes(entrez gene id).
#' @param universe Count of universal background protein-coding genes.
#' For example, background of homo sapiens is 21306.
#' @param padj P-values adjusted methods, one of "holm", "hochberg", "hommel",
#' "bonferroni", "BH", "BY", "fdr", "none".
#' @param pcutoff 0.05.
#' @param ont c("BP", "MF", "CC", "KEGG", "WikiPathways", "Reactome", "DO").
#' @details Use \code{ssim_search()} to search a set of compounds with similar
#' structure for target compound while length of \code{compound} is 1 and
#' \code{use.ssim_search} is TRUE.
#' @seealso \code{\link{ssim_search}}, \code{\link{stitch_network}}
#' @return list contained weights of each edge, annotations of each node
#' result of disease-related genes ORA analysis and result of enrichment
#' analysis as following:
#'
#' the tibble \code{edges} with fields:
#' * \code{node1, node2} protein- or compound-nodes in the network
#' * \code{type} type of interaction edges in network, e.g., cpi
#' (compound-protein), cci (compound-compound) and ppi (protein-protein)
#' * \code{score} combined score
#' * \code{experiment} experimental score
#' * \code{neighborhood} gene neighborhood score
#' * \code{fusion} gene fusion score
#' * \code{prediction} prediction score
#' * \code{cooccurence} cooccurence score
#' * \code{coexpression} coexpression score
#' * \code{database} database score
#' * \code{textmining} textmining score
#' * \code{evidence} more convincing interaction, the rank are: 1) experiment,
#' 2) neighborhood/fusion, 3) prediction/cooccurence/coexpression,
#' 4) database/textmining
#'
#' the tibble \code{nodes} with fields:
#' * \code{id} STRING id (compound) or STTICH id (protein)
#' * \code{type} type of nodes in network, protein or metabolite
#' * \code{name} common synonyms of compound or protein
#' * \code{external} external id of compound (PubChem cid) or protein
#' (entrez gene id)
#' * \code{description} description of compound or protein
#' @references Huimin Zheng, et al. (2022). In silico method to maximise the
#' biological potential of understudied metabolomic biomarkers: a study in
#' pre-eclampsia.
#' @export
medam_ssim_metabolites <- function(medam,
                                   compound,
                                   score = 700,
                                   search_ssim = TRUE,
                                   drg = NULL,
                                   universe = 21306,
                                   padj = "BH",
                                   pcutoff = 0.05,
                                   ont = c("BP", "MF", "CC", "KEGG",
                                           "WikiPathways", "Reactome", "DO")) {
  c2cid <- compound2cid(medam, compound)
  if (nrow(c2cid) > 0) {
    if (search_ssim) {
      c2ssimcid <- ssim_search(medam, unique(c2cid$cid))
      stitchid <- pull(c2ssimcid, stitch)
    } else {
      stitchid <- pull(c2cid, stitch)
    }
  } else {
    stitchid <- c()
  }
  if (length(stitchid) > 0) {
    res <- stitch_network(medam, stitchid, score = score)
    res$drg <- drgene_ora(res$nodes, drg, universe)
    res$enrichment <- enrichment_analysis(res$nodes, padj, pcutoff, ont)
  } else {
    res <- NULL
  }
  return(res)
}