#' @title MeDAM batch search
#' @description Batch search for potential metabolite-disease associations using
#' MeDAM pipeline. Ref: fig 1.C in Huimin Zheng, et al. (2022).
#' @param medam A Pool object connected MeDAM.db.
#' @param disease Aliases of disease (ignored case),  e.g. "pre-eclampsia".
#' @param abundance Abundance table, a matrix or data frame in which columns are
#' metabolites and rows ar samples.
#' @param resp Response (i.e, group) to be modelled, a factor (same length as
#' \code{abundance} row number).
#' @param compared Define the condition A and condition B in "response" for
#' calculating fold change.
#' For example, set "compared" to c("case", "control"), then the fold change is
#' "case/control".
#' @param ortho FALSE (for PLS-DA case) or TRUE (for OPLS-DA).
#' @param log10L Whether the 'Abudnace Table' table needs to be
#' log10-transformed.
#' @param scaling Scaling methods for "abundance", one of "none", "center",
#' "pareto", and "standard" (default).
#' @param test Use t.test or wilcox.test for differential abundance analysis.
#' @param padj_method  P-values adjusted methods, one of "holm", "hochberg",
#' "hommel", "bonferroni", "BH", "BY", "fdr" and "none".
#' @param vip_cutoff Threshold of VIP score of the (O)PLS-DA model.
#' @param padj_cutoff Threshold of P-value adjusted for t.test (or wilcox.test).
#' @param auc_cutoff Threshold of auc of ROC curve.
#' @param log2fc_cutoff Threshold of absolute value of log2 fold change.
#' @param score Threshold of significance to include a interaction, a number
#' between 0 and 1000 (default: 700).
#' @param ecpi_score Threshold of significance for metabilite-protein
#' experimental interactions (from 0 to 1000, default: 800). See
#' \code{stitch_expt_cpi()}.
#' @param networkType Network type, one of "signed", "unsigned" and
#' "signed hybrid".
#' @param corMethod Correlation algorithm for network construction, e.g.,
#' "spearman" or "pearson".
#' @param minModuleSize Minimum size of module or cluster, default is 5.
#' @param cutHeight Maximum dissimilarity (i.e., 1-correlation) that qualifies
#' modules for merging.
#' @param universe Count of universal background protein-coding genes.
#' For example, background of homo sapiens is 21306.
#' @return a list contained:
#' * \code{daa} Result of the differential abundance analysis for metabolomics.
#' See \code{\link{diff_metabolites}} and \code{\link{wgcna_analysis}}.
#' * \code{network} Interaction network. See \code{\link{string_network}} and
#' \code{\link{stitch_network}}.
#' * \code{drgora} Disease-related genes ORA. See \code{\link{drgene_ora}}
#' * \code{disease} Common name, DOID and related (entrez) genes of disease.
#' @references Huimin Zheng, et al. (2022). In silico method to maximise the
#' biological potential of understudied metabolomic biomarkers: a study in
#' pre-eclampsia.
#' @export

medam_batch <- function(medam,
                        disease,
                        abundance,
                        resp,
                        compared,
                        ortho = FALSE,
                        log10L = FALSE,
                        scaling = "standard",
                        vip_cutoff = 1,
                        test = "t.test",
                        padj_method = "none",
                        padj_cutoff = 0.01,
                        auc_cutoff = 0.5,
                        log2fc_cutoff = 1.5,
                        score = 700,
                        ecpi_score = 800,
                        networkType = "signed",
                        corMethod = "spearman",
                        minModuleSize = 5,
                        cutHeight = 0.25,
                        universe = 21306) {
  daa <- diff_metabolites(abundance = abundance,
                          resp = resp,
                          compared = compared,
                          ortho = ortho,
                          log10L = log10L,
                          scaling = scaling, test = test,
                          padj_method = padj_method,
                          vip_cutoff = vip_cutoff,
                          padj_cutoff = padj_cutoff,
                          auc_cutoff = auc_cutoff,
                          log2fc_cutoff = log2fc_cutoff)
  wgcna <- wgcna_analysis(abundance, networkType, corMethod,
                          minModuleSize, cutHeight)
  daa <- daa |>
    left_join(wgcna, by = "metabolite")
  c2cid <- compound2cid(medam, colnames(abundance))
  daa <- daa |>
    left_join(c2cid, by = c("metabolite" = "compound"))

  tgtp_nw <- target_proteins_batch(medam, daa, score, ecpi_score)
  coabm_nw <- coab_metabolites_batch(medam, daa, score)
  ssimm_nw <- ssim_metabolites_batch(medam, daa, score)
  all_network <- list(target_proteins = tgtp_nw,
                      coab_metabolites = coabm_nw,
                      ssim_metabolites = ssimm_nw)
  if (grepl(disease, "DOID:\\d+")) {
    doid <- disease
  } else {
    doid <- disease2doid(medam, disease) |> pull(doid)
  }
  drg <- drgene_search(medam, doid)
  eg <- pull(drg, ENTREZID)
  tgtp_drgora <- drgora_batch(tgtp_nw, "tgt_proteins_", eg, universe)
  coabm_drgora <- drgora_batch(coabm_nw, "co_metabolites_", eg, universe)
  ssimm_drgora <- drgora_batch(ssimm_nw, "ss_metabolites_", eg, universe)
  all_drgora <- tgtp_drgora |>
    left_join(coabm_drgora, by = "metabolite") |>
    left_join(ssimm_drgora, by = "metabolite")
  disease <- list(name = disease, doid = doid, drg = drg)
  res <- list(daa = daa, network = all_network, drgora = all_drgora,
              disease = disease)
  return(res)
}

## target proteins batch
#' @importFrom dplyr group_map
#' @keywords internal
target_proteins_batch <- function(medam, dat, score, ecpi_score) {
  dat <- dat |> filter(significant == 1)
  stitchid <- dat |> pull(stitch)
  ecpi <- stitch_expt_cpi(medam, stitchid, ecpi_score)
  arm1_search <- ecpi |>
    group_by(stitchv5) |>
    group_map(~ list(key = .y, nw = string_network(medam, .x$stringv12, score)))
  key <- lapply(arm1_search, "[[", "key")
  nw <- lapply(arm1_search, "[[", "nw")
  names(nw) <- unlist(key)
  res <- lapply(stitchid, function(x) {
    nw[[x]]
  })
  names(res) <- pull(dat, metabolite)
  return(res)
}


## co-abundance metabolites batch
#' @keywords internal
coab_metabolites_batch <- function(medam, dat, score) {
  biomk2mod <- dat |>
    filter(significant == 1)
  modmk <- biomk2mod |> pull(module) |> unique()
  arm2b1_search <- dat |>
    filter(module %in% modmk) |>
    group_by(module) |>
    group_map(~ list(key = .y, nw = stitch_network(medam, .x$stitch, score)))
  key <- lapply(arm2b1_search, "[[", "key")
  nw <- lapply(arm2b1_search, "[[", "nw")
  names(nw) <- unlist(key)
  res <- pull(biomk2mod, module) |>
    lapply(function(m) {
      nw[[m]]
    })
  names(res) <- pull(biomk2mod, metabolite)
  return(res)
}

## structural similar metabolites metabolites batch
#' @keywords internal
ssim_metabolites_batch <- function(medam, dat, score) {
  dat <- dat |> filter(significant == 1)
  cid <- dat |> pull(cid)
  c2ssimcid <- ssim_search(medam, cid)
  arm2b2_search <- c2ssimcid |>
    group_by(cid) |>
    group_map(~ list(key = .y, nw = stitch_network(medam, .x$stitch, score)))
  key <- lapply(arm2b2_search, "[[", "key")
  nw <- lapply(arm2b2_search, "[[", "nw")
  names(nw) <- unlist(key)
  res <- lapply(cid, function(x) {
    nw[[x]]
  })
  names(res) <- dat |> pull(metabolite)
  return(res)
}

## drgora batch
#' @importFrom tibble rownames_to_column as_tibble
#' @importFrom dplyr rename_at vars
#' @keywords internal
drgora_batch <- function(network, prefix, drg, universe) {
  drgora <- lapply(network, function(x) {
    ora <- drgene_ora(x$nodes, drg, universe)
    as.data.frame(ora)
  })
  drgora <- do.call(rbind, drgora) |>
    rownames_to_column("metabolite") |>
    as_tibble() |>
    mutate(padj = p.adjust(pvalue, method = "BH")) |>
    select(1:4, 6, 5) |>
    rename_at(vars(-metabolite), ~ paste0(prefix, .x))
  return(drgora)
}