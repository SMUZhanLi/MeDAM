#' @title MeDAM batch search
#' @description Batch search for potential metabolite-disease associations using
#' MeDAM pipeline. Ref: fig 1.C in Huimin Zheng, et al. (2022).
#' @param medam A Pool object connected MeDAM.db.
#' @param abundance Abundance table, a matrix or data frame in which columns are
#' metabolites and rows ar samples.
#' @param disease Aliases of disease (ignored case),  e.g. "pre-eclampsia".
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
#' @param universe Count of universal background protein-coding genes.
#' For example, background of homo sapiens is 21306.
#' @param wgcna Result of WGCNA analysis for metabolite abundance. A dataframe
#' contains metabolite and module columns.
#' @param wgcna_params Define the parameters in WGCNA analysis when the param
#' \code{wgcna} is NULL and use \code{wgcna_analysis()} in medam batch.
#' @importFrom stats setNames
#' @importFrom utils modifyList
#' @return a list contained:
#' * \code{daa} Result of the differential abundance analysis for metabolomics.
#' See \code{\link{diff_metabolites}} and \code{\link{wgcna_analysis}}.
#' * \code{network} Interaction network. See \code{\link{string_network}} and
#' \code{\link{stitch_network}}.
#' * \code{drgora} Disease-related genes ORA. See \code{\link{drgene_ora}}
#' * \code{disease} Common name, DOID and related (entrez) genes of disease.
#' * \code{enrichment} NULL by default. Use \code{\link{medam_add_enrichment}}
#' to excute enriment analysis and add result into `$enrichment`
#' @references Huimin Zheng, et al. (2022). In silico method to maximise the
#' biological potential of understudied metabolomic biomarkers: a study in
#' pre-eclampsia.
#' @export

medam_batch <- function(medam,
                        abundance,
                        disease = NULL,
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
                        universe = 21306,
                        wgcna = NULL,
                        wgcna_params = list(abundance = abundance)) {
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
  if (!is.null(wgcna)) {
    wgcna <- setNames(wgcna, c("metabolite", "module"))
  } else {
    def_wgcna_params <- list(networkType = "signed",
                             corMethod = "spearman",
                             minModuleSize = 5,
                             cutHeight = 0.25,
                             nthreads = 10)
    def_wgcna_params <- modifyList(def_wgcna_params, wgcna_params)
    wgcna <- do.call(wgcna_analysis, def_wgcna_params)
  }
  daa <- daa |>
    left_join(wgcna, by = "metabolite")
  c2cid <- compound2cid(medam, colnames(abundance))
  daa <- daa |>
    left_join(c2cid, by = c("metabolite" = "compound"))
  tgtp_nw <- medam_tgtp_batch(medam, daa, score, ecpi_score)
  ssimm_nw <- medam_ssimm_batch(medam, daa, score)
  coabm_nw <- medam_coabm_batch(medam, daa, score)
  all_network <- get_all_network(daa, tgtp_nw, ssimm_nw, coabm_nw)
  if (grepl("DOID:\\d+", disease)) {
    doid <- disease
  } else {
    doid <- disease2doid(medam, disease) |> pull(doid)
  }
  drg <- drgene_search(medam, doid)
  eg <- pull(drg, ENTREZID)
  tgtp_drgora <- drgora_batch(tgtp_nw, "tgt_proteins_", eg, universe)
  ssimm_drgora <- drgora_batch(ssimm_nw, "ss_metabolites_", eg, universe)
  coabm_drgora <- drgora_batch(coabm_nw, "co_metabolites_", eg, universe)
  all_drgora <- tgtp_drgora |>
    left_join(ssimm_drgora, by = "metabolite") |>
    left_join(coabm_drgora, by = "metabolite")
  disease <- list(name = disease, doid = doid, drg = drg)
  res <- list(daa = daa, network = all_network, drgora = all_drgora,
              disease = disease, enrichment = NULL)
  return(res)
}

#' @title MeDAM batch search
#' @param medam A Pool object connected MeDAM.db.
#' @param metabolites A biomarker metabolites list.
#' @param disease Aliases of disease (ignored case),  e.g. "pre-eclampsia".
#' @param score Threshold of significance to include a interaction, a number
#' between 0 and 1000 (default: 700).
#' @param ecpi_score Threshold of significance for metabilite-protein
#' experimental interactions (from 0 to 1000, default: 800). See
#' \code{stitch_expt_cpi()}..
#' @param universe Count of universal background protein-coding genes.
#' For example, background of homo sapiens is 21306.
#' @param wgcna Result of WGCNA analysis for metabolite abundance. A dataframe
#' contains metabolite and module columns.
#' @param wgcna_params Define the parameters in WGCNA analysis when the
#' \code{wgcna} is NULL and use \code{wgcna_analysis()} in medam batch.
#' @rdname medam_batch
#' @export
medam_batch_manual <- function(medam,
                               metabolites,
                               disease = NULL,
                               score = 700,
                               ecpi_score = 800,
                               universe = 21306,
                               wgcna = NULL,
                               wgcna_params = list(abundance = NULL)) {
  def_wgcna_params <- list(networkType = "signed",
                           corMethod = "spearman",
                           minModuleSize = 5,
                           cutHeight = 0.25,
                           nthreads = 10)
  def_wgcna_params <- modifyList(def_wgcna_params, wgcna_params)
  daa <- tibble(metabolite = metabolites, significant = 1)
  if (is.null(wgcna) || is.null(def_wgcna_params$abundance)) {
    if (!is.null(wgcna)) {
      wgcna <- setNames(wgcna, c("metabolite", "module"))
    } else {
      wgcna <- do.call(wgcna_analysis, def_wgcna_params)
    }
    daa <- wgcna |>
      left_join(daa, by = "metabolite") |>
      mutate(significant = if_else(is.na(significant), 0, 1))
  }
  daa <- daa |>
    left_join(compound2cid(medam, daa$metabolite),
              by = c("metabolite" = "compound"))
  tgtp_nw <- medam_tgtp_batch(medam, daa, score, ecpi_score)
  ssimm_nw <- medam_ssimm_batch(medam, daa, score)
  if (!is.null(wgcna)) {
    coabm_nw <- medam_coabm_batch(medam, daa, score)
  } else {
    coabm_nw <- NULL
  }
  all_network <- get_all_network(daa, tgtp_nw, ssimm_nw, coabm_nw)
  if (grepl("DOID:\\d+", disease)) {
    doid <- disease
  } else {
    doid <- pull(disease2doid(medam, disease), doid)
  }
  drg <- drgene_search(medam, doid)
  eg <- pull(drg, ENTREZID)
  tgtp_drgora <- drgora_batch(tgtp_nw, "tgt_proteins_", eg, universe)
  ssimm_drgora <- drgora_batch(ssimm_nw, "ss_metabolites_", eg, universe)
  if (!is.null(wgcna)) {
    coabm_drgora <- drgora_batch(coabm_nw, "co_metabolites_", eg, universe)
  } else {
    coabm_drgora <- tibble(metabolite = daa$metabolite, ngene = 0, indrg = NA,
                           pvalue = NA, padj = NA, overlap = NA) |>
      rename_at(vars(-metabolite), ~ paste0("co_metabolites_", .x))
  }
  all_drgora <- all_drgora <- tgtp_drgora |>
    left_join(ssimm_drgora, by = "metabolite") |>
    left_join(coabm_drgora, by = "metabolite")
  disease <- list(name = disease, doid = doid, drg = drg)
  res <- list(daa = daa, network = all_network, drgora = all_drgora,
              disease = disease, enrichment = NULL)
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

# target proteins batch
#' @importFrom dplyr right_join
#' @keywords internal
medam_tgtp_batch <- function(medam, dat, score, ecpi_score) {
  sdat <- dat |> filter(significant == 1)
  c2tgtp <- sdat |>
    pull(stitch) |>
    stitch_expt_cpi(medam, stitchid = _, ecpi_score) |>
    rename(stitch = 1) |>
    right_join(sdat, by = "stitch", relationship = "many-to-many") |>
    filter(!is.na(stringv12))
  stringid <- c2tgtp |> pull(stringv12) |> unique()
  nw <- string_network(medam, stringid, score)
  c2tgtp <- split(x = c2tgtp$stringv12, f = c2tgtp$metabolite)
  c2nw <- lapply(c2tgtp, function(tgtpstring) {
    subppi1 <- nw$edges |>
      filter(node1 %in% !!tgtpstring | node2 %in% !!tgtpstring)
    pnode1 <- subppi1 |> pull(node1)
    pnode2 <- subppi1 |> pull(node2)
    allnode <- c(pnode1, pnode2) |> unique()
    nottgtp <- allnode[!allnode %in% tgtpstring]
    subppi2 <- nw$edges |>
      filter(node1 %in% !!nottgtp, node2 %in% !!nottgtp)
    edges <- bind_rows(subppi1, subppi2)
    nodes <- nw$nodes |>
      filter(id %in% !!allnode) |>
      mutate(network = if_else(id %in% tgtpstring, "input", "output"))
    list(edges = edges, nodes = nodes)
  })
  res <- vector("list", nrow(sdat))
  names(res) <- sdat |> pull(metabolite)
  for (x in names(c2nw)) {
    res[[x]] <- c2nw[[x]]
  }
  return(res)
}

## structural similar metabolites metabolites batch
#' @keywords internal
medam_ssimm_batch <- function(medam, dat, score) {
  sdat <- dat |> filter(significant == 1)
  c2ssimm <- sdat |>
    pull(cid) |>
    ssim_search(medam, cid = _) |>
    rename(ssimstitch = 3) |>
    right_join(sdat, by = "cid", relationship = "many-to-many") |>
    select(metabolite, cid, ssimcid, ssimstitch) |>
    filter(!is.na(ssimstitch))
  stitchid <- c2ssimm |> pull(ssimstitch) |> unique()
  nw <- stitch_network(medam, stitchid, score)
  c2ssimm <- split(x = c2ssimm$ssimstitch, f = c2ssimm$metabolite)
  c2nw <- lapply(c2ssimm, function(ssimstitch) {
    subcpi <- nw$edges |> filter(node1 %in% !!ssimstitch, type == "cpi")
    cnode <- subcpi |> pull(node1) |> unique()
    subcci <- nw$edges |>
      filter(node1 %in% !!cnode, node2 %in% !!cnode, type == "cci")
    pnode <- subcpi |> pull(node2) |> unique()
    subppi <- nw$edges |>
      filter(node1 %in% !!pnode, node2 %in% !!pnode, type == "ppi")
    edges <- bind_rows(subcpi, subcci, subppi)
    allnode <- c(cnode, pnode)
    nodes <- nw$nodes |> filter(id %in% !!allnode)
    list(edges = edges, nodes = nodes)
  })
  res <- vector("list", nrow(sdat))
  names(res) <- sdat |> pull(metabolite)
  for (x in names(c2nw)) {
    res[[x]] <- c2nw[[x]]
  }
  return(res)
}

## co-abundance metabolites batch
#' @keywords internal
medam_coabm_batch <- function(medam, dat, score) {
  sdat <- dat |> filter(significant == 1)
  mod <- sdat |> pull(module) |> unique()
  stitchid <- dat |>
    filter(module %in% mod) |>
    pull(stitch) |>
    na.omit() |>
    unique()
  nw <- stitch_network(medam, stitchid, score)
  mod2coabm <- dat |>
    filter(!is.na(stitch), module %in% mod) |>
    (\(d) split(x = d$stitch, f = d$module))()
  mod2nw <- lapply(mod2coabm, function(coabm) {
    subcpi <- nw$edges |> filter(node1 %in% !!coabm, type == "cpi")
    cnode <- subcpi |> pull(node1) |> unique()
    subcci <- nw$edges |>
      filter(node1 %in% !!cnode, node2 %in% !!cnode, type == "cci")
    pnode <- subcpi |> pull(node2) |> unique()
    subppi <- nw$edges |>
      filter(node1 %in% !!pnode, node2 %in% !!pnode, type == "ppi")
    edges <- bind_rows(subcpi, subcci, subppi)
    allnode <- c(cnode, pnode)
    nodes <- nw$nodes |> filter(id %in% !!allnode)
    list(edges = edges, nodes = nodes)
  })
  res <- sdat |>
    pull(module) |>
    lapply(function(m) mod2nw[[m]]) 
  names(res) <- sdat |> pull(metabolite)
  return(res)
}

## Combine all network of three branch in MeDAM batch search
#' @keywords internal
get_all_network <- function(dat, tgtp_nw, ssimm_nw, coabm_nw) {
  diffmetabo <- dat |> filter(significant == 1) |> pull(metabolite)
  all_network <- lapply(diffmetabo, function(x) {
    list(target_proteins = tgtp_nw[[x]],
         ssim_metabolites = ssimm_nw[[x]],
         coab_metabolites = coabm_nw[[x]])
  }) |>
  setNames(diffmetabo)
  return(all_network)
}