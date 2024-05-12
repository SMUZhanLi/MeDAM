#' @title STITCH experimental CPI
#' @description STITCH experimental metabilite-protein interactions.
#' @param medam A Pool object connected MeDAM.db.
#' @param stitchid STITCH id, e.g. CIDs00002244.
#' @param score Threshold of experimental significance(default: 900).
#' @details The synonyms of compound can be converted to PubChem cid using
#' \code{\link{compound2cid}}.
#' @return the tibble with fields:
#' * \code{stitchv5} query metabilite(STITCH id in version V5)
#' * \code{experiment} experimental-interaction score
#' * \code{stringv10p5} experimental-interaction protein(STRING id in V10.5)
#' * \code{stringv12} covert STRING id from V10.5 to V12
#' * \code{description} description for protein
#' @export
stitch_expt_cpi <- function(medam, stitchid, score = 900) {
  ecpi <- medam |>
    dbquery("experimentcpi",
            node1 %in% !!stitchid,
            experiment >= !!score)
  stringid <- pull(ecpi, node2)
  eprotinfo <- medam |>
    dbquery("ecpiproteininfo", stringv10p5 %in% !!stringid)
  ecpi <- ecpi |>
    left_join(eprotinfo, by = c("node2" = "stringv10p5")) |>
    rename(stitchv5 = 1) |>
    select(-node2)
  return(ecpi)
}

#' @title STRING network
#' @description STRING interactions network.
#' @param medam A Pool object connected MeDAM.db.
#' @param stringid STRING id, e.g. ENSP00000269305.
#' @param score Threshold of significance to include a interaction, a number
#' between 0 and 1000 (default: 700).
#' @importFrom dplyr rename if_else bind_rows across where
#' @importFrom tibble add_column
#' @details The synonyms of protein can be converted to STRING id using
#' \code{\link{protein2stringid}}.
#' @return list contained weights of each edge and the annotations of each node
#' as following:
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
#' * \code{annotion} description of compound or protein
#' @export
string_network <- function(medam, stringid, score = 700) {
  ppi <- medam |>
    dbquery("stringv12ppi",
            node1 %in% !!stringid | node2 %in% !!stringid,
            score >= !!score)
  pnode1 <- pull(ppi, node1)
  pnode2 <- pull(ppi, node2)
  allprot <- unique(c(pnode1, pnode2))
  newprot <- allprot[!allprot %in% stringid]
  ppi2 <- medam |>
    dbquery("stringv12ppi",
            node1 %in% !!newprot,
            node2 %in% !!newprot,
            score >= !!score)
  ppi <- rbind(ppi, ppi2) |>
    add_column(prediction = NA, .before = "cooccurence") |>
    mutate(type = "ppi",
           across(where(is.integer), \(x) if_else(is.na(x), 0, x)))

  protinfo <- medam |>
    dbquery("stringv12proteininfo", stringv12 %in% allprot) |>
    rename(id = 1, external = 3) |>
    mutate(type = "protein",
           network = if_else(id %in% stringid, "input", "output"))
  res <- list(edges = ppi, nodes = protinfo)
  return(res)
}

#' @title STITCH network
#' @description STITCH interactions network
#' @param medam a Pool object connected MeDAM.db
#' @param stitchid STITCH id, e.g. CIDs00002244
#' @param score threshold of significance to include a interaction, a number
#' between 0 and 1000 (default: 700).
#' @details The synonyms of compound can be converted to PubChem cid using
#' \code{\link{compound2cid}}.
#' @return list contained weights of each edge and annotations of each node
#' as following:
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
#' * \code{annotion} description of compound or protein
#' @export
stitch_network <- function(medam, stitchid, score = 700) {
  cpi <- medam |>
    dbquery("stitchv5cpi",
            node1 %in% !!stitchid,
            score >= !!score) |>
    mutate(type = "cpi")
  stitchid <- cpi |> pull(node1) |> unique()
  cci <- medam |>
    dbquery("stitchv5cci",
            node1 %in% !!stitchid,
            node2 %in% !!stitchid,
            score >= !!score) |>
    mutate(type = "cci")

  stringid <- cpi |> pull(node2) |> unique()
  ppi <- medam |>
    dbquery("stringv10p5ppi",
            node1 %in% !!stringid,
            node2 %in% !!stringid,
            score >= !!score) |>
    mutate(type = "ppi")

  cheminfo <- medam |>
    dbquery("stitchv5compoundinfo", stitchv5 %in% !!stitchid) |>
    rename(id = 1, external = 3) |>
    mutate(type = "metabolite", network = "input")

  protinfo <- medam |>
    dbquery("stringv10p5proteininfo", stringv10p5 %in% !!stringid) |>
    rename(id = 1, external = 3) |>
    mutate(type = "protein", network = "output")

  edges <- bind_rows(cpi, cci, ppi) |>
    select("node1", "node2", "experiment", "fusion", "neighborhood",
           "prediction", "cooccurence", "coexpression", "database",
           "textmining", "score", "evidence", "type") |>
    mutate(across(where(is.integer), \(x) if_else(is.na(x), 0, x)))
  nodes <- bind_rows(cheminfo, protinfo)
  res <- list(edges = edges, nodes = nodes)
  return(res)
}