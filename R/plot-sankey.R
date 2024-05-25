# transf_tgtp_data <- function(medamres) {
#   tgtpedges <- medamres$edges
#   tgtpnodes <- medamres$nodes
#   tgtp <- tgtpnodes |> filter(network == "input") |> pull(id)
#   tgtpedges2 <- tgtpedges |>
#     select(2, 1, 3:12) |>
#     rename(node1 = 1, node2 = 2) |>
#     bind_rows(tgtpedges) |>
#     filter(node1 %in% tgtp) |>
#     add_row(node1 = tgtp,
#             node2 = tgtp,
#             evidence = "itself") |>
#     mutate(node1 = paste0(node1, "-tgtp"), type = "ppi")
#   tgtpnodes <- tgtpnodes |> mutate(network = "output")
#   tgtpnodes2 <- tgtpnodes |>
#     filter(id %in% tgtp) |>
#     mutate(id = paste0(id, "-tgtp"),
#            name = paste0(name, "[tgtp]"),
#            network = "input") |>
#     bind_rows(tgtpnodes)
#   medamres$edges <- tgtpedges2
#   medamres$nodes <- tgtpnodes2
#   return(medamres)
# }

# ## Get data for sankey plot
# #' @importFrom dplyr add_row
# #' @keywords internal
# get_sankey_data <- function(medamres, metabolite, branch) {
#   edges <- medamres$edges[[metabolite]][[branch]]
#   nodes <- medamres$nodes[[metabolite]][[branch]]
#   enrichres <- medamres$enrichres[[metabolite]]
#   drg <- enrichres$geneID |>
#     strsplit(split = "/") |>
#     unlist()
#   drg_nodes <- nodes |>
#     filter(type == "protein", name %in% drg)
#   edges <- edges |>
#     filter(node2 %in% drg_nodes$id)
#   allid <- c(edges$node1, edges$node2) |> unique()
#   nodes <- nodes |> filter(id %in% allid)

#   desc_col <- c("geneID", "p.adjust", "GeneRatio", "BgRatio", "Count")
#   new_desc <- paste0(desc_col, ": ", enrichres[, desc_col])
#   new_desc <- paste0(new_desc, collapse = "; ")
#   nodes <- nodes |>
#     mutate(level = if_else(network == "input", 0, 1)) |>
#     add_row(id = enrichres$ID,
#             name = enrichres$Description,
#             description = new_desc,
#             type = "pathway",
#             level = 2) |>
#     mutate(branch = branch)
#   edges <- edges |>
#     add_row(node1 = drg_nodes$id,
#             node2 = enrichres$ID,
#             evidence = enrichres$Description) |>
#     mutate(branch = branch)
#   res <- list(edges = edges, nodes = nodes)
# }
