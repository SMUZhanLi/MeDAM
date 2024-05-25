#' @keywords internal
transf_tgtp_data <- function(tgtpnw) {
  tgtpedges <- tgtpnw$edges
  tgtpnodes <- tgtpnw$nodes
  tgtp <- tgtpnodes |> filter(network == "input") |> pull(id)
  tgtpedges2 <- tgtpedges |>
    select(2, 1, 3:12) |>
    rename(node1 = 1, node2 = 2) |>
    bind_rows(tgtpedges) |>
    filter(node1 %in% tgtp) |>
    add_row(node1 = tgtp,
            node2 = tgtp,
            evidence = "itself") |>
    mutate(node1 = paste0(node1, "-tgtp"), type = "ppi")
  tgtpnodes <- tgtpnodes |> mutate(network = "output")
  tgtpnodes2 <- tgtpnodes |>
    filter(id %in% tgtp) |>
    mutate(id = paste0(id, "-tgtp"),
           name = paste0(name, "[tgtp]"),
           network = "input") |>
    bind_rows(tgtpnodes)
  new_tgtpnw <- list(edges = tgtpedges2, nodes = tgtpnodes2)
  return(new_tgtpnw)
}

## Get data for sankey plot
#' @importFrom dplyr add_row
#' @keywords internal
get_sankey_data <- function(medamnw, enrichda, branch) {
  edges <- medamnw$edges
  nodes <- medamnw$nodes
  drg <- enrichda$geneID |>
    strsplit(split = "/") |>
    unlist()
  drg_nodes <- nodes |>
    filter(type == "protein", name %in% drg)
  edges <- edges |>
    filter(node2 %in% drg_nodes$id)
  allid <- c(edges$node1, edges$node2) |> unique()
  nodes <- nodes |> filter(id %in% allid)
  desc_col <- c("geneID", "p.adjust", "GeneRatio", "BgRatio", "Count")
  new_desc <- paste0(desc_col, ": ", enrichda[, desc_col])
  new_desc <- paste0(new_desc, collapse = "; ")
  nodes <- nodes |>
    mutate(level = if_else(network == "input", 0, 1)) |>
    add_row(id = enrichda$ID,
            name = enrichda$Description,
            description = new_desc,
            type = "pathway",
            level = 2) |>
    mutate(branch = branch)
  edges <- edges |>
    add_row(node1 = drg_nodes$id,
            node2 = enrichda$ID,
            evidence = enrichda$Description) |>
    mutate(branch = branch)
  new_medamnw <- list(edges = edges, nodes = nodes)
  return(new_medamnw)
}

#' @title Sankey plot
#' @description Show interactions between: (1) the potential proteins and the
#' query nodes in MeDAM branch (e.g. the target proteins, the structural similar
#' metabolites and the co-abundant metabolites); (2) the potential proteins and
#' functional terms.
#' @param medamnw The medam network contained edges and nodes data.
#' @param enrichres The result of enrichment analysis.
#' @param branch One of "target_proteins", "ssim_metabolites" and
#' "coab_metabolites"
#' @param termsda The disease related functional items or pathways.
#' @param max_metabo_len A numeric value sets wrap length for labels of
#' metabolites. It will wrap names longer that 30 characters by default.
#' @param max_desc_len A numeric value sets wrap length for labels of items.
#' It will wrap names longer that 30 characters by default.
#' @importFrom aplot insert_left insert_right
#' @export
medam_sankey_plot <- function(medamnw,
                              enrichres,
                              branch,
                              termsda,
                              max_metabo_len = 30,
                              max_desc_len = 30) {
  medam_branch <- c("target_proteins", "ssim_metabolites", "coab_metabolites")
  branch <- match.arg(branch, medam_branch)
  if (branch == "target_proteins") {
    edges.type <- "ppi"
    branch <- "target proteins"
    medamnw <- transf_tgtp_data(medamnw)
  } else {
    edges.type <- "cpi"
    if (branch == "ssim_metabolites") {
      branch <- "structural similar\nmetabolites"
    } else {
      branch <- "co-abundant\nmetabolites"
    }
  }
  medamnw$edges <- medamnw$edges |>
    filter(type == edges.type)
  medamnw$nodes <- medamnw$nodes |>
    mutate(name = metabo_wrap(name, max_metabo_len))
  enrichda <- lapply(enrichres, as.data.frame)
  sankyda <- lapply(seq_len(nrow(termsda)), function(n) {
    ont <- termsda[n, "ont"]
    i <- termsda[n, "index"]
    eda <- enrichda[[ont]][i, ]
    eda$Description <- desc_wrap(eda$Description, max_desc_len)
    get_sanky_data(medamnw, eda, branch)
  })
  new_edges <- lapply(sankyda, `[[`, "edges") |>
    do.call("rbind", args = _) |>
    distinct(node1, node2, .keep_all = TRUE)

  new_nodes <- lapply(sankyda, `[[`, "nodes") |>
    do.call("rbind", args = _) |>
    distinct(id, .keep_all = TRUE)
  freq <- split(new_nodes, new_nodes$level) |> sapply(nrow)
  max_freq <- max(freq)
  degreeda <- c(new_edges$node1, new_edges$node2) |>
    table() |>
    as.data.frame() |>
    setNames(c("id", "degree")) |>
    filter(grepl("^ENSP|^CIDs", id))
  term_nodes <- new_nodes |>
    filter(level == 2) |>
    select(id, description)
  term_nodes_desc <- term_nodes |> pull(description) |> strsplit(";|:")
  degreeda <- term_nodes |>
    mutate(count = as.integer(unlist(lapply(term_nodes_desc, `[`, 10))),
           p.adjust = as.numeric(unlist(lapply(term_nodes_desc, `[`, 4)))) |>
    arrange(count, desc(p.adjust)) |>
    mutate(degree = seq_len(n())) |>
    select(c("id", "degree")) |>
    bind_rows(degreeda)

  node_type <- c(branch, "interaction protein", "functional term")
  new_nodes <- new_nodes |>
    left_join(degreeda, by = "id") |>
    group_by(level) |>
    arrange(degree) |>
    mutate(y = rm_head_tail(seq(1, max_freq, length.out = if_else(n() == max_freq, max_freq + 2, n() + 2)), 1)) |>
    ungroup() |>
    mutate(x = level, `Node type` = node_type[level + 1])
  index1 <- match(new_edges$node1, new_nodes$id)
  index2 <- match(new_edges$node2, new_nodes$id)
  new_edges <- new_edges |>
    mutate(x = new_nodes$x[index1], y = new_nodes$y[index1],
           xend = new_nodes$x[index2], yend = new_nodes$y[index2],
           circular = FALSE, edge.id = seq_len(n()),
           evidence = factor(evidence, levels = link_type_levels))
  query2prot_nodes <- new_nodes |> filter(level != 2)
  query2prot_edges <- new_edges |>
    filter(node1 %in% query2prot_nodes$id,
           node2 %in% query2prot_nodes$id)
  p1 <- internal_sankey_plot(query2prot_edges, query2prot_nodes, branch, TRUE,
                             max(freq) + 1)
  prot2pw_nodes <- new_nodes |> filter(level != 0)
  prot2pw_edges <- new_edges |>
    filter(node1 %in% prot2pw_nodes$id,
           node2 %in% prot2pw_nodes$id)
  p2 <- internal_sankey_plot(prot2pw_edges, prot2pw_nodes, branch, FALSE,
                             max(freq) + 1)
  p12 <- insert_left(p2[[1]], p1, 0.8) |> insert_right(p2[[2]], 0.3)
  return(p12)
}

## Add left of right network plot in sanky plot
#' @importFrom scales hue_pal
#' @importFrom ggplot2 scale_fill_manual sec_axis
#' @importFrom ggraph scale_edge_colour_manual
#' @keywords internal
internal_sankey_plot <- function(new_edges, new_nodes, branch, left = TRUE, limit_ymax) {
  point_shapes <- c(23, 22, 21) |>
    setNames(c(branch, "interaction protein", "functional term"))
  point_cols <- rev(hue_pal()(3)) |>
    setNames(c(branch, "interaction protein", "functional term"))
  string_edges_cols <- c(
    "#ADF2AD", "#4075C1", "#ffff33", "#c2a5cf", "#ff7f00",
    "#810F7C", "#F4AF6E", "#0fb451", "#E58585"
  )
  string_edges_levels <- c(
    "itself", "experiment", "neighborhood", "fusion",
    "prediction", "cooccurence", "coexpression", "database", "textmining"
  )
  term_edges_levels <- new_nodes |> filter(level == 2) |> pull(name)
  term_edges_cols <- hue_pal()(length(term_edges_levels))
  link_type_levels <- c(string_edges_levels, term_edges_levels)
  link_type_cols <- c(string_edges_cols, term_edges_cols) |>
    setNames(link_type_levels)
  psankey <- ggplot(new_nodes) +
    geom_edge_diagonal(aes(x = x, y = y, xend = xend, yend = yend,
                           circular = circular, edge_colour = evidence),
                       data = new_edges, edge_width = 0.5, flipped = TRUE,
                       edge_alpha = ifelse(left, 1, 0.5)) +
    scale_edge_colour_manual(values = link_type_cols) +
    geom_point(aes(x = x, y = y, shape = `Node type`, color = `Node type`,
                   fill = `Node type`),
               data = new_nodes, size = 3, show.legend = FALSE) +
    scale_shape_manual(values = point_shapes) +
    scale_fill_manual(values = point_cols) +
    scale_color_manual(values = point_cols)

  if (left) {
    psankey <- psankey +
      scale_y_continuous(
        expand = c(0.02, 0.02),
        # limits = c(0, limit_ymax),
        breaks = new_nodes$y[new_nodes$x == 0],
        labels = new_nodes$name[new_nodes$x == 0],
        sec.axis = sec_axis(transform = ~.,
                            breaks = new_nodes$y[new_nodes$x == 1],
                            labels = new_nodes$name[new_nodes$x == 1])) +
      scale_x_continuous(
          expand = c(0.02, 0.02),
          breaks = c(0, 1),
          labels = c(branch, "interaction proteins")) +
      guides(edge_color = guide_legend(
        title = paste0(
          ifelse(branch == "target proteins", "proteins", "metabolotes"),
          "-proteins\ninteractions"))) +
      theme(axis.text.y.right = element_text(hjust = 0.5, colour = "black",
                                             size = 12),
            axis.text.y = element_text(colour = "black", size = 12))
  } else {
    psankey <- psankey +
      scale_y_continuous(expand = c(0.02, 0.02), position = "right") +
      scale_x_continuous(expand = c(0.02, 0.02), breaks = 2,
                         labels = "functional terms") +
      theme(legend.position = "none",
            axis.text.y = element_blank())
  }
  psankey <- psankey +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_rect(fill = "white", colour = "white"),
          panel.border = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_text(colour = "black", size = 12, angle = 30,
                                     vjust = 1, hjust = 1),
          axis.title = element_blank(),
          legend.title = element_text(colour = "black", size = 13),
          legend.text = element_text(colour = "black", size = 12),
          plot.margin = margin(),
          plot.title = element_text(size = 16, hjust = 0.5))
  if (!left) {
    term_nodes <- new_nodes |> filter(level == 2)
    term_nodes_desc <- term_nodes |> pull(description) |> strsplit(";|:")
    GeneRatio <- lapply(term_nodes_desc, `[`, 6) |>
      unlist() |>
      strsplit("/") |>
      lapply(function(x) {
        x <- as.numeric(x)
        x[1] / x[2]
      }) |>
      unlist()
    term_nodes <- term_nodes |>
      mutate(count = as.integer(unlist(lapply(term_nodes_desc, `[`, 10))),
             GeneRatio = GeneRatio,
             p.adjust = as.numeric(unlist(lapply(term_nodes_desc, `[`, 4))))
    psankey <- list(psankey, internal_dotplot(term_nodes))
  }
  return(psankey)
}

## Add dotplot in sanky plot
#' @importFrom enrichplot set_enrichplot_color
#' @importFrom ggplot2 theme_bw margin
#' @keywords internal
internal_dotplot <- function(dat) {
  pdot <- ggplot(dat, aes(x = GeneRatio, y = y, size = GeneRatio,
                          fill = p.adjust)) +
    geom_point(shape = 21) +
    set_enrichplot_color(type = "fill", name = "p.adjust") +
    scale_y_continuous(expand = c(0.02, 0.02), breaks = dat$y,
                       labels = dat$name, position = "right") +
    scale_size(range = c(3, 8)) +
    scale_x_continuous(expand = c(0.02, 0.02)) +
    theme_bw() +
    theme(axis.text.x = element_text(colour = "black", size = 12, angle = 90),
          axis.text.y.right = element_text(colour = "black", size = 12),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          legend.title = element_text(colour = "black", size = 12),
          legend.text = element_text(colour = "black", size = 12),
          plot.margin = margin(),
          plot.title = element_text(size = 16, hjust = 0.5))
  return(pdot)
}