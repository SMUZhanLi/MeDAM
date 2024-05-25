#' @title Circle plot
#' @description Show statistical measure of metabolites which have
#' differentially abundance in the two groups.
#' @param batchres Result of MeDAM batch search.
#' @param sort.by Order the metabolites by selected variable, one of "vip",
#' "abs(log2fc)", "-log10(pvalue)" and "auc".
#' @param max_metabo_len A numeric value sets wrap length for metabolite labels.
#' It will wrap names longer that 30 characters by default.
#' @param start Offset of starting point from 12 o'clock in radians. Offset is
#' applied clockwise.
#' @param .offset Offset of each statistical measure in different circular ring.
#' @importFrom stats median
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot geom_point geom_text geom_tile scale_color_manual
#' scale_shape_manual scale_fill_gradient scale_size coord_polar labs guides
#' guide_legend scale_x_continuous scale_y_continuous expansion theme
#' element_blank element_text unit
#' @importFrom ggraph geom_edge_diagonal
#' @importFrom ggnewscale new_scale_color new_scale_fill
#' @export
medam_circle_plot <- function(batchres,
                              sort.by = "vip",
                              max_metabo_len = 30,
                              start = 0,
                              .offset = 0.1) {
  rad2deg <- start * 180 / pi
  measures <- c("auc", "abs(log2fc)", "-log10(pvalue)", "vip")
  sort.by <- match.arg(sort.by, measures)
  measures <- measures[measures != sort.by]
  pdat <- batchres$daa |>
    filter(significant == 1) |>
    left_join(batchres$drgora, by = "metabolite")
  g <- unique(pdat$fc_direct) |>
    strsplit("/") |>
    unlist()
  pdat <- pdat |>
    mutate(
      `abs(log2fc)` = abs(log2fc),
      `-log10(pvalue)` = -log10(pvalue),
      bar.h = abs(.data[[sort.by]])) |>
    arrange(bar.h * sign(log2fc)) |>
    mutate(group = if_else(log2fc > 0, g[1], g[2]),
           edge.id = seq(1, n()),
           edge.x = edge.id + max(sign(log2fc), 0),
           edge.y = 2,
           edge.yend = -1,
           circular = FALSE,
           label_len = unlist(lapply(strsplit(metabolite, ""), length)),
           label_size = 1 - (label_len %/% (max_metabo_len + 1)) / 5,
           label = metabo_wrap(metabolite, max_metabo_len),
           angle = 90 - rad2deg - (360 * (edge.x) / (n() + 2)),
           hjust = if_else(angle < - 90 - rad2deg, 1, 0),
           angle = if_else(angle < - 90 - rad2deg, angle + 180, angle),
           tgtp.y = if_else(is.na(tgt_proteins_padj), NA, 2 + .offset),
           ssimm.y = if_else(is.na(ss_metabolites_padj), NA, 2.5 + .offset * 2),
           coabm.y = if_else(is.na(co_metabolites_padj), NA, 3 + .offset * 3),
           tile1.y = 4 + .offset * 4,
           tile2.y = 5 + .offset * 5,
           tile3.y = 6 + .offset * 6,
           norm.bar.h = norm_min_max(bar.h) * 1.5 + 1,
           bar.y = 6.5 + .offset * 7 + norm.bar.h / 2
    ) |>
    group_by(group) |>
    mutate(edge.xend = median(edge.x)) |>
    ungroup()
  grdat <- pdat |>
    select(edge.xend, edge.yend, group) |>
    distinct()
  medamdat <- pdat |>
    select(edge.x, tgtp.y, ssimm.y, coabm.y) |>
    pivot_longer(!edge.x, names_to = "branch", values_to = "medam.y",
                 values_drop_na = TRUE) |>
    mutate(branch = factor(sub("\\.y$", "", branch),
                           levels = c("tgtp", "ssimm", "coabm")))
  gg <- ggplot(pdat) +
    geom_edge_diagonal(aes(x = edge.x, y = edge.y, xend = edge.xend,
                           yend = edge.yend, circular = circular),
                       data = NULL, color = "gray", linewidth = 0.5) +
    geom_point(aes(x = edge.xend, y = edge.yend, color = group),
               data = grdat, size = 20, show.legend = FALSE) +
    scale_color_manual(values = c("#00AED7", "#FD9347")) +
    new_scale_color() +
    geom_text(aes(x = edge.xend, y = edge.yend, label = group),
              data = grdat, size = 4) +
    geom_point(aes(x = edge.x, y = medam.y, shape = branch, color = branch),
               data = medamdat, size = 2.5) +
    scale_shape_manual(values = setNames(c(16, 15, 17),
                                         c("tgtp", "ssimm", "coabm"))) +
    scale_color_manual(values = setNames(c("#9bd1df", "#f0d091", "#f9aea8"),
                                         c("tgtp", "ssimm", "coabm"))) +
    geom_tile(aes(x = edge.x, y = tile1.y, fill = .data[[measures[1]]]),
              width = 1, height = 1) +
    scale_fill_gradient(low = "#f5fbfb", high = "#51b9b2") +
    new_scale_fill() +
    geom_tile(aes(x = edge.x, y = tile2.y, fill = .data[[measures[2]]]),
              width = 1, height = 1) +
    scale_fill_gradient(low = "#cbdeef", high = "#1f69da") +
    new_scale_fill() +
    geom_tile(aes(x = edge.x, y = tile3.y, fill = .data[[measures[3]]]),
              width = 1, height = 1) +
    scale_fill_gradient(low = "#e0e4f2", high = "#9175c5") +
    new_scale_fill() +
    geom_tile(aes(x = edge.x, y = bar.y, height = norm.bar.h, fill = bar.h),
              width = 0.9, alpha = 0.5) +
    scale_fill_gradient(low = "#fbe6d1", high = "#fe621f") +
    labs(fill = sort.by) +
    new_scale_fill() +
    geom_text(aes(x = edge.x, y = bar.y + norm.bar.h / 2 + 0.05, label = label,
                  angle = angle, hjust = hjust, size = label_size),
              fontface = "bold", show.legend = FALSE) +
    scale_size(range = c(2.5, 3)) +
    coord_polar(start = start) +
    scale_x_continuous(expand = expansion(mult = c(0, 0)),
                       limits = c(0, max(pdat$edge.x) + 1)) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    guides(color = guide_legend(order = 1),
           shape = guide_legend(order = 1, override.aes = list(size = 3))) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.background = element_rect(fill = "white", colour = "white"),
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      legend.title = element_text(colour = "black", size = 11),
      legend.text = element_text(colour = "black", size = 9),
      plot.margin = unit(c(0, 0, 0, 0), "pt")
    )
  return(gg)
}
