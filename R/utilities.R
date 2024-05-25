## color palette
#' @importFrom grDevices colorRampPalette
#' @keywords internal
.set_default_cols <- function(n, type = 1) {
  col <- c(
    "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3",
    "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd",
    "#ccebc5", "#ffed6f"
  )
  col2 <- c(
    "#1f78b4", "#ffff33", "#c2a5cf", "#ff7f00", "#810f7c",
    "#a6cee3", "#006d2c", "#4d4d4d", "#8c510a", "#d73027",
    "#78c679", "#7f0000", "#41b6c4", "#e7298a", "#54278f"
  )
  col3 <- c(
    "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99",
    "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a",
    "#ffff99", "#b15928"
  )
  cols <- list(col, col2, col3)
  colorRampPalette(cols[[type]])(n)
}

## Suppress the print and cat
#' @keywords internal
suppressPrintAndCat <- function(x) {
  sink(tempfile(), type = "out")
  on.exit(sink())
  invisible(force(x))
}

#' @keywords internal
#' @importFrom utils head tail
rm_head_tail <- function(x, n) {
  x <- tail(head(x, -n), -n)
  return(x)
}

## Metabolite name break
#' @keywords internal
metabo_wrap <- function(metabo, max_len = 40, wrap = "\n") {
  metabo |>
    strsplit(split = "") |>
    lapply(function(x) {
      divisible <- length(x) %/% max_len
      lapply(seq_len(divisible + 1), function(i) {
        x[seq_len(max_len) + max_len * (i - 1)] |>
          na.omit() |>
          paste0(collapse = "")
      }) |>
        paste0(collapse = "\n")
    }) |>
    unlist()
}

## Descripton (funtional item) break
#' @keywords internal
desc_wrap <- function(desc, max_len = 40, wrap = "\n") {
  desc |>
    strsplit(split = "") |>
    lapply(function(x) {
      which_space <- grep(" ", x)
      if (length(which_space) > 0 && length(x) > max_len) {
        index <- cut(which_space, seq(1, length(x), max_len)) |>
          table() |>
          cumsum() |>
          (\(i) which_space[i])()
        x[index] <- wrap
      }
      paste0(x, collapse = "")
    }) |>
    unlist()
}

## Normalization
#' @keywords internal
norm_min_max <- function(x) {
  y <- na.omit(x)
  new.x <- (x - min(y)) / (max(y) - min(y))
  return(new.x)
}