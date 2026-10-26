#' S3 Methods and WGCNA Compatibility for rotatR
#'
#' @name compat
NULL

#' @export
print.rotatR <- function(x, ...) {
  mod_ids <- sort(setdiff(unique(x$labels), 0L))
  n_assigned <- sum(x$labels > 0L)
  p <- length(x$labels)
  cat(sprintf("rotatR fit (method = '%s', k = %d, beta = %s)\n",
              x$method, x$k, format(x$beta)))
  cat(sprintf("  Modules resolved : %d (min_size = %d, merge_cor = %.2f)\n",
              length(mod_ids), x$min_size, x$merge_cor))
  cat(sprintf("  Genes assigned   : %d / %d (%.1f%%)\n",
              n_assigned, p, 100 * n_assigned / max(1L, p)))
  if (length(mod_ids) > 0L) {
    sizes <- sort(table(x$labels[x$labels > 0L]), decreasing = TRUE)
    top_str <- paste(sprintf("M%s:%d", names(sizes), as.integer(sizes)), collapse = ", ")
    cat(sprintf("  Module sizes     : %s\n", top_str))
  }
  invisible(x)
}

#' @export
summary.rotatR <- function(object, ...) {
  mod_ids <- sort(setdiff(unique(object$labels), 0L))
  rows <- lapply(mod_ids, function(m) {
    idx <- which(object$labels == m)
    best_col <- max.col(abs(object$loadings[idx, , drop = FALSE]), ties.method = "first")
    l_vec <- object$loadings[cbind(idx, best_col)]
    top_ord <- order(abs(l_vec), decreasing = TRUE)
    hub_genes <- object$gene_names[idx[top_ord[seq_len(min(3L, length(idx)))]]]
    data.frame(
      module = m,
      size = length(idx),
      mean_abs_loading = mean(abs(l_vec)),
      max_abs_loading = max(abs(l_vec)),
      top_hubs = paste(hub_genes, collapse = ", "),
      stringsAsFactors = FALSE
    )
  })
  df <- if (length(rows) > 0L) {
    do.call(rbind, rows)
  } else {
    data.frame(
      module = integer(0), size = integer(0),
      mean_abs_loading = numeric(0), max_abs_loading = numeric(0),
      top_hubs = character(0), stringsAsFactors = FALSE
    )
  }
  structure(
    list(
      method = object$method,
      k = object$k,
      beta = object$beta,
      n_genes = length(object$labels),
      n_assigned = sum(object$labels > 0L),
      modules = df
    ),
    class = "summary.rotatR"
  )
}

#' @export
print.summary.rotatR <- function(x, ...) {
  cat(sprintf("Summary of rotatR (method = '%s', k = %d, beta = %s)\n",
              x$method, x$k, format(x$beta)))
  cat(sprintf("Assigned %d / %d genes across %d modules:\n\n",
              x$n_assigned, x$n_genes, nrow(x$modules)))
  if (nrow(x$modules) > 0L) {
    print(x$modules, row.names = FALSE)
  }
  invisible(x)
}

#' Project New Samples onto Discovered Modules
#'
#' @param object A fitted `rotatR` object.
#' @param newdata Numeric matrix or data.frame of new samples with the same gene
#'   columns as the training data.
#' @param ... Unused.
#' @return Numeric matrix of projected module scores.
#' @export
predict.rotatR <- function(object, newdata, ...) {
  if (missing(newdata) || is.null(newdata)) {
    return(object$eigengenes)
  }
  newdata <- as.matrix(newdata)
  if (!is.null(colnames(newdata)) && all(object$gene_names %in% colnames(newdata))) {
    newdata <- newdata[, object$gene_names, drop = FALSE]
  } else if (ncol(newdata) != length(object$gene_names)) {
    stop("`newdata` must have the same columns as the training matrix `x`.")
  }
  z_new <- standardize(newdata)
  mod_ids <- sort(setdiff(unique(object$labels), 0L))
  if (length(mod_ids) == 0L) {
    return(matrix(numeric(0), nrow = nrow(newdata), ncol = 0L))
  }
  scores <- matrix(0, nrow = nrow(newdata), ncol = length(mod_ids))
  colnames(scores) <- paste0("ME", mod_ids)
  rownames(scores) <- rownames(newdata)
  for (j in seq_along(mod_ids)) {
    m <- mod_ids[j]
    idx <- which(object$labels == m)
    best_col <- max.col(abs(object$loadings[idx, , drop = FALSE]), ties.method = "first")
    w <- object$loadings[cbind(idx, best_col)]
    w_norm <- sqrt(sum(w^2))
    if (w_norm > 0) w <- w / w_norm
    s <- as.numeric(z_new[, idx, drop = FALSE] %*% w)
    s_norm <- sqrt(sum(s^2))
    if (s_norm > 0) s <- s / s_norm
    scores[, j] <- s
  }
  scores
}

#' Extract Gene-Level Module Assignment and Hubness Table
#'
#' @param x A fitted `rotatR` object.
#' @param row.names Optional row names.
#' @param optional Logical.
#' @param ... Unused.
#' @return A `data.frame` with columns `gene`, `module`, `loading`, and
#'   `abs_loading`.
#' @export
as.data.frame.rotatR <- function(x, row.names = NULL, optional = FALSE, ...) {
  p <- length(x$labels)
  best_col <- max.col(abs(x$raw_loadings), ties.method = "first")
  loading <- x$raw_loadings[cbind(seq_len(p), best_col)]
  loading[x$labels == 0L] <- 0
  data.frame(
    gene = x$gene_names,
    module = x$labels,
    loading = loading,
    abs_loading = abs(loading),
    row.names = row.names,
    stringsAsFactors = FALSE
  )
}

#' Convert a rotatR Fit to a WGCNA-Compatible Network List
#'
#' @param x A fitted `rotatR` object.
#' @param use_color_names Logical; if `TRUE` and `WGCNA` is installed, converts
#'   integer labels to WGCNA color strings via `WGCNA::labels2colors()`.
#' @return A list with `colors`, `MEs`, `unmergedColors`, `MEsOK`, `dendrograms`,
#'   `blockGenes`, and `goodGenes`.
#' @export
as_wgcna <- function(x, use_color_names = FALSE) {
  if (!inherits(x, "rotatR")) {
    stop("`x` must be a `rotatR` object.")
  }
  colors <- x$labels
  if (use_color_names && requireNamespace("WGCNA", quietly = TRUE)) {
    colors <- WGCNA::labels2colors(x$labels)
    names(colors) <- x$gene_names
  }
  list(
    colors = colors,
    unmergedColors = colors,
    MEs = as.data.frame(x$eigengenes),
    goodSamples = rep(TRUE, nrow(x$eigengenes)),
    goodGenes = rep(TRUE, length(x$labels)),
    dendrograms = list(),
    TOMFiles = character(0),
    blockGenes = list(seq_along(x$labels)),
    blocks = rep(1L, length(x$labels)),
    MEsOK = TRUE
  )
}
