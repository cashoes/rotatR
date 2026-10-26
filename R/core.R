#' Core Varimax and Spectral Co-expression Module Discovery
#'
#' @name core
NULL

#' Standardize Columns of an Expression Matrix
#'
#' Centers and scales each gene column to mean 0 and unit variance, setting
#' zero-variance columns to 0.
#'
#' @param x Numeric matrix (samples in rows, genes in columns).
#' @return Standardized numeric matrix of the same dimensions.
#' @export
standardize <- function(x) {
  z <- scale(x)
  z[is.na(z)] <- 0
  z
}

#' Scale Columns to Unit Euclidean Norm
#'
#' @param x Numeric matrix (samples in rows, genes in columns).
#' @return Numeric matrix whose columns have inner products equal to sample
#'   correlations.
#' @export
unit_genes <- function(x) {
  standardize(x) / sqrt(nrow(x) - 1)
}

#' Fast Powered Absolute Value
#'
#' Computes `abs(r)^beta` using repeated squaring for integer exponents.
#'
#' @param r Numeric matrix or vector.
#' @param beta Numeric exponent.
#' @return Numeric matrix or vector `abs(r)^beta`.
#' @export
abs_pow <- function(r, beta) {
  a <- abs(r)
  if (beta < 1 || beta != round(beta)) {
    return(a^beta)
  }
  out <- NULL
  b <- as.integer(beta)
  while (b > 0L) {
    if (b %% 2L == 1L) {
      out <- if (is.null(out)) a else out * a
    }
    b <- b %/% 2L
    if (b > 0L) a <- a * a
  }
  out
}

#' Exact Powered-Correlation Kernel Matrix
#'
#' @param x Numeric matrix (samples in rows, genes in columns).
#' @param beta Numeric soft-threshold power.
#' @return Symmetric \eqn{p \times p} matrix \eqn{|R|^\beta}.
#' @export
kernel_gram <- function(x, beta) {
  r <- stats::cor(x)
  r[is.na(r)] <- 0
  abs_pow(r, beta)
}

#' Spectral Factor Decompositions
#'
#' @param x Numeric matrix (samples in rows, genes in columns).
#' @param rank Integer target factor rank.
#' @param oversample Integer oversampling parameter for randomized SVD.
#' @param n_iter Integer power iterations for randomized SVD.
#' @return A list with `v` (\eqn{p \times \text{rank}} eigenvectors) and `d`
#'   (singular values).
#' @export
factors_randomized <- function(x, rank = 20, oversample = 10, n_iter = 2) {
  rank <- min(rank, min(dim(x)) - 1)
  s <- rsvd::rsvd(standardize(x), k = rank, nu = 0, p = oversample, q = n_iter)
  list(v = s$v, d = s$d)
}

#' @rdname factors_randomized
#' @param gram Symmetric positive semi-definite \eqn{p \times p} matrix.
#' @export
factors_gram <- function(gram, rank = 20) {
  e <- RSpectra::eigs_sym(gram, rank, which = "LA")
  list(v = e$vectors, d = sqrt(pmax(e$values, 0)))
}

#' @rdname factors_randomized
#' @param phi Numeric \eqn{p \times D} feature matrix.
#' @export
factors_features <- function(phi, rank = 20) {
  rank <- min(rank, min(dim(phi)) - 1)
  s <- rsvd::rsvd(phi, k = rank, nv = 0, p = 10, q = 4)
  list(v = s$u, d = s$d)
}

#' Leverage-Weighted Nystrom Features of \eqn{|R|^\beta}
#'
#' Approximates the \eqn{p \times p} powered-correlation kernel \eqn{|R|^\beta}
#' using `m` landmark genes sampled proportional to their randomized SVD
#' leverage scores (or uniformly when `leverage_rank = 0`).
#'
#' @param x Numeric matrix (samples in rows, genes in columns).
#' @param beta Numeric soft-threshold power (default `2`).
#' @param m Integer number of landmarks (default `256`).
#' @param seed Integer random seed (default `0`).
#' @param leverage_rank Integer rank for leverage-score sampling (default `20`).
#' @return Numeric \eqn{p \times m'} feature matrix \eqn{\Phi} such that
#'   \eqn{\Phi \Phi^\top \approx |R|^\beta}.
#' @export
nystrom <- function(x, beta = 2, m = 256, seed = 0, leverage_rank = 20) {
  set.seed(seed)
  y <- unit_genes(x)
  m <- min(m, ncol(y))
  landmarks <- if (leverage_rank > 0) {
    leverage <- rowSums(factors_randomized(x, leverage_rank)$v^2)
    sample.int(ncol(y), m, prob = leverage)
  } else {
    sample.int(ncol(y), m)
  }
  kernel_cols <- abs(crossprod(y, y[, landmarks]))^beta
  e <- eigen(kernel_cols[landmarks, ], symmetric = TRUE)
  keep <- e$values > max(e$values) * 1e-8
  kernel_cols %*%
    sweep(e$vectors[, keep, drop = FALSE], 2, sqrt(e$values[keep]), "/")
}

#' @rdname nystrom
#' @export
nystrom_features <- nystrom

#' Kernel Spectral Factors (Exact or Nystrom)
#'
#' @param x Numeric matrix (samples in rows, genes in columns).
#' @param beta Numeric soft-threshold power.
#' @param rank Integer target factor rank.
#' @param m Optional integer Nystrom landmark count (`NULL` for exact kernel).
#' @param leverage_rank Integer rank for Nystrom leverage sampling.
#' @return A list with `v` and `d`.
#' @export
kernel_factors <- function(x, beta, rank, m = NULL, leverage_rank = rank) {
  if (is.null(m)) {
    factors_gram(kernel_gram(x, beta), rank)
  } else {
    factors_features(nystrom(x, beta, m, leverage_rank = leverage_rank), rank)
  }
}

#' Linear Varimax-Rotated PCA Loadings (vsp)
#'
#' Computes rank-`k` Varimax-rotated gene loadings on standardized expression
#' data using `vsp::vsp()` (Rohe & Zeng 2023) with degree normalization and
#' double centering disabled.
#'
#' @param x Numeric matrix (samples in rows, genes in columns).
#' @param k Integer number of components (\eqn{k \ge 2}).
#' @return Numeric \eqn{p \times k} loading matrix.
#' @export
vsp_loadings <- function(x, k) {
  stopifnot(k >= 2)
  fit <- vsp::vsp(
    standardize(x),
    rank = k,
    center = FALSE,
    degree_normalize = FALSE
  )
  unname(as.matrix(fit$Y))
}

#' Unrotated PCA Gene Loadings
#'
#' Computes the leading `k` right singular vectors of the standardized data,
#' sign-oriented by third-moment skewness to match `vsp`.
#'
#' @param x Numeric matrix (samples in rows, genes in columns).
#' @param k Integer number of components.
#' @return Numeric \eqn{p \times k} loading matrix.
#' @export
pca_loadings <- function(x, k) {
  z <- standardize(x)
  v <- base::svd(z, nu = 0, nv = min(k, min(dim(z))))$v
  skew <- apply(v, 2, \(col) sign(sum((col - mean(col))^3)))
  sweep(v, 2, ifelse(skew == 0, 1, skew), `*`)
}

#' Kernel Varimax-Rotated Loadings (kvsp / kvsp_nys)
#'
#' Applies raw Varimax rotation (`normalize = FALSE`) and skewness sign
#' alignment to the top `k` eigenvectors of \eqn{|R|^\beta} (either from a
#' precomputed factor list `fac` or directly from an expression matrix `x`).
#'
#' @param fac Either a factor list `list(v = ..., d = ...)` from
#'   \code{\link{kernel_factors}}, or a numeric expression matrix `x`.
#' @param k Integer number of components (\eqn{k \ge 2}).
#' @param beta Numeric soft-threshold power (used only when `fac` is a matrix).
#' @param m Optional integer Nystrom landmark count (used only when `fac` is a
#'   matrix).
#' @return Numeric \eqn{p \times k} loading matrix.
#' @export
kvsp_loadings <- function(fac, k, beta = 2, m = NULL) {
  if (is.matrix(fac) || is.data.frame(fac)) {
    fac <- kernel_factors(as.matrix(fac), beta = beta, rank = max(k, 2), m = m)
  }
  v <- fac$v[, seq_len(k), drop = FALSE]
  y <- v %*% stats::varimax(v, normalize = FALSE)$rotmat
  skew <- apply(y, 2, \(col) sign(sum((col - mean(col))^3)))
  sweep(y, 2, ifelse(skew == 0, 1, skew), `*`)
}

#' Threshold Loadings by Component-Wise Robust Z-Score
#'
#' Zeroes out loadings whose deviation from their component's median is less
#' than `z` median absolute deviations (MAD).
#'
#' @param loadings Numeric \eqn{p \times k} loading matrix.
#' @param z Numeric robust z-score cutoff.
#' @return Thresholded \eqn{p \times k} loading matrix.
#' @export
robust_cut_loadings <- function(loadings, z) {
  centre <- apply(loadings, 2, stats::median)
  spread <- apply(loadings, 2, stats::mad)
  dev <- abs(sweep(loadings, 2, centre))
  loadings[dev < rep(z * spread, each = nrow(loadings))] <- 0
  loadings
}

#' Threshold Loadings by False Discovery Rate on Robust Z-Scores
#'
#' Converts each loading's robust z-score into a two-sided normal p-value and
#' applies Benjamini-Hochberg FDR control at level `q` across all entries.
#'
#' @param loadings Numeric \eqn{p \times k} loading matrix.
#' @param q Numeric target FDR level (default `0.01`).
#' @return Thresholded \eqn{p \times k} loading matrix.
#' @export
fdr_cut_loadings <- function(loadings, q = 0.01) {
  centre <- apply(loadings, 2, stats::median)
  spread <- apply(loadings, 2, stats::mad)
  z <- sweep(sweep(loadings, 2, centre), 2, spread, `/`)
  keep <- matrix(
    stats::p.adjust(2 * stats::pnorm(-abs(z)), "BH") <= q,
    nrow(loadings)
  )
  keep[is.na(keep)] <- FALSE
  loadings[!keep] <- 0
  loadings
}

#' Assign Genes to Winner-Take-All Module Labels
#'
#' Assigns each gene with at least one non-zero loading to the component where
#' it has the largest absolute loading (`0L` if all loadings are zero).
#'
#' @param loadings Numeric \eqn{p \times k} loading matrix.
#' @return Integer vector of length \eqn{p} with `0L` for unassigned genes.
#' @export
loadings_to_labels <- function(loadings) {
  has_loading <- rowSums(abs(loadings) > 0) > 0
  ifelse(has_loading, max.col(abs(loadings), ties.method = "first"), 0L)
}

#' Data-Driven Gene Budget from Linear Varimax
#'
#' Estimates the number of module genes in `x` as the count of genes with at
#' least one non-zero loading after a strict robust z-score cut (`cut = 4`) on
#' rank-`k` linear Varimax loadings.
#'
#' @param x Numeric matrix (samples in rows, genes in columns).
#' @param k Integer number of components.
#' @param cut Numeric robust z-score cutoff (default `4`).
#' @return Integer count of assigned genes.
#' @export
gene_budget <- function(x, k, cut = 4) {
  y <- vsp_loadings(x, max(k, 2))
  sum(rowSums(robust_cut_loadings(y, cut) != 0) > 0)
}

#' First Left Singular Vector (Module Eigengene)
#'
#' @param x Numeric matrix (samples in rows, genes in columns).
#' @param keep Logical or integer index of genes in the module.
#' @return Numeric vector of length `nrow(x)`.
#' @export
module_eigengene <- function(x, keep) {
  z <- scale(x[, keep, drop = FALSE])
  z[is.na(z)] <- 0
  if (ncol(z) == 1) z[, 1] else base::svd(z, nu = 1, nv = 0)$u[, 1]
}

#' Merge Modules Whose Eigengenes Correlate Above a Threshold
#'
#' Iteratively clusters module eigengenes by average linkage on \eqn{1 - |r|}
#' and merges modules below height `1 - merge_cor` (matching WGCNA's
#' `mergeCloseModules` at `mergeCutHeight = 0.15`).
#'
#' @param x Numeric matrix (samples in rows, genes in columns).
#' @param labels Integer vector of module labels (`0` = unassigned).
#' @param merge_cor Numeric correlation threshold (default `0.85`).
#' @return Integer vector of merged module labels.
#' @export
merge_modules <- function(x, labels, merge_cor = 0.85) {
  repeat {
    ids <- sort(setdiff(unique(labels), 0L))
    if (length(ids) < 2) {
      break
    }
    eigengenes <- vapply(ids, \(id) module_eigengene(x, labels == id), numeric(nrow(x)))
    r <- abs(suppressWarnings(stats::cor(eigengenes)))
    r[is.na(r)] <- 0
    diag(r) <- 1
    groups <- stats::cutree(stats::hclust(stats::as.dist(1 - r), "average"), h = 1 - merge_cor)
    if (length(unique(groups)) == length(ids)) {
      break
    }
    labels[labels > 0] <- groups[match(labels[labels > 0], ids)]
  }
  labels
}

#' Drop Modules Smaller Than Minimum Size
#'
#' @param labels Integer vector of module labels (`0` = unassigned).
#' @param min_size Integer minimum module size (default `20`).
#' @return Integer vector of pruned module labels.
#' @export
drop_small_modules <- function(labels, min_size = 20) {
  sizes <- table(labels[labels > 0])
  labels[labels %in% as.integer(names(sizes)[sizes < min_size])] <- 0L
  labels
}

#' Renumber Module Labels Contiguously
#'
#' @param labels Integer vector of module labels (`0` = unassigned).
#' @return Integer vector with modules numbered `1, 2, ...` and `0` preserved.
#' @export
relabel <- function(labels) {
  match(labels, c(0L, sort(setdiff(unique(labels), 0L)))) - 1L
}

#' Resolve Module Labels via Merging and Minimum-Size Pruning
#'
#' @param x Numeric matrix (samples in rows, genes in columns).
#' @param labels Integer vector of initial module labels (`0` = unassigned).
#' @param merge_cor Numeric eigengene correlation merge threshold (default `0.85`).
#' @param min_size Integer minimum module size (default `20`).
#' @return Integer vector of resolved module labels.
#' @export
resolve_labels <- function(x, labels, merge_cor = 0.85, min_size = 20) {
  drop_small_modules(merge_modules(x, labels, merge_cor), min_size)
}

#' Resolve Loadings into Co-expression Modules
#'
#' Performs winner-take-all assignment (`loadings_to_labels`), optional mask
#' filtering (`keep` from stability selection), eigengene merging
#' (`merge_modules`), minimum module size pruning (`drop_small_modules`), and
#' contiguous relabeling.
#'
#' @param x Numeric matrix (samples in rows, genes in columns).
#' @param loadings Numeric \eqn{p \times k} loading matrix.
#' @param keep Optional logical vector of length \eqn{p} from stability selection.
#' @param merge Logical whether to merge close modules and drop small modules
#'   (default `TRUE`).
#' @param merge_cor Numeric eigengene correlation threshold for merging
#'   (default `0.85`).
#' @param min_size Integer minimum module size (default `20`).
#' @return Integer vector of module labels (`0` = unassigned, `1..M` = modules).
#' @export
resolve_modules <- function(x, loadings, keep = NULL, merge = TRUE,
                            merge_cor = 0.85, min_size = 20) {
  labels <- loadings_to_labels(loadings)
  if (!is.null(keep)) {
    labels[!keep] <- 0L
  }
  labels <- relabel(labels)
  if (merge) {
    labels <- relabel(resolve_labels(x, labels, merge_cor, min_size))
  }
  labels
}

#' Default Robust Z-Score Grids
#'
#' @export
robust_z_grid <- c(2, 3, 4, 5, 6, 8, 10)

#' @rdname robust_z_grid
#' @export
kernel_robust_z_grid <- c(4, 6, 8, 10, 15, 20, 30, 50, 100)

#' Fit a rotatR Co-expression Module Model
#'
#' Discovers co-expression modules via truncated SVD or powered-correlation
#' kernel spectral decomposition followed by Varimax rotation, loading
#' thresholding (analytical FDR control for linear `vsp`/`pca`, or subsampling
#' stability selection for `kvsp`/`kvsp_nys`), and module resolution.
#'
#' @param x Numeric matrix or data.frame (samples in rows, genes in columns).
#' @param method Character method: `"vsp"` (default linear Varimax PCA),
#'   `"kvsp_nys"` (Nystrom kernel Varimax PCA), `"kvsp"` (exact kernel Varimax
#'   PCA), or `"pca"` (unrotated PCA).
#' @param k Integer number of components, or `"eigcv"` (default) to estimate
#'   via 2-fold cross-projected eigenvalues capped at `k_max`.
#' @param k_max Integer maximum rank for `eigcv` (default `40L`).
#' @param alpha Numeric significance level for `eigcv` (default `0.01`).
#' @param beta Numeric soft-threshold power for kernel methods (default `2`).
#' @param m Integer Nystrom landmark count for `"kvsp_nys"` (default `256L`).
#' @param fdr Numeric BH-FDR cutoff for `"vsp"` and `"pca"` when `z` is `NULL`
#'   (default `0.01`).
#' @param z Optional numeric robust z-score cutoff. If `NULL`, `"vsp"` and
#'   `"pca"` use `fdr_cut_loadings(loadings, fdr)` while `"kvsp"` and
#'   `"kvsp_nys"` select `z` via stability selection over `grid`.
#' @param grid Numeric vector of candidate robust z-score thresholds for kernel
#'   stability selection (default `sort(kernel_robust_z_grid, decreasing = TRUE)`).
#' @param n_sub Integer number of subsamples for stability selection (default `5L`).
#' @param frac Numeric subsample fraction for stability selection (default `0.63`).
#' @param merge_cor Numeric eigengene correlation threshold for merging modules
#'   (default `0.85`).
#' @param min_size Integer minimum module size (default `20L`).
#' @param seed Integer random seed (default `0L`).
#' @return An S3 object of class `"rotatR"` containing:
#'   \itemize{
#'     \item `labels`: Integer module assignments (`0` = unassigned, `1..M` = modules)
#'     \item `eigengenes`: Numeric `n x M` matrix of module eigengenes
#'     \item `loadings`: Thresholded `p x k` loading matrix
#'     \item `raw_loadings`: Unthresholded `p x k` loading matrix
#'     \item `k`: Selected number of components
#'     \item `beta`: Soft-threshold power used (`1` for linear, `beta` for kernel)
#'     \item `param`: Selected sparsity parameter (`fdr` or `z`)
#'     \item `method`: Method string
#'     \item `gene_names`: Column names of `x`
#'   }
#' @export
rotatR <- function(
  x,
  method = c("vsp", "kvsp_nys", "kvsp", "pca"),
  k = "eigcv",
  k_max = 40L,
  alpha = 0.01,
  beta = 2,
  m = 256L,
  fdr = 0.01,
  z = NULL,
  grid = sort(kernel_robust_z_grid, decreasing = TRUE),
  n_sub = 5L,
  frac = 0.63,
  merge_cor = 0.85,
  min_size = 20L,
  seed = 0L
) {
  method <- match.arg(method)
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  n <- nrow(x)
  p <- ncol(x)
  gene_names <- colnames(x) %||% paste0("G", seq_len(p))
  colnames(x) <- gene_names

  beta_used <- if (method %in% c("vsp", "pca")) 1 else beta

  k_linear <- if (identical(k, "eigcv") || is.null(k)) {
    eigcv_k(x, cap = k_max, seed = seed, alpha = alpha)
  } else {
    as.integer(k)
  }

  use_stability <- is.null(z) && method %in% c("kvsp", "kvsp_nys")
  k_used <- if (use_stability) {
    min(k_linear, max(2L, floor(0.5 * frac * n)))
  } else {
    max(2L, k_linear)
  }

  compute_comp <- function(dat) {
    switch(
      method,
      vsp = vsp_loadings(dat, max(k_used, 2L)),
      pca = pca_loadings(dat, k_used),
      kvsp = kvsp_loadings(kernel_factors(dat, beta_used, max(k_used, 2L), NULL), max(k_used, 2L)),
      kvsp_nys = kvsp_loadings(kernel_factors(dat, beta_used, max(k_used, 2L), m, max(k_used, 2L)), max(k_used, 2L))
    )
  }

  raw_loadings <- compute_comp(x)
  keep <- NULL

  if (!is.null(z)) {
    loadings <- robust_cut_loadings(raw_loadings, z)
    param <- z
  } else if (method %in% c("vsp", "pca")) {
    loadings <- fdr_cut_loadings(raw_loadings, fdr)
    param <- fdr
  } else {
    budget <- gene_budget(x, k_linear, cut = 4)
    default_cap <- sqrt(0.05 * p * (2 * 0.8 - 1) * p)
    q_max <- if (is.na(budget) || budget == 0) {
      default_cap
    } else if (n < 100L) {
      max(default_cap, 1.45 * budget)
    } else {
      1.45 * budget
    }

    set.seed(seed)
    sub_comps <- lapply(seq_len(n_sub), \(i) {
      rows <- sample.int(n, round(frac * n))
      compute_comp(x[rows, , drop = FALSE])
    })

    eval_z <- function(val) {
      sel <- do.call(cbind, lapply(sub_comps, \(comp) {
        rowSums(abs(robust_cut_loadings(comp, val)) > 0) > 0
      }))
      list(q = mean(colSums(sel)), prob = rowMeans(sel))
    }

    max_prob <- numeric(p)
    last_ok <- NULL
    first_bad <- NULL
    for (val in grid) {
      s <- eval_z(val)
      if (s$q > q_max) {
        first_bad <- val
        break
      }
      max_prob <- pmax(max_prob, s$prob)
      last_ok <- val
    }
    if (!is.null(first_bad) && !is.null(last_ok)) {
      lo <- log10(last_ok)
      hi <- log10(first_bad)
      for (iter in 1:4) {
        mid <- 10^((lo + hi) / 2)
        s <- eval_z(mid)
        if (s$q > q_max) {
          hi <- log10(mid)
        } else {
          max_prob <- pmax(max_prob, s$prob)
          last_ok <- mid
          lo <- log10(mid)
        }
      }
    }
    if (is.null(last_ok)) {
      loadings <- matrix(0, nrow = p, ncol = ncol(raw_loadings))
      keep <- rep(FALSE, p)
      param <- NA_real_
    } else {
      loadings <- robust_cut_loadings(raw_loadings, last_ok)
      keep <- max_prob >= 0.8
      param <- last_ok
    }
  }

  labels <- resolve_modules(
    x, loadings, keep = keep, merge = TRUE,
    merge_cor = merge_cor, min_size = min_size
  )
  names(labels) <- gene_names
  rownames(loadings) <- gene_names
  rownames(raw_loadings) <- gene_names

  mod_ids <- sort(setdiff(unique(labels), 0L))
  eigengenes <- if (length(mod_ids) > 0L) {
    me <- vapply(mod_ids, \(id) module_eigengene(x, labels == id), numeric(n))
    me <- matrix(me, nrow = n, ncol = length(mod_ids))
    colnames(me) <- paste0("ME", mod_ids)
    rownames(me) <- rownames(x)
    me
  } else {
    matrix(numeric(0), nrow = n, ncol = 0L)
  }

  structure(
    list(
      labels = labels,
      eigengenes = eigengenes,
      loadings = loadings,
      raw_loadings = raw_loadings,
      k = k_used,
      k_linear = k_linear,
      beta = beta_used,
      param = param,
      method = method,
      min_size = min_size,
      merge_cor = merge_cor,
      gene_names = gene_names
    ),
    class = "rotatR"
  )
}

#' @rdname rotatR
#' @export
spec_wgcna <- rotatR

`%||%` <- function(a, b) if (is.null(a)) b else a
