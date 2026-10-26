#' 9-Panel Diagnostic Model Card for rotatR
#'
#' @name diagnostics
#' @import ggplot2
#' @import patchwork
NULL

diag_eigcv <- function(x, seed = 0, cap = 40, alpha = 0.01, n_splits = 2L, show = 60) {
  prof <- eigcv_profile(
    x, cap = max(show, cap), seed = seed,
    alpha = alpha, n_splits = n_splits
  )
  cap_val <- if (is.finite(cap)) cap else Inf
  list(
    k = as.integer(min(prof$k, cap_val)),
    k_uncapped = prof$k,
    cap = cap_val,
    table = tibble::tibble(
      component = seq_along(prof$lambda),
      eigenvalue = prof$lambda,
      null_lo = 2 - prof$threshold,
      null_hi = prof$threshold
    )
  )
}

diag_survival <- function(x, k, beta, n_perm = 20, seed = 0) {
  null <- vapply(
    seq_len(n_perm),
    \(i) {
      set.seed(seed + 104729L * i)
      kernel_survival_scores(apply(x, 2, sample), k, beta)
    },
    numeric(max(k, 2))
  )
  list(
    beta = beta,
    score = kernel_survival_scores(x, k, beta),
    null = as.numeric(null),
    null_bar = stats::quantile(null, 0.95, na.rm = TRUE, names = FALSE)
  )
}

diag_cut <- function(loadings, q = 0.01) {
  centre <- apply(loadings, 2, stats::median)
  spread <- apply(loadings, 2, stats::mad)
  z <- sweep(sweep(loadings, 2, centre), 2, spread, `/`)
  keep <- matrix(
    stats::p.adjust(2 * stats::pnorm(-abs(z)), "BH") <= q,
    nrow(loadings)
  )
  keep[is.na(keep)] <- FALSE
  list(z = z, keep = keep, q = q, threshold = min(abs(z[keep])))
}

diag_eigengenes <- function(x, labels) {
  ids <- sort(setdiff(unique(labels), 0L))
  e <- vapply(
    ids,
    \(id) {
      z <- scale(x[, labels == id, drop = FALSE])
      z[is.na(z)] <- 0
      if (ncol(z) == 1) z[, 1] else base::svd(z, nu = 1, nv = 0)$u[, 1]
    },
    numeric(nrow(x))
  )
  colnames(e) <- ids
  e
}

diag_heatmap <- function(x, loadings, raw_labels, merged, n_unassigned = 120,
                         max_per_module = Inf, seed = 0) {
  labels <- raw_labels
  own <- loadings[cbind(seq_along(raw_labels), pmax(raw_labels, 1L))]
  set.seed(seed)
  unassigned <- which(labels == 0)
  shown <- sample(unassigned, min(n_unassigned, length(unassigned)))
  rank_in_module <- stats::ave(-abs(own), labels, FUN = rank)
  assigned <- which(labels > 0 & rank_in_module <= max_per_module)
  genes <- c(assigned[order(labels[assigned], -own[assigned])], shown)
  r <- stats::cor(x[, genes])
  r[is.na(r)] <- 0
  list(
    r = unname(r),
    genes = tibble::tibble(
      pos = seq_along(genes),
      module = labels[genes],
      merged = merged[genes]
    ),
    n_unassigned = length(unassigned),
    n_assigned = sum(labels > 0),
    eigengene_cor = abs(stats::cor(diag_eigengenes(x, labels)))
  )
}

diag_stability <- function(x, method = "kvsp", k, beta = 2, q_max, seed = 0,
                           n_sub = 5L, frac = 0.63, threshold = 0.8, m = 256L) {
  prep_fn <- switch(
    method,
    kvsp = \(z) kvsp_loadings(factors_gram(kernel_gram(z, beta), max(k, 2)), max(k, 2)),
    kvsp_nys = \(z) {
      phi <- nystrom(z, beta, min(m, ncol(z)), leverage_rank = max(k, 2))
      kvsp_loadings(factors_features(phi, max(k, 2)), max(k, 2))
    },
    stop("Unknown method in diag_stability: ", method)
  )
  grid <- sort(kernel_robust_z_grid, decreasing = TRUE)
  n <- nrow(x)
  p <- ncol(x)
  set.seed(seed)
  subsamples <- purrr::map(seq_len(n_sub), \(i) {
    sample.int(n, round(frac * n))
  })
  states <- purrr::map(subsamples, \(rows) prep_fn(x[rows, , drop = FALSE]))

  selection <- function(value) {
    sel <- purrr::map(states, \(st) rowSums(abs(robust_cut_loadings(st, value)) > 0) > 0) |>
      do.call(what = cbind)
    list(q = mean(colSums(sel)), prob = rowMeans(sel))
  }

  max_prob <- numeric(p)
  last_ok <- NULL
  first_bad <- NULL
  for (value in grid) {
    sel <- selection(value)
    if (sel$q > q_max) {
      first_bad <- value
      break
    }
    max_prob <- pmax(max_prob, sel$prob)
    last_ok <- value
  }
  if (!is.null(first_bad) && !is.null(last_ok)) {
    lo <- log10(last_ok)
    hi <- log10(first_bad)
    for (iter in 1:4) {
      mid <- 10^((lo + hi) / 2)
      sel <- selection(mid)
      if (sel$q > q_max) {
        hi <- log10(mid)
      } else {
        max_prob <- pmax(max_prob, sel$prob)
        last_ok <- mid
        lo <- log10(mid)
      }
    }
  }

  full <- prep_fn(x)
  selected <- vapply(
    states,
    \(st) rowSums(abs(robust_cut_loadings(st, last_ok)) > 0) > 0,
    logical(p)
  )
  z_full <- abs(sweep(sweep(full, 2, apply(full, 2, stats::median)), 2, apply(full, 2, stats::mad), `/`))
  keep <- max_prob >= threshold
  lbls <- loadings_to_labels(robust_cut_loadings(full, last_ok))
  lbls[!keep] <- 0L
  labels <- relabel(lbls)

  list(
    method = method,
    cut = last_ok,
    labels = labels,
    loadings = full,
    genes = tibble::tibble(
      z = apply(z_full, 1, max),
      frequency = rowMeans(selected),
      label = labels
    )
  )
}

diag_landmarks <- function(x, m, seed = 0, leverage_rank = 20) {
  set.seed(seed)
  leverage <- rowSums(factors_randomized(x, leverage_rank)$v^2)
  sample.int(ncol(x), min(m, ncol(x)), prob = leverage)
}

match_loadings <- function(a, b) {
  r <- stats::cor(a, b)
  perm <- as.integer(clue::solve_LSAP(abs(r), maximum = TRUE))
  signs <- sign(r[cbind(seq_along(perm), perm)])
  sweep(b[, perm, drop = FALSE], 2, signs, `*`)
}

diag_nystrom <- function(x, k, beta, labels, m = 256,
                         m_grid = c(32, 64, 128, 256, 512, 1024)) {
  m_grid <- m_grid[m_grid >= 2 * k]
  exact <- kvsp_loadings(factors_gram(kernel_gram(x, beta), k), k)
  sketch <- function(m_val) {
    phi <- nystrom(x, beta, min(m_val, ncol(x)), leverage_rank = k)
    match_loadings(exact, kvsp_loadings(factors_features(phi, k), k))
  }
  nys <- sketch(m)
  landmarks <- diag_landmarks(x, m, leverage_rank = k)
  list(
    m = m,
    loadings = tibble::tibble(
      component = rep(seq_len(k), each = nrow(exact)),
      exact = as.numeric(exact),
      nystrom = as.numeric(nys)
    ),
    by_m = purrr::map_dfr(m_grid, \(m_val) {
      tibble::tibble(
        m = m_val,
        component = seq_len(k),
        r = diag(stats::cor(exact, sketch(m_val)))
      )
    }),
    coverage = tibble::tibble(module = labels) |>
      dplyr::mutate(landmark = seq_along(labels) %in% landmarks) |>
      dplyr::summarise(
        size = dplyr::n(),
        landmarks = sum(landmark),
        .by = module
      ) |>
      dplyr::mutate(uniform = m * size / length(labels)) |>
      dplyr::arrange(module)
  )
}

#' Compute Diagnostic Model Card Data
#'
#' Computes all 9 diagnostic panels for an expression dataset `sim` (either a
#' `list(x = ..., labels = ...)` or a numeric matrix `x`).
#'
#' @param sim A list with `x` (samples x genes matrix) and optional `labels`,
#'   or a numeric matrix `x` directly.
#' @param seed Integer random seed (default `0`).
#' @param fdr Numeric BH-FDR cutoff for linear `vsp` (default `0.01`).
#' @param merge_cor Numeric eigengene merge correlation cutoff (default `0.85`).
#' @param min_size Integer minimum module size (default `20`).
#' @param kernel_beta Numeric soft-threshold power for `kvsp` (default `2`).
#' @param survival_beta Numeric soft-threshold power for survival profile (default `2`).
#' @param m Integer Nystrom landmark count (default `256`).
#' @param pa_cap Integer maximum rank cap for `eigcv` (default `40`).
#' @param max_per_module Numeric maximum genes per module shown in the heatmap (default `Inf`).
#' @return A list of diagnostic structures consumed by \code{\link{plot_model_card}}.
#' @export
diag_compute <- function(sim, seed = 0, fdr = 0.01, merge_cor = 0.85,
                         min_size = 20, kernel_beta = 2,
                         survival_beta = 2, m = 256L,
                         pa_cap = 40, max_per_module = Inf) {
  if (is.matrix(sim) || is.data.frame(sim)) {
    sim <- list(x = as.matrix(sim), labels = NULL)
  }
  x <- sim$x
  pa <- diag_eigcv(x, seed = seed, cap = pa_cap)
  k <- pa$k
  y <- vsp_loadings(x, max(k, 2))
  cut <- diag_cut(y, fdr)
  raw_labels <- loadings_to_labels(fdr_cut_loadings(y, fdr))
  vsp_labels <- raw_labels
  merged <- relabel(resolve_labels(x, relabel(vsp_labels), merge_cor, min_size))
  budget <- gene_budget(x, k, cut = 4)
  default_cap <- sqrt(0.05 * ncol(x) * (2 * 0.8 - 1) * ncol(x))
  q_max <- if (is.na(budget) || budget == 0) {
    default_cap
  } else if (nrow(x) < 100L) {
    max(default_cap, 1.45 * budget)
  } else {
    1.45 * budget
  }
  stability <- diag_stability(x, "kvsp", k, kernel_beta, q_max, seed, m = m)
  list(
    n = nrow(x),
    p = ncol(x),
    truth = sim$labels,
    pa = pa,
    survival = list(
      kvsp = diag_survival(x, k, kernel_beta, seed = seed),
      guided = diag_survival(x, k, survival_beta, seed = seed)
    ),
    cut = cut,
    vsp_labels = vsp_labels,
    merged = merged,
    merge_cor = merge_cor,
    heatmap = diag_heatmap(
      x, y, raw_labels, merged, max_per_module = max_per_module, seed = seed
    ),
    stability = stability,
    stab_threshold = 0.8,
    nystrom = diag_nystrom(x, max(k, 2), kernel_beta, stability$labels, m)
  )
}

diag_npg <- ggsci::pal_npg()(10)
diag_ink <- "grey20"
diag_null <- "grey82"
diag_kept <- diag_npg[4]
diag_dropped <- diag_npg[1]

mname <- function(id) sprintf("M%02d", as.integer(id))

module_palette <- function(ids, option = "viridis") {
  ids <- sort(setdiff(ids, 0L))
  colours <- if (length(ids) <= 9) {
    diag_npg[c(4, 1, 3, 5, 2, 7, 9, 6, 10)][seq_along(ids)]
  } else {
    viridisLite::viridis(length(ids), begin = 0.08, end = 0.92, option = option)
  }
  c("0" = "grey85", purrr::set_names(colours, ids))
}

small_labels <- function(n, axis = c("y", "x"), limit = 12) {
  if (n <= limit) {
    return(NULL)
  }
  text <- element_text(size = 4.5)
  switch(
    match.arg(axis),
    y = theme(axis.text.y = text),
    x = theme(axis.text.x = element_text(size = 4.5, angle = 90, hjust = 1, vjust = 0.5))
  )
}

panel_scree <- function(d, labels = c("counted", "not counted")) {
  pa <- d$pa$table |>
    dplyr::mutate(counted = component <= d$pa$k)
  half <- max(2L, round(0.2 * nrow(pa)))
  centres <- seq(half + 1L, nrow(pa) - half)
  label_y <- vapply(centres, \(i) max(pa$null_hi[(i - half):(i + half)]) * 1.15, 0)
  clearance <- vapply(seq_along(centres), \(j) {
    window <- (centres[j] - half):(centres[j] + half)
    min(abs(log(pa$eigenvalue[window]) - log(label_y[j])))
  }, 0)
  best <- which.max(clearance)
  band_at <- tibble::tibble(component = centres[best], y = label_y[best])
  ggplot(pa, aes(component)) +
    geom_ribbon(aes(ymin = null_lo, ymax = null_hi), fill = diag_null) +
    geom_line(aes(y = eigenvalue), colour = diag_ink, linewidth = 0.3) +
    geom_point(aes(y = eigenvalue, fill = counted), shape = 21, size = 1.7) +
    annotate(
      "text", x = band_at$component, y = band_at$y,
      label = "null, 1% level", size = 2.6, colour = "grey40", vjust = 0
    ) +
    scale_fill_manual(
      NULL,
      values = c(`TRUE` = diag_kept, `FALSE` = "white"),
      labels = c(`TRUE` = labels[1], `FALSE` = labels[2]),
      breaks = c("TRUE", "FALSE")
    ) +
    scale_x_continuous(
      "Component",
      breaks = \(limits) c(1, seq(if (limits[2] > 35) 10 else 5, limits[2], if (limits[2] > 35) 10 else 5))
    ) +
    scale_y_log10(
      "Out-of-sample eigenvalue",
      breaks = \(lim) {
        decades <- log10(lim[2]) - log10(lim[1])
        if (decades < 1.3) {
          candidates <- c(20, 30, 50, 70, 100, 200, 300, 500)
          brks <- candidates[candidates >= lim[1] & candidates <= lim[2]]
          if (length(brks) >= 2) return(brks)
        }
        10^seq(ceiling(log10(lim[1])), floor(log10(lim[2])))
      },
      labels = \(b) {
        vapply(b, \(val) {
          if (abs(log10(val) - round(log10(val))) < 1e-6) {
            paste0("10^", round(log10(val)))
          } else {
            as.character(val)
          }
        }, character(1)) |>
          parse(text = _)
      }
    ) +
    labs(title = sprintf("Cross-validated eigenvalues: %d components", d$pa$k))
}

panel_survival <- function(d, survival = c("guided", "kvsp"), threshold = 1.1,
                           compact = FALSE) {
  survival <- d$survival[[match.arg(survival)]]
  scores <- tibble::tibble(
    component = seq_along(survival$score),
    score = survival$score
  ) |>
    dplyr::filter(!is.na(score)) |>
    dplyr::mutate(
      survives = score > threshold,
      component = forcats::fct_reorder(factor(component), -score)
    )
  bars <- tibble::tibble(
    y = c(survival$null_bar, threshold),
    label = c("permuted data, 95%", sprintf("threshold %.1f", threshold))
  )
  ggplot(scores, aes(component, score)) +
    geom_hline(aes(yintercept = y, linetype = label), bars, colour = "grey40") +
    geom_segment(aes(xend = component, yend = 1), colour = diag_ink, linewidth = 0.3) +
    geom_point(
      aes(fill = survives), shape = 21, size = if (nrow(scores) > 12) 1.3 else 2.2,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = c(`TRUE` = diag_kept, `FALSE` = "white")) +
    scale_linetype_manual(NULL, values = c(2, 1)) +
    scale_y_log10() +
    guides(linetype = guide_legend(nrow = if (compact) 2 else 1)) +
    labs(
      x = "Rotated component",
      y = bquote("Kernel survival score, " * beta == .(survival$beta)),
      title = sprintf("Kernel survival: %d of %d", sum(scores$survives), nrow(scores))
    ) +
    small_labels(nrow(scores), "x")
}

panel_cut <- function(d, z_max = 12, binwidth = 0.25, show_cut = TRUE) {
  z <- abs(as.numeric(d$cut$z))
  thr <- d$cut$threshold

  steps_down <- seq(thr, 0, by = -binwidth)
  if (utils::tail(steps_down, 1) > 0) steps_down <- c(steps_down, 0)
  breaks_left <- sort(unique(steps_down))

  steps_up <- seq(thr, z_max, by = binwidth)
  if (utils::tail(steps_up, 1) < z_max) steps_up <- c(steps_up, z_max)
  breaks_right <- sort(unique(steps_up))

  z_left <- z[z < thr]
  h_left <- graphics::hist(z_left, breaks = breaks_left, include.lowest = TRUE, right = FALSE, plot = FALSE)
  df_left <- tibble::tibble(
    xmin = h_left$breaks[-length(h_left$breaks)],
    xmax = h_left$breaks[-1],
    count = h_left$counts,
    kept = FALSE
  )

  z_right <- pmin(z[z >= thr], z_max)
  h_right <- graphics::hist(z_right, breaks = breaks_right, include.lowest = TRUE, right = TRUE, plot = FALSE)
  df_right <- tibble::tibble(
    xmin = h_right$breaks[-length(h_right$breaks)],
    xmax = h_right$breaks[-1],
    count = h_right$counts,
    kept = TRUE
  )

  df_all <- dplyr::bind_rows(df_left, df_right)

  grid <- seq(binwidth / 2, z_max - binwidth / 2, binwidth)
  n_loadings <- length(z)
  null <- tibble::tibble(
    z = grid,
    count = n_loadings * (stats::pnorm(grid + binwidth / 2) - stats::pnorm(grid - binwidth / 2)) * 2
  ) |>
    dplyr::filter(count >= 0.5)

  genes <- sum(rowSums(d$cut$keep) > 0)
  y_max <- max(c(null$count, df_all$count))

  y_breaks <- if (y_max > 50000) {
    c(0, 10000, 50000, 100000)
  } else if (y_max > 15000) {
    c(0, 5000, 15000, 30000)
  } else {
    c(0, 1000, 2500, 5000)
  }

  format_k <- function(b) {
    ifelse(b == 0, "0", ifelse(b >= 1000, paste0(b / 1000, "k"), as.character(b)))
  }

  ggplot(df_all) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = count, fill = kept), colour = NA) +
    geom_line(aes(x = z, y = count, linetype = "normal null"), null, colour = diag_ink, linewidth = 0.4) +
    geom_vline(xintercept = thr, linewidth = 0.4) +
    annotate(
      "text", x = thr + 0.15, y = Inf, hjust = 0, vjust = 1.4, size = 2.8,
      label = if (show_cut) {
        sprintf("FDR %g%%: |z| \u2265 %.2f", 100 * d$cut$q, thr)
      } else {
        sprintf("FDR %g%%", 100 * d$cut$q)
      }
    ) +
    scale_fill_manual(
      NULL,
      values = c(`TRUE` = diag_kept, `FALSE` = "grey75"),
      labels = c(`TRUE` = "assigned", `FALSE` = "not assigned"),
      breaks = c("TRUE", "FALSE")
    ) +
    scale_linetype_manual(NULL, values = 2) +
    scale_x_continuous(
      breaks = seq(0, z_max, 2),
      labels = \(b) ifelse(b == z_max, paste0("\u2265 ", b), b),
      expand = expansion(mult = c(0, 0.01))
    ) +
    scale_y_continuous(
      breaks = y_breaks,
      labels = format_k,
      expand = expansion(mult = c(0, 0.05))
    ) +
    coord_trans(y = "sqrt") +
    labs(x = "Robust z-score of loading, |z|", y = "Loadings (sqrt scale)", title = sprintf(
      "%s genes assigned of %s",
      format(genes, big.mark = ","), format(d$p, big.mark = ",")
    ))
}

panel_heatmap <- function(d, strips = TRUE, legend = "bottom", fixed = TRUE,
                          scale_name = "Correlation") {
  h <- d$heatmap
  g <- nrow(h$r)
  genes <- h$genes
  colours <- module_palette(c(d$vsp_labels, d$merged))
  cells <- tibble::tibble(
    row = rep(seq_len(g), times = g),
    col = rep(seq_len(g), each = g),
    r = as.numeric(h$r)
  )
  runs <- function(column) {
    genes |>
      dplyr::summarise(
        from = min(pos) - 0.5,
        to = max(pos) + 0.5,
        .by = c(module, dplyr::all_of(column))
      ) |>
      dplyr::mutate(fill = colours[as.character(.data[[column]])])
  }
  fitted <- runs("module") |>
    dplyr::mutate(
      label = ifelse(module == 0, "none", mname(module)),
      mid = (from + to) / 2
    )
  if (nrow(fitted) > 12) {
    fitted$label[fitted$module > 0 & seq_len(nrow(fitted)) %% 2 == 0] <- ""
  }
  after <- runs("merged")
  w <- g * 0.025
  edges <- fitted$to[-nrow(fitted)]
  shown <- sum(genes$module > 0)
  n_assigned <- sum(d$vsp_labels > 0)
  ggplot(cells) +
    geom_raster(aes(col, row, fill = r)) +
    geom_hline(yintercept = edges, colour = "white", linewidth = 0.25) +
    geom_vline(xintercept = edges, colour = "white", linewidth = 0.25) +
    (if (strips) {
      annotate(
        "rect", xmin = -2.2 * w, xmax = -1.2 * w, ymin = fitted$from, ymax = fitted$to,
        fill = fitted$fill
      )
    }) +
    (if (strips) {
      annotate(
        "rect", xmin = -1.1 * w, xmax = -0.1 * w, ymin = after$from, ymax = after$to,
        fill = after$fill
      )
    }) +
    scale_fill_gradient2(
      scale_name, low = diag_npg[4], mid = "white", high = diag_npg[8],
      limits = c(-1, 1), breaks = c(-1, 0, 1)
    ) +
    scale_y_reverse(
      NULL, breaks = fitted$mid, labels = fitted$label, expand = c(0, 0)
    ) +
    scale_x_continuous(
      NULL,
      breaks = if (strips) c(-1.7 * w, -0.6 * w) else NULL,
      labels = if (strips) c("fitted", "merged") else NULL,
      expand = c(0, 0)
    ) +
    (if (fixed) coord_fixed()) +
    labs(title = if (shown < n_assigned) {
      sprintf(
        "Correlation: the %d strongest of %s module genes, %d of %s unassigned",
        shown, format(n_assigned, big.mark = ","), sum(genes$module == 0),
        format(h$n_unassigned, big.mark = ",")
      )
    } else {
      sprintf(
        "Correlation of %d assigned genes and %d of %s unassigned",
        shown, sum(genes$module == 0), format(h$n_unassigned, big.mark = ",")
      )
    }) +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      legend.position = legend,
      legend.key.width = unit(1.6, "lines")
    ) +
    small_labels(nrow(fitted))
}

panel_eigengenes <- function(d) {
  r <- d$heatmap$eigengene_cor
  ids <- colnames(r)
  eig <- tibble::tibble(
    a = rep(ids, times = length(ids)),
    b = rep(ids, each = length(ids)),
    r = as.numeric(r)
  ) |>
    dplyr::mutate(
      merges = r > d$merge_cor & a != b,
      a = factor(mname(a), mname(ids)),
      b = factor(mname(b), rev(mname(ids)))
    )
  few <- length(ids) <= 10
  ggplot(eig, aes(a, b)) +
    geom_tile(aes(fill = r), colour = "white", linewidth = if (few) 0.4 else 0.1) +
    geom_tile(data = \(e) dplyr::filter(e, merges), fill = NA, colour = "black", linewidth = 0.6) +
    (if (few) geom_text(aes(label = sprintf("%.2f", r), colour = r > 0.6), size = 2.3)) +
    scale_fill_gradient(low = "white", high = diag_npg[8], limits = c(0, 1), guide = "none") +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey20"), guide = "none") +
    scale_x_discrete(NULL, expand = c(0, 0), guide = guide_axis(check.overlap = TRUE)) +
    scale_y_discrete(NULL, expand = c(0, 0), guide = guide_axis(check.overlap = TRUE)) +
    coord_fixed() +
    labs(title = sprintf("Eigengene |r| (merge above %.2f)", d$merge_cor)) +
    theme(axis.line = element_blank(), axis.ticks = element_blank()) +
    small_labels(length(ids), "x") +
    small_labels(length(ids), "y")
}

panel_sizes <- function(d, both = TRUE, option = "viridis") {
  ids <- sort(setdiff(unique(d$vsp_labels), 0L))
  colours <- module_palette(c(d$vsp_labels, d$merged), option)
  sizes <- dplyr::bind_rows(
    fitted = tibble::tibble(module = d$vsp_labels),
    merged = tibble::tibble(module = d$merged),
    .id = "stage"
  ) |>
    dplyr::filter(module > 0) |>
    dplyr::count(stage, module) |>
    dplyr::mutate(name = factor(mname(module), rev(mname(ids)))) |>
    dplyr::filter(both | stage == "merged")
  ggplot(sizes, aes(n, name)) +
    geom_col(aes(fill = name, alpha = stage), position = position_dodge(0.8), width = 0.75) +
    (if (length(ids) <= 12) {
      geom_text(
        aes(label = n, group = stage), position = position_dodge(0.8),
        hjust = -0.15, size = 2.3
      )
    }) +
    scale_fill_manual(
      values = purrr::set_names(colours[as.character(ids)], mname(ids)), guide = "none"
    ) +
    scale_alpha_manual(
      NULL, values = c(fitted = 1, merged = if (both) 0.45 else 1),
      guide = if (both) "legend" else "none"
    ) +
    scale_x_continuous(NULL, expand = expansion(mult = c(0, 0.18))) +
    scale_y_discrete(guide = guide_axis(check.overlap = TRUE)) +
    labs(y = NULL, title = "Genes per module") +
    small_labels(length(ids))
}

panel_frequency <- function(d, compact = FALSE, counts = TRUE, quadrants = FALSE) {
  s <- d$stability
  set.seed(1)
  genes <- s$genes |>
    dplyr::mutate(
      status = dplyr::case_when(
        label > 0 ~ "assigned",
        z >= s$cut ~ "past the cut, unstable",
        frequency >= d$stab_threshold ~ "stable, short of the cut",
        .default = "not assigned"
      ),
      status = factor(status, c(
        "assigned", "past the cut, unstable", "stable, short of the cut", "not assigned"
      )),
      jittered = frequency + stats::runif(dplyr::n(), -0.03, 0.03)
    ) |>
    dplyr::arrange(dplyr::desc(status))
  n_status <- dplyr::count(genes, status)
  corners <- tibble::tibble(
    status = factor(levels(genes$status), levels(genes$status)),
    label = c("assigned", "unstable", "short of cut", "not assigned"),
    x = s$cut * c(1.15, 1.15, 1 / 1.15, 1 / 1.15),
    y = c(1.13, -0.13, 1.13, -0.13),
    hjust = c(0, 0, 1, 1)
  ) |>
    dplyr::filter(status %in% n_status$status)
  ggplot(genes, aes(z, jittered)) +
    geom_hline(yintercept = d$stab_threshold - 0.05, linetype = 2, colour = "grey40") +
    geom_vline(xintercept = s$cut, linetype = 2, colour = "grey40") +
    geom_point(aes(colour = status), size = 0.5, alpha = 0.6) +
    scale_colour_manual(
      NULL,
      values = purrr::set_names(c(diag_kept, diag_dropped, diag_npg[5], "grey65"), levels(genes$status)),
      limits = as.character(n_status$status),
      labels = if (counts) {
        sprintf("%s (%s)", n_status$status, format(n_status$n, big.mark = ",", trim = TRUE))
      } else {
        sub("past the cut, |stable, ", "", n_status$status)
      }
    ) +
    scale_x_log10() +
    (if (quadrants) {
      geom_text(
        aes(x, y, label = label, colour = status, hjust = hjust), corners,
        size = 2.3, show.legend = FALSE
      )
    }) +
    scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    (if (quadrants) theme(legend.position = "none")) +
    labs(
      x = "Robust z-score of loading, all samples",
      y = "Selection frequency, 5 subsamples"
    ) +
    guides(colour = guide_legend(
      nrow = if (compact) 2 else 1, override.aes = list(size = 2, alpha = 1)
    )) +
    labs(title = sprintf(
      "%s stability: |z| \u2265 %.3g in %g%% of subsamples",
      s$method, s$cut, 100 * d$stab_threshold
    ))
}

component_colours <- function(d, option = "viridis") {
  comps <- sort(unique(d$nystrom$loadings$component))
  purrr::set_names(module_palette(comps, option)[-1], comps)
}

panel_nystrom_by_m <- function(d, option = "viridis") {
  ny <- d$nystrom
  ggplot(ny$by_m, aes(m, abs(r))) +
    geom_vline(xintercept = ny$m, linetype = 2, colour = "grey40") +
    geom_line(aes(colour = factor(component)), linewidth = 0.4) +
    geom_point(aes(colour = factor(component)), size = 1) +
    scale_colour_manual(values = component_colours(d, option), guide = "none") +
    scale_x_log10("Landmarks", breaks = unique(ny$by_m$m)) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(y = "Correlation with exact loadings", title = "Nystr\u00f6m agreement by landmarks")
}

panel_coverage <- function(d, option = "viridis") {
  cover <- d$nystrom$coverage |>
    dplyr::mutate(
      name = ifelse(module == 0, "none", sprintf("K%02d", module)),
      name = factor(name, rev(name))
    )
  fills <- purrr::set_names(
    module_palette(cover$module, option)[as.character(cover$module)], cover$name
  )
  ggplot(cover, aes(landmarks, name)) +
    geom_col(aes(fill = name), width = 0.7, show.legend = FALSE) +
    geom_point(
      aes(x = uniform, shape = "expected if drawn uniformly"),
      size = if (nrow(cover) > 12) 0.7 else 1.6, stroke = if (nrow(cover) > 12) 0.4 else 0.6
    ) +
    guides(shape = guide_legend(override.aes = list(size = 1.6, stroke = 0.6))) +
    (if (nrow(cover) <= 12) geom_text(aes(label = landmarks), hjust = -0.3, size = 2.3)) +
    scale_fill_manual(values = fills) +
    scale_shape_manual(NULL, values = 5) +
    scale_x_continuous("Landmark genes", expand = expansion(mult = c(0, 0.15))) +
    scale_y_discrete(guide = guide_axis(check.overlap = TRUE)) +
    labs(y = NULL, title = "Landmarks per kvsp module") +
    small_labels(nrow(cover))
}

card_header <- function(d, name) {
  stats <- tibble::tibble(
    label = c(
      "samples", "genes", "components", "modules", "genes in modules",
      "kvsp modules", "kvsp genes"
    ),
    value = c(
      d$n, d$p, d$pa$k, length(setdiff(unique(d$merged), 0L)), sum(d$merged > 0),
      length(setdiff(unique(d$stability$labels), 0L)), sum(d$stability$labels > 0)
    )
  ) |>
    dplyr::mutate(
      x = seq_along(label),
      value = format(value, big.mark = ",", trim = TRUE)
    )
  ggplot(stats, aes(x)) +
    geom_text(aes(y = 1, label = value), hjust = 0, vjust = 0, size = 4, fontface = "bold") +
    geom_text(aes(y = 0.7, label = label), hjust = 0, vjust = 1, size = 2.6, colour = "grey35") +
    scale_x_continuous(limits = c(1, nrow(stats) + 0.95), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 2)) +
    labs(title = name) +
    theme_void(base_size = 9) +
    theme(plot.title = element_text(face = "bold", size = 10, hjust = 0, margin = margin(b = 4)))
}

card_panel_list <- function(d, name) {
  inside <- function(x, y) {
    theme(
      legend.position = "inside",
      legend.position.inside = c(x, y),
      legend.justification = c(x, y),
      legend.direction = "vertical",
      legend.background = element_blank()
    )
  }
  axes <- theme(
    axis.line = cowplot::theme_cowplot(9)$axis.line,
    axis.ticks = cowplot::theme_cowplot(9)$axis.ticks
  )
  list(
    patchwork::free(card_header(d, name)),
    panel_scree(d, labels = c("counted", "not counted")) +
      labs(title = sprintf("Cross-validated eigenvalues: %d", d$pa$k)) +
      guides(fill = guide_legend(override.aes = list(size = 1.2))) +
      inside(1, 1),
    panel_survival(d, compact = TRUE) +
      labs(y = bquote("Survival score, " * beta == .(d$survival$guided$beta))) +
      guides(linetype = guide_legend(ncol = 1)) +
      inside(1, 1),
    panel_cut(d, show_cut = FALSE) +
      labs(
        y = "Loadings (sqrt scale)",
        title = sprintf(
          "Cut: %s of %s genes",
          format(sum(rowSums(d$cut$keep) > 0), big.mark = ","), format(d$p, big.mark = ",")
        )
      ) +
      inside(1, 1),
    panel_heatmap(
      d, strips = FALSE, legend = "right", fixed = FALSE, scale_name = expression(italic(r))
    ) +
      labs(title = "Gene correlation by module") +
      guides(fill = guide_colourbar(barwidth = unit(1.5, "mm"), barheight = unit(18, "mm"))) +
      theme(
        aspect.ratio = 1,
        legend.position = "inside",
        legend.position.inside = c(1.02, 0.5),
        legend.justification = c(0, 0.5),
        legend.background = element_blank()
      ) +
      axes,
    panel_eigengenes(d) +
      labs(title = "Eigengene |r|") +
      theme(aspect.ratio = 1, axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      axes,
    panel_sizes(d, both = FALSE, option = "viridis"),
    panel_frequency(d, compact = TRUE, quadrants = TRUE) +
      labs(
        title = sprintf("%s stability", d$stability$method),
        x = "Robust z-score of loading", y = "Selection frequency"
      ),
    panel_nystrom_by_m(d, option = "plasma") +
      labs(title = "Nystr\u00f6m agreement", y = "Correlation with exact"),
    panel_coverage(d, option = "mako") +
      labs(title = "Landmarks per module") +
      inside(1, 0.04) +
      theme(legend.text = element_text(size = 6.5, colour = "grey15"))
  )
}

card_design <- function() c(
  patchwork::area(1, 1, 1, 3),
  patchwork::area(2, 1), patchwork::area(2, 2), patchwork::area(2, 3),
  patchwork::area(3, 1), patchwork::area(3, 2), patchwork::area(3, 3),
  patchwork::area(4, 1), patchwork::area(4, 2), patchwork::area(4, 3)
)

card_theme <- function() {
  theme(
    legend.margin = margin(0, 0, 0, 0),
    legend.box.spacing = unit(3, "pt"),
    legend.text = element_text(size = 6, colour = "grey35"),
    legend.title = element_text(size = 6, colour = "grey35"),
    legend.key.size = unit(0.5, "lines"),
    legend.key.width = unit(0.7, "lines"),
    legend.key.spacing.y = unit(1, "pt")
  )
}

card_panels <- function(d, name = "") {
  purrr::set_names(card_panel_list(d, name), c("header", letters[1:9]))
}

card_rows <- list(count = c("a", "b", "c"), modules = c("d", "e", "f"), kernel = c("g", "h", "i"))

#' Assemble the 9-Panel Diagnostic Model Card
#'
#' @param d Diagnostic object computed by \code{\link{diag_compute}}.
#' @param name Character title for the card header.
#' @return A `patchwork` plot object.
#' @export
plot_model_card <- function(d, name = "") {
  patchwork::wrap_plots(card_panels(d, name), design = card_design(), heights = c(0.3, 1, 1.15, 1)) +
    patchwork::plot_annotation(tag_levels = list(c("", letters[1:9]))) &
    card_theme()
}

#' Plot a Single Row of the Diagnostic Model Card
#'
#' @param d Diagnostic object computed by \code{\link{diag_compute}}.
#' @param row Character row name: `"count"`, `"modules"`, or `"kernel"`.
#' @return A `patchwork` plot object.
#' @export
plot_card_row <- function(d, row = names(card_rows)) {
  tags <- card_rows[[match.arg(row)]]
  panels <- card_panels(d)[tags]
  if ("d" %in% tags) {
    panels$d <- panels$d + theme(plot.margin = margin(5.5, 26, 5.5, 5.5))
  }
  patchwork::wrap_plots(panels, nrow = 1) +
    patchwork::plot_annotation(tag_levels = list(tags)) &
    card_theme()
}

#' Plot a Single Panel of the Diagnostic Model Card
#'
#' @param d Diagnostic object computed by \code{\link{diag_compute}}.
#' @param panel Character panel tag (`"a"` through `"i"`).
#' @return A `ggplot` object.
#' @export
plot_card_panel <- function(d, panel = letters[1:9]) {
  panel <- match.arg(panel)
  p <- card_panels(d)[[panel]] + card_theme()
  if (panel == "d") {
    p <- p + theme(plot.margin = margin(5.5, 30, 5.5, 5.5))
  }
  p
}

#' Render and Save the 9-Panel Diagnostic Model Card
#'
#' @param d Diagnostic object computed by \code{\link{diag_compute}}.
#' @param name Character title for the card header.
#' @param out_name Optional output filename stem.
#' @param prefix Default filename prefix when `out_name` is `NULL`.
#' @param rows Logical; if `TRUE`, also saves each row separately.
#' @param panels Logical; if `TRUE`, also saves each panel separately.
#' @param save_fn Function with signature `(plot, name, width, height)` used to
#'   save figures (defaults to `save_figure` in the global environment if defined).
#' @export
draw_model_card <- function(d, name, out_name = NULL, prefix = "pkg", rows = FALSE,
                            panels = FALSE, save_fn = NULL) {
  if (is.null(save_fn)) {
    if (exists("save_figure", mode = "function", envir = .GlobalEnv)) {
      save_fn <- get("save_figure", mode = "function", envir = .GlobalEnv)
    } else {
      stop("`save_fn` must be provided or `save_figure()` must exist in .GlobalEnv.")
    }
  }
  card_file <- if (!is.null(out_name)) out_name else paste0(prefix, "_model_card")
  save_fn(plot_model_card(d, name), card_file, 180, 196)
  if (rows) {
    save_fn(plot_card_row(d, "count"), paste0(prefix, "_card_row1_count"), 180, 62)
    save_fn(plot_card_row(d, "modules"), paste0(prefix, "_card_row2_modules"), 180, 70)
    save_fn(plot_card_row(d, "kernel"), paste0(prefix, "_card_row3_kernel"), 180, 62)
  }
  if (panels) {
    for (panel in letters[1:9]) {
      save_fn(
        plot_card_panel(d, panel), paste0(prefix, "_card_panel_", panel),
        if (panel == "d") 70 else 62, 62
      )
    }
  }
}

#' Plot Method for rotatR Objects
#'
#' @param x A fitted `rotatR` object.
#' @param data Optional expression matrix `x` used to compute the full 9-panel
#'   Diagnostic Model Card via \code{\link{diag_compute}}. If omitted, plots the
#'   module size bar chart and loading distribution from the fitted object.
#' @param name Title for the diagnostic card.
#' @param ... Passed to \code{\link{diag_compute}}.
#' @export
plot.rotatR <- function(x, data = NULL, name = sprintf("rotatR (%s)", x$method), ...) {
  if (!is.null(data)) {
    d <- diag_compute(list(x = as.matrix(data), labels = NULL), ...)
    return(plot_model_card(d, name = name))
  }
  mod_ids <- sort(setdiff(unique(x$labels), 0L))
  sizes <- tibble::tibble(module = x$labels) |>
    dplyr::filter(module > 0L) |>
    dplyr::count(module) |>
    dplyr::mutate(name = factor(mname(module), rev(mname(mod_ids))))
  ggplot(sizes, aes(n, name)) +
    geom_col(fill = diag_kept, width = 0.75) +
    labs(x = "Genes", y = NULL, title = sprintf("%s: %d modules", name, length(mod_ids))) +
    cowplot::theme_cowplot(9)
}
