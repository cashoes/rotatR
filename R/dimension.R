#' Dimension Selection for Spectral and Varimax Co-expression Discovery
#'
#' @name dimension
NULL

#' Cross-Projected Eigenvalue Profile (eigcv)
#'
#' Evaluates out-of-sample eigenvalues via 2-fold sample splitting on the linear
#' spectrum (Chen et al. 2021, arXiv:2108.03336). Under the null hypothesis of
#' independent/isotropic features, out-of-sample quadratic forms have closed-form
#' CLT standard error \eqn{\sqrt{(2 / n) / \text{valid\_splits}}}.
#'
#' @param x Numeric matrix (samples in rows, genes in columns).
#' @param cap Numeric/integer maximum rank to evaluate (default `40`). Use `Inf`
#'   for uncapped half-sample rank.
#' @param seed Integer random seed (default `0`).
#' @param alpha Numeric one-sided significance level (default `0.01`).
#' @param n_splits Integer number of balanced 2-fold splits (default `2`).
#' @return A list with `k` (number of surviving components, at least `2L`),
#'   `lambda` (cross-projected eigenvalue vector), and `threshold` (null cutoff).
#' @export
eigcv_profile <- function(x, cap = 40, seed = 0, alpha = 0.01, n_splits = 2) {
  empty <- list(k = 2L, lambda = numeric(0), threshold = NA_real_)
  n <- nrow(x)
  p <- ncol(x)
  if (n < 6L || p < 2L) return(empty)

  r <- min(if (is.finite(cap)) as.integer(cap) else .Machine$integer.max, floor(n / 2) - 1L, p - 1L)
  if (r < 2L) return(empty)

  lambda_accum <- numeric(r)
  valid_splits <- 0L

  for (split_idx in seq_len(n_splits)) {
    set.seed(seed + 10007L * split_idx)
    n1 <- floor(n / 2)
    n2 <- n - n1
    perm <- sample.int(n)
    idx1 <- perm[seq_len(n1)]
    idx2 <- perm[(n1 + 1L):n]

    z1 <- standardize(x[idx1, , drop = FALSE])
    z2 <- standardize(x[idx2, , drop = FALSE])

    # Fold 1 -> 2
    g1 <- tcrossprod(z1) / (n1 - 1)
    eg1 <- eigen(g1, symmetric = TRUE)
    pos1 <- which(eg1$values[seq_len(r)] > 1e-10)
    if (length(pos1) == 0L) next
    r1 <- min(r, length(pos1))
    u1 <- eg1$vectors[, seq_len(r1), drop = FALSE]
    d1 <- sqrt(eg1$values[seq_len(r1)])
    v1 <- crossprod(z1, sweep(u1, 2, d1 * sqrt(n1 - 1), "/"))
    y2 <- (z2 %*% v1) / sqrt(n2 - 1)
    lambda_12 <- colSums(y2^2)

    # Fold 2 -> 1
    g2 <- tcrossprod(z2) / (n2 - 1)
    eg2 <- eigen(g2, symmetric = TRUE)
    pos2 <- which(eg2$values[seq_len(r)] > 1e-10)
    if (length(pos2) == 0L) next
    r2 <- min(r1, length(pos2))
    u2 <- eg2$vectors[, seq_len(r2), drop = FALSE]
    d2 <- sqrt(eg2$values[seq_len(r2)])
    v2 <- crossprod(z2, sweep(u2, 2, d2 * sqrt(n2 - 1), "/"))
    y1 <- (z1 %*% v2) / sqrt(n1 - 1)
    lambda_21 <- colSums(y1^2)

    lambda_accum[seq_len(r2)] <- lambda_accum[seq_len(r2)] + 0.5 * (lambda_12[seq_len(r2)] + lambda_21[seq_len(r2)])
    valid_splits <- valid_splits + 1L
  }

  if (valid_splits == 0L) return(empty)
  lambda_cv <- lambda_accum / valid_splits

  se_null <- sqrt((2 / n) / valid_splits)
  z_scores <- (lambda_cv - 1) / se_null
  p_vals <- stats::pnorm(z_scores, lower.tail = FALSE)

  pass_seq <- (p_vals <= alpha) & (lambda_cv > 1)
  k_seq <- if (all(pass_seq)) length(lambda_cv) else which(!pass_seq)[1L] - 1L

  list(
    k = max(2L, k_seq),
    lambda = lambda_cv,
    threshold = 1 + stats::qnorm(1 - alpha) * se_null
  )
}

#' Select Latent Dimension via Cross-Projected Eigenvalues (eigcv)
#'
#' @inheritParams eigcv_profile
#' @return Integer selected rank \eqn{k \ge 2}.
#' @export
eigcv_k <- function(x, cap = 40, seed = 0, alpha = 0.01, n_splits = 2) {
  eigcv_profile(x, cap, seed, alpha, n_splits)$k
}

#' Analytical Moment of Absolute Sample Correlation Under Null
#'
#' Computes \eqn{E[|r|^\beta]} for two independent genes on \code{n} samples,
#' where \eqn{r^2 \sim \text{Beta}(1/2, (n - 2)/2)}.
#'
#' @param beta Numeric exponent (\eqn{\beta > 0}).
#' @param n Integer sample size (\eqn{n > 3}).
#' @return Numeric expected value \eqn{E[|r|^\beta]}.
#' @export
null_abs_r_moment <- function(beta, n) {
  exp(
    lgamma((beta + 1) / 2) + lgamma((n - 1) / 2) -
      lgamma(1 / 2) - lgamma((n - 1 + beta) / 2)
  )
}

#' How Strongly Each Linear Component Survives in the Kernel
#'
#' For each of the \code{k} varimax-rotated linear components, takes the genes
#' where it has the largest absolute loading (up to \code{max_genes}) and compares
#' its kernel energy \eqn{v' |R|^\beta v} with the null expectation for
#' independent genes \eqn{1 + ((\sum v)^2 - 1) E[|r|^\beta]}.
#'
#' @param x Numeric matrix (samples in rows, genes in columns).
#' @param k Integer number of linear components.
#' @param beta Numeric soft-threshold power.
#' @param max_genes Integer maximum genes per component (default `500`).
#' @return Numeric vector of survival ratios of length \code{max(k, 2)}.
#' @export
kernel_survival_scores <- function(x, k, beta, max_genes = 500) {
  y <- abs(vsp_loadings(x, max(k, 2)))
  owner <- max.col(y, ties.method = "first")
  moment <- null_abs_r_moment(beta, nrow(x))
  vapply(seq_len(ncol(y)), \(j) {
    own <- which(owner == j)
    if (length(own) < 2) {
      return(NA_real_)
    }
    own <- own[order(y[own, j], decreasing = TRUE)]
    own <- own[seq_len(min(max_genes, length(own)))]
    v <- y[own, j] / sqrt(sum(y[own, j]^2))
    r <- stats::cor(x[, own, drop = FALSE])
    r[is.na(r)] <- 0
    energy <- sum(v * as.numeric((abs(r)^beta) %*% v))
    energy / (1 + (sum(v)^2 - 1) * moment)
  }, 0)
}

#' Select Kernel Dimension Among Surviving Linear Components
#'
#' @inheritParams kernel_survival_scores
#' @param threshold Numeric minimum survival score (default `1.1`).
#' @return Integer number of surviving components, between `2` and `k`.
#' @export
kernel_survival_k <- function(x, k, beta, threshold = 1.1, max_genes = 500) {
  scores <- kernel_survival_scores(x, k, beta, max_genes)
  survivors <- sum(scores > threshold, na.rm = TRUE)
  as.integer(min(k, max(min(k, 2L), survivors)))
}
