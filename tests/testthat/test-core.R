make_synthetic_modules <- function(n = 80, p = 300, k_true = 3, mod_size = 60, snr = 2.0, seed = 42) {
  set.seed(seed)
  x <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(x) <- paste0("Gene", seq_len(p))
  rownames(x) <- paste0("Sample", seq_len(n))
  true_labels <- integer(p)
  for (m in seq_len(k_true)) {
    idx <- ((m - 1) * mod_size + 1):(m * mod_size)
    z <- rnorm(n)
    x[, idx] <- x[, idx] + snr * outer(z, rep(1, length(idx)))
    true_labels[idx] <- m
  }
  list(x = x, labels = true_labels, k = k_true)
}

test_that("vsp recovers planted modules with high accuracy", {
  sim <- make_synthetic_modules(n = 80, p = 300, k_true = 3, mod_size = 60, snr = 2.0)
  fit <- rotatR(sim$x, method = "vsp", k = 3, min_size = 15, seed = 1)

  expect_s3_class(fit, "rotatR")
  expect_equal(length(fit$labels), ncol(sim$x))
  expect_equal(ncol(fit$eigengenes), length(setdiff(unique(fit$labels), 0L)))
  expect_gte(length(setdiff(unique(fit$labels), 0L)), 3L)

  for (m in 1:3) {
    idx <- which(sim$labels == m)
    tab <- table(fit$labels[idx])
    tab <- tab[names(tab) != "0"]
    expect_gt(max(tab) / length(idx), 0.85)
  }
})

test_that("kvsp and kvsp_nys recover modules and agree closely", {
  sim <- make_synthetic_modules(n = 80, p = 300, k_true = 3, mod_size = 60, snr = 2.0)
  fit_exact <- rotatR(sim$x, method = "kvsp", k = 3, beta = 2, z = 2.5, min_size = 15, seed = 1)
  fit_nys <- rotatR(sim$x, method = "kvsp_nys", k = 3, beta = 2, m = 150, z = 2.5, min_size = 15, seed = 1)

  expect_s3_class(fit_exact, "rotatR")
  expect_s3_class(fit_nys, "rotatR")
  expect_gte(sum(fit_exact$labels > 0L), 120L)
  expect_gte(sum(fit_nys$labels > 0L), 120L)

  both_assigned <- (fit_exact$labels > 0L) & (fit_nys$labels > 0L)
  expect_gt(mean(both_assigned), 0.5)
})

test_that("dimension selection eigcv_k detects planted rank", {
  sim <- make_synthetic_modules(n = 100, p = 300, k_true = 4, mod_size = 50, snr = 2.2)
  k_cv <- eigcv_k(sim$x, cap = 15, seed = 1)

  expect_gte(k_cv, 3L)
  expect_lte(k_cv, 6L)
})

test_that("module resolution merges redundant axes and drops tiny clusters", {
  sim <- make_synthetic_modules(n = 80, p = 200, k_true = 2, mod_size = 50, snr = 2.5)
  v <- vsp_loadings(sim$x, k = 4)
  v_cut <- fdr_cut_loadings(v, q = 0.01)

  labels <- resolve_modules(sim$x, v_cut, min_size = 20, merge_cor = 0.85)
  mod_ids <- setdiff(unique(labels), 0L)
  expect_true(all(table(labels[labels > 0L]) >= 20L))
  expect_lte(length(mod_ids), 3L)
})
