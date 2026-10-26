test_that("S3 methods and as_wgcna() work cleanly", {
  set.seed(101)
  n <- 60
  p <- 180
  x <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(x) <- paste0("G", seq_len(p))
  rownames(x) <- paste0("S", seq_len(n))
  for (m in 1:2) {
    idx <- ((m - 1) * 50 + 1):(m * 50)
    x[, idx] <- x[, idx] + 2.2 * outer(rnorm(n), rep(1, length(idx)))
  }

  fit <- rotatR(x, method = "vsp", k = 2, min_size = 15, seed = 1)

  expect_output(print(fit), "rotatR fit")

  s <- summary(fit)
  expect_s3_class(s, "summary.rotatR")
  expect_output(print(s), "Summary of rotatR")
  expect_true(is.data.frame(s$modules))

  df <- as.data.frame(fit)
  expect_equal(nrow(df), p)
  expect_true(all(c("gene", "module", "loading", "abs_loading") %in% colnames(df)))

  preds <- predict(fit, x[1:10, ])
  expect_equal(nrow(preds), 10L)
  expect_equal(ncol(preds), ncol(fit$eigengenes))

  w <- as_wgcna(fit, use_color_names = FALSE)
  expect_equal(w$colors, fit$labels)
  expect_equal(nrow(w$MEs), n)

  p_obj <- plot(fit)
  expect_s3_class(p_obj, "ggplot")
})
