# rotatR

**Spectral and Varimax Co-expression Module Discovery**

`rotatR` provides fast, parameter-free co-expression module discovery via truncated singular value decomposition (SVD), leverage-score Nyström kernel approximation, and varimax rotation, with a drop-in compatibility interface for `WGCNA` workflows.

## Installation

```r
# install.packages("remotes")
remotes::install_github("cashoes/rotatR")
```

## Overview

`rotatR` implements a unified 4-step spectral module discovery pipeline for sample-by-gene expression matrices (`x`, $n \times p$):

1. **Similarity / Kernel Construction (`find_beta`)**:
   - Linear standardized correlation (`vsp_modules`)
   - Exact powered correlation kernel $|R|^\beta$ (`kvsp_modules`, default $\beta = 2$)
   - Leverage-score Nyström kernel approximation (`kvsp_nys_modules`, $O(p m^2)$ memory/time without forming the $p \times p$ matrix)
2. **Dimension Selection (`find_k`, `eigcv_k`)**:
   - Bi-cross-validated eigenvalue significance testing (Gabriel / Owen–Perry split-sample test) to separate biological module rank from diffuse background correlation.
3. **Structure & Loading Sparsification (`fdr_cut_loadings`, `stability_select`)**:
   - Analytical Benjamini–Hochberg false discovery rate (1% BH-FDR) thresholding on robust median/MAD $z$-scores for linear Varimax PCA (`vsp_modules`), or Meinshausen–Bühlmann subsampling stability selection for kernel Varimax PCA (`kvsp_modules`, `kvsp_nys_modules`).
4. **Module Resolution (`resolve_modules`)**:
   - Winner-take-all (`argmax`) assignment of surviving loadings, iterative module eigengene merging ($r \ge 0.85$), and minimum module size pruning ($\ge 20$ genes).

## Quick Start

```r
library(rotatR)

# Linear Varimax PCA module discovery (1% BH-FDR cut, O(n^2 p) without forming R)
fit <- vsp_modules(X)
table(fit$labels)

# Leverage-Nystrom Kernel Varimax PCA (beta = 2, m = 256 landmarks, stability-selected)
fit_nys <- kvsp_nys_modules(X, beta = 2, m = 256)

# Drop-in WGCNA-compatible output ($colors, $MEs, $kME, $powerEstimate)
net <- blockwiseModules_vsp(X)
head(net$MEs)

# 9-panel diagnostic model card
diag <- diag_compute(X, method = "kvsp_nys")
p <- plot_model_card(diag)
```

## License

MIT © Casey Shannon
