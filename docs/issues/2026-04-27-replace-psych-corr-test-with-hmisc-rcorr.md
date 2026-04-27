# Issue: Replace `psych::corr.test` with `Hmisc::rcorr`

> Status: superseded on April 27, 2026. The implementation was reverted to
> `psych::corr.test` to avoid unnecessary all-pairs correlation work in `rcorr`.

## Summary
Switch the gene-wise correlation and p-value computation in
`signature_projection_contributors()` from `psych::corr.test` to `Hmisc::rcorr`.

## Motivation
- Align with a widely used correlation utility in R/Bioconductor workflows.
- Keep a single function call that returns both correlation coefficients and p-values.
- Simplify dependency direction for this package's core correlation path.

## Scope
- Update correlation implementation in `R/signature_projection_contributors.R`.
- Update imports and dependency metadata:
  - remove `psych::corr.test`,
  - add `Hmisc::rcorr`.
- Update user-facing documentation for `cor_method`.

## Implementation Details
- Build one matrix with signature score + gene expression:
  - `corr_input <- cbind(sig_score = eset$sig_score, t(exprs(eset)))`
- Run:
  - `COR <- Hmisc::rcorr(corr_input, type = cor_method)`
- Extract score-vs-gene outputs:
  - `COR$r["sig_score", featureNames(eset)]` for correlations
  - `COR$P["sig_score", featureNames(eset)]` for p-values
- Persist into `fData(eset)$score_cor` and `fData(eset)$pval_cor`.

## Acceptance Criteria
- `signature_projection_contributors()` still returns:
  - `score_cor` column with gene-vs-signature correlations,
  - `pval_cor` column with corresponding p-values.
- Sorting and downstream heatmap/KS behavior remains unchanged.
- Package metadata reflects dependency switch from `psych` to `Hmisc`.

## Changed Files
- `R/signature_projection_contributors.R`
- `DESCRIPTION`
- `NAMESPACE`
- `man/signature_projection_contributors.Rd`
