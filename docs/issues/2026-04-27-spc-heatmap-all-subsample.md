# Issue: Add `subsample` Feature Filtering to `spc_heatmap_all`

## Summary
Add a `subsample` argument to `spc_heatmap_all` (default `NULL`) so users can reduce plotted features when generating the all-genes SPC heatmap.

When `subsample` is provided, rows should be restricted to:
- the top `subsample` features ranked by median absolute deviation (MAD), and
- all features in the input signature,

using the union of those two sets.

## Motivation
For large expression matrices, plotting all genes can be slow and visually noisy. A variation-based filter keeps informative genes while guaranteeing signature genes remain visible.

## Scope
- Add and document `subsample` in `spc_heatmap_all`.
- Implement subsampling in heatmap prep logic.
- Validate `subsample` input:
  - single numeric integer value,
  - greater than `0`,
  - strictly less than `nrow(eset)`.
- Keep existing ordering semantics from SPC output.
- Add tests for signature retention under subsampling.

## Implementation Notes
- Compute MAD per feature from `Biobase::exprs(eset)`.
- Select top `subsample` features by descending MAD.
- Build `keep_genes <- union(top_mad_genes, signature_genes)`.
- Restrict heatmap rows to `keep_genes` while preserving SPC-derived row order.

## Acceptance Criteria
- `spc_heatmap_all(eset, spc_out, subsample = NULL)` behaves as before.
- `spc_heatmap_all(..., subsample = k)` returns a heatmap built from `topMAD(k) ∪ signature`.
- Signature genes are included even if not in top MAD set.
- Invalid `subsample` values fail fast.

## Verification
- Added regression test in `tests/testthat/test-signature_projection_contributors.R`:
  - confirms signature genes are retained with `subsample`,
  - checks row count reflects union behavior.

## Related Files
- `R/spc_heatmaps.R`
- `man/spc_heatmap_all.Rd`
- `tests/testthat/test-signature_projection_contributors.R`
