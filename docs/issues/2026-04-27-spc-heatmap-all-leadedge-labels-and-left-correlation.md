# Issue: Improve `spc_heatmap_all` Annotation Readability

## Summary
Update `spc_heatmap_all()` to improve annotation readability by:
- adding optional labels for top leading-edge genes, and
- moving the correlation barplot to the left side of the heatmap.

## Motivation
- Leading-edge tick marks (`yes/no`) are informative but can be hard to interpret without explicit gene labels.
- Correlation bars are easier to scan when placed on the left margin before entering the heatmap body.

## Scope
- Add optional controls to label top leading-edge genes:
  - `leadedge_label_n` (default `0`, disabled)
  - `leadedge_label_side` (`"left"` or `"right"`, default `"right"`)
- Render labels with offset + connector lines for readability.
- Move correlation row annotation from right-side annotation block to `left_annotation`.
- Preserve existing defaults/behavior when new options are not used.

## Implementation Notes
- In `spc_heatmap_all`, compute candidate genes from `prep$ks_hits`.
- Rank by absolute correlation (`abs(score_cor)`) and keep top `leadedge_label_n`.
- Add labels with `ComplexHeatmap::anno_mark(...)`:
  - `link_width = unit(4, "mm")`
  - `extend = unit(3, "mm")`
- Keep `genes` and `leadedge` categorical annotations on the right.
- Move correlation bars into `left_row_ha` passed to `Heatmap(..., left_annotation = left_row_ha)`.

## Acceptance Criteria
- With defaults, `spc_heatmap_all()` renders as before except correlation bars are shown on the left.
- With `leadedge_label_n > 0`, top leading-edge genes are labeled with non-overlapping offset marks and connector lines.
- Label side is configurable via `leadedge_label_side`.

## Changed Files
- `R/spc_heatmaps.R`
- `man/spc_heatmap_all.Rd`
- `NAMESPACE`
