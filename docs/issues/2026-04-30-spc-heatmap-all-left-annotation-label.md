# Issue: spc_heatmap_all left annotation label should read "gene correlation"

## Summary
In `spc_heatmap_all`, the left row annotation barplot was not presenting the intended descriptive label and users only saw the right-side categorical annotation label (`genes`).

## Requested behavior
Display `gene correlation` for the left-side annotation in the all-genes SPC heatmap.

## Root cause
The left annotation used a generic internal name and had annotation names hidden (`show_annotation_name = FALSE`), so the intended left-side descriptor was not visible.

## Fix implemented
- Renamed the left row annotation key to `gene correlation`.
- Enabled annotation name display for the left row annotation so that label renders.

## Files changed
- `R/spc_heatmaps.R`

## Notes
This change is scoped to `spc_heatmap_all` and does not alter correlation values or row ordering logic.
