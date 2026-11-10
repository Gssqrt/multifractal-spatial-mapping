# What changed relative to the research script

- Kept an unchanged copy of `Metotologia_tese.R` in `legacy/`.
- Split independent tasks into functions with explicit inputs and outputs.
- Retained FracLac-generated spectra as external input rather than substituting Python implementations.
- Retained alpha-span, signed vertical f-alpha difference, alpha0, D0/D1/D2 and r definitions.
- Made numeric record count per partition configurable (the study used 219 after removal of textual records).
- Added strict checks for record count, q sequence, duplicated/missing IDs and source CRS.
- Replaced ambiguous empty-tile classification in original script with explicit background vs NA logic.
- No automatic assignment of CRS to a raster missing georeference.
- Writes geostatistical model diagnostics even when SSErr is the chosen criterion.
- New default model selection: LOOCV RMSE; original research script used SSErr. To recreate original behavior, set `criterion='SSErr'`.
- Matrices require aligned raster geometries; they are not silently resampled.
- Optional cluster step retains Ward.D2, but no longer hardcodes a dendrogram display inconsistent with assigned clusters.
- **No claim of full numerical equivalence or article-result reproduction** until tested against supplied original FracLac spreadsheets, georeferenced imagery, and output maps.
