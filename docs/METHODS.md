# Methods / Descrição metodológica

## Scope
This repository does **not** compute multifractal spectra. Those are computed outside the repository by FracLac for ImageJ. Scripts in R spatially index results, calculate spectral descriptors, and interpolate values using ordinary kriging. The raw Java FracLac implementation is not redistributed.

## 1. Partitioning
Georeferenced single-band masks are divided into an nx × ny grid; each tile receives an identifier and its centroid coordinates in the raster's original CRS. When `drop_border=TRUE`, remainder pixels at the right/bottom edges are excluded. The original analysis used a configurable grid, commonly 17 × 17; these choices affect the support of descriptor measurements.

## 2. FracLac output and spectrum indexing
Each partition is analyzed separately in FracLac/ImageJ. The processed numeric records per partition are assigned the corresponding partition ID in the same sequence as the image processing order. The study's particular cleaned export used 219 records per partition: **not a general FracLac output guarantee**. Textual/metadata records must be explicitly removed before running the indexing stage, with the cleaning procedure archived and verified for each export.

## 3. Extracted descriptors
For each partition, using FracLac's alpha(q), f(alpha(q)), and D(q):
- `d.alpha` = max(alpha) - min(alpha)
- `df.alpha` = f(alpha_at_max_alpha) - f(alpha_at_min_alpha); signed convention preserved from the research script
- `alpha.0` = alpha at argmax f(alpha)
- `D0`, `D1`, `D2` = D(q=0,1,2)
- `r` = (max(alpha) - alpha.0)/(alpha.0 - min(alpha))

Nonfinite spectra are rejected; absent or repeated matching q values yield NA for the corresponding D(q), not silently chosen values. A degenerate r denominator yields NA. These definitions should be reported explicitly in the paper.

## 4. Spatial attribution
Join descriptor values with the centroid manifest on the *partition key*, never row order. Reproject centroids from source CRS to EPSG:4326 for longitude/latitude. Cardinality and missing keys are checked. A tile descriptor is an *areal summary* assigned to the tile centroid and should not be interpreted as a directly measured pixel-scale descriptor.

## 5. Geostatistical modeling
Point-based empirical semivariograms are constructed in a projected metric coordinate system (UTM by default). Candidate isotropic models: spherical, exponential, Gaussian and Matérn. Selection can be based on leave-one-out cross-validation RMSE (new default) or variogram SSErr (original script default). Ordinary kriging estimates continuous descriptor fields on a defined grid. Both prediction and kriging variance are exported. The suggested 100 × 100 grid is a visualization/interpolation resolution, not an increase in the spatial resolution of the original multifractal observations.

## 6. Limitations
- No assertion that ordinary kriging assumptions (spatial stationarity, isotropy) hold without investigation.
- Tile centroids simplify areal support to point locations; change-of-support effects must be considered.
- LOOCV across spatially correlated neighbors may be optimistic; spatial-block validation is preferable for stronger predictive claims.
- Interpolated pixels are correlated; clustering them is descriptive, not independent replication.
- **No end-to-end validation against original satellite scenes/FracLac exports has been run.** Provide anonymized input/output pairs to establish numerical reproducibility.
