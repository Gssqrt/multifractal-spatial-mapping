# Case study: spatial attribution of image-partition descriptors (2021-05-16)

This directory provides a real, already-processed example showing how descriptors calculated for numbered image sections are connected to their spatial centroids.

**Provenance:** These supplied tables contain fields such as `selected_position`, `foreground` and `status` associated with an experimental analysis workflow. They must **not** be described as unmodified, direct FracLac/ImageJ exports. They demonstrate the indexing, descriptor aggregation, and spatial attribution data structure. The original satellite image, partition TIFFs and original FracLac output are **not** included.

## Files

| File | Role |
|---|---|
| `data/centroids.csv` | 289 numbered sections of a 17 × 17 partition grid, their centroids and spatial extents |
| `data/indexed_spectra.csv` | 50,820 spectral rows for 242 sections with numerical results (210 rows per valid section) |
| `data/multifractal_descriptors.csv` | One row per section, with descriptors and processing status (289 records) |
| `data/spatial_descriptors.csv` | Join of section coordinates and descriptor/status columns (289 records) |

`partition` is the shared join key. All four files use semicolon separators and comma decimal marks (`read.csv2` in R). Keep partition IDs as strings (`01`, `02`, ...), not integers. The `arquivo` values in the spatial tables have been changed from private Windows paths to *relative illustrative paths* under `partitions/`; those TIFF files are not included.

## Data integrity checks

- 289 unique partition IDs in the centroid, descriptor, and joined tables.
- 242 valid sections in `indexed_spectra.csv`, with 210 numerical spectral rows each.
- The remaining 47 sections have `status=ERRO` in the parameter tables and no numeric spectra. Do not interpret their zero-filled descriptors as valid measured zeros; exclude such records from modeling using `status == "OK"`.
- These sample files are intended to show *post-processing and linkage*, not to validate equivalence of any spectral-computation implementation with FracLac.

## Inspect in R

```r
folder <- "examples/case_20210516/data"
centroids <- read.csv2(file.path(folder, "centroids.csv"), colClasses = c(partition = "character"))
spectra <- read.csv2(file.path(folder, "indexed_spectra.csv"), colClasses = c(partition = "character"))
params <- read.csv2(file.path(folder, "multifractal_descriptors.csv"), colClasses = c(partition = "character"))
spatial <- read.csv2(file.path(folder, "spatial_descriptors.csv"), colClasses = c(partition = "character"))

stopifnot(nrow(centroids) == 289, nrow(spectra) == 50820,
          length(unique(spectra$partition)) == 242,
          identical(sort(centroids$partition), sort(params$partition)),
          identical(sort(centroids$partition), sort(spatial$partition)))

valid <- subset(spatial, status == "OK")
stopifnot(nrow(valid) == 242)
```

## Scientific cautions

Spatially interpolated descriptors are **estimates for locations between partition centroids**, not direct pixel-level fractal measurements. The tables alone are not sufficient to reproduce the original satellite segmentation or FracLac calculation; those steps require the original images, settings, and processing environment.

Before public release, confirm that redistribution of the provided scientific data and their coordinates is permitted under the project's data-sharing and authorship arrangements.
