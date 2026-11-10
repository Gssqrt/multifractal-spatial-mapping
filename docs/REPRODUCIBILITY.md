# Reproducibility checklist

1. Archive original georeferenced raster metadata: sensor, date, pixel size, CRS, extent, segmentation algorithm and no-data rules.
2. Record nx, ny, `drop_border`, missing-mask semantics, identifier ordering, and complete centroid manifest.
3. Record FracLac and ImageJ version, scanning/box-counting settings, q range and spacing, grids/offsets, and export format.
4. Retain the **untouched** FracLac output AND the **cleaned numerical** dataset; document which 220-row-export metadata records were removed and why 219 remain in the study's dataset.
5. Check that the number of numeric records is `number_of_processed_partitions × rows_per_partition` and that each partition has the same ordered q sequence when appropriate.
6. Verify a one-to-one `partition` join; also store coordinates, bounds, original CRS and tile ID.
7. Archive fitted empirical/theoretical semivariograms, model comparison and LOOCV errors; inspect model assumptions and spatial extrapolation.
8. Save prediction and kriging-variance GeoTIFFs, software versions, timestamps and execution configuration.
9. Validate numerically against one original real-data FracLac export and published map before claiming equivalence with the article pipeline.

The source script uses `SSErr` as its variogram selection criterion; the new module defaults to `RMSE_CV` and must be configured to `SSErr` when reproducing historical outputs.
