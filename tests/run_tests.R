# Run: Rscript tests/run_tests.R (base R only; spatial integration requires packages)
root <- normalizePath(file.path(dirname(sub('^--file=', '', grep('^--file=', commandArgs(), value = TRUE)[1])), '..'), mustWork = TRUE)
for (f in sort(list.files(file.path(root,'R'), pattern = '\\.R$', full.names = TRUE))) source(f)
# 219 numerical records per partition, matching original study convention.
q <- seq_len(219) - 110
s <- data.frame(q = rep(q,2), a = rep(seq(0.8,1.8,length.out=219),2),
                f.a. = rep(1 - (q/200)^2, 2), d.q. = rep(1.3 + q/500, 2))
indexed <- index_spectra(s, c('01','02'), rows_per_partition=219)
stopifnot(nrow(indexed)==438, all(indexed$partition[1:219]=='1'), all(indexed$partition[220:438]=='2'))
d <- extract_descriptors(indexed)
stopifnot(nrow(d)==2, all(abs(d$d.alpha - 1) < 1e-12),
          all(abs(d$D0 - 1.3) < 1e-12), all(is.finite(d$alpha.0)))
wrong <- try(index_spectra(s[-1, ], c('01','02'), 219), silent=TRUE)
stopifnot(inherits(wrong, 'try-error'))
# Mismatched q sequence should be rejected despite correct row count.
s_bad <- s; s_bad$q[220] <- 1000
stopifnot(inherits(try(index_spectra(s_bad,c('01','02'),219),silent=TRUE),'try-error'))
# Sign of Delta f follows original, even when nonzero.
small <- data.frame(partition='01', q=c(0,1,2), a=c(1,2,3), f.a.=c(1,4,2), d.q.=c(1.2,1.1,1.0))
stopifnot(abs(extract_descriptors(small)$df.alpha - 1) < 1e-12)
cat('PASS: indexing, bounds, descriptor extraction, sign, q-grid consistency\n')
if (requireNamespace('sf',quietly=TRUE)) {
 cent <- data.frame(partition=c('01','02'), centroid_x=c(-39,-38),centroid_y=c(-10,-10),
                    crs='EPSG:4326', salvo=TRUE)
 merged <- join_descriptors_to_centroids(cent, d)
 stopifnot(nrow(merged)==2, all(merged$longitude==c(-39,-38)))
 stopifnot(inherits(try(join_descriptors_to_centroids(cent[1,,drop=FALSE],d),silent=TRUE),'try-error'))
 cat('PASS: spatial join and missing-partition rejection\n')
} else cat('SKIP: sf not installed (spatial tests)\n')
