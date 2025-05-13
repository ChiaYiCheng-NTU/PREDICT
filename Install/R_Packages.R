cran_packages <- c(
  "askpass", "cli", "colorspace", "crayon", "fansi", "farver", "ggplot2",
  "glue", "gtable", "httr", "isoband", "jsonlite", "labeling", "lifecycle",
  "magrittr", "matrixStats", "mime", "munsell", "openssl", "pillar",
  "pkgconfig", "R6", "RColorBrewer", "Rcpp", "RcppThread", "rlang", "scales",
  "sys", "tibble", "utf8", "vctrs", "viridisLite", "withr", "yaml"
)

bioc_packages <- c(
  "BiocGenerics", "BiocManager", "BiocVersion", "Biostrings", "GenomeInfoDb",
  "GenomeInfoDbData", "IRanges", "MatrixGenerics", "S4Vectors", "XVector", "universalmotif",
  "UCSC.utils", "zlibbioc"
)

# Set CRAN mirror
options(repos = "https://cloud.r-project.org")

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}
invisible(sapply(cran_packages, install_if_missing))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
invisible(sapply(bioc_packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE)
  }
}))
