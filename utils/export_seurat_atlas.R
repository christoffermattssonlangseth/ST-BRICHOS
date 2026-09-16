#!/usr/bin/env Rscript
# Export a Seurat .RDS reference to a Matrix-Market bundle that AnnData can
# ingest, WITHOUT SeuratDisk/h5Seurat (which are fragile across Seurat
# versions). Works with Seurat v4 and v5 assay layouts.
#
# Usage:
#   Rscript utils/export_seurat_atlas.R <in.RDS> <outdir> [assay] [layer]
#
# Defaults: assay="RNA", layer="counts" (raw integer counts — what
# cell2location's reference regression model needs).
#
# Writes into <outdir>/:
#   matrix.mtx.gz      genes x cells sparse counts (Matrix::writeMM)
#   features.tsv.gz    gene symbols (rownames of the assay)
#   barcodes.tsv.gz    cell barcodes (colnames)
#   metadata.csv.gz    full cell metadata (obs) incl. cluster/state labels
# and prints the metadata column names + a preview so you can pick the
# annotation column to deconvolve on.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("usage: export_seurat_atlas.R <in.RDS> <outdir> [assay=RNA] [layer=counts]")
}
in_rds <- args[[1]]
outdir <- args[[2]]
assay  <- ifelse(length(args) >= 3, args[[3]], "RNA")
layer  <- ifelse(length(args) >= 4, args[[4]], "counts")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
message("reading ", in_rds, " ...")
obj <- readRDS(in_rds)
message("class: ", paste(class(obj), collapse = ", "))

# Some deposits store a list or an older object; coerce if needed.
if (!inherits(obj, "Seurat")) {
  stop("object is not a Seurat object; got class ", paste(class(obj), collapse=","))
}

message("assays available: ", paste(Assays(obj), collapse = ", "))
if (!assay %in% Assays(obj)) {
  assay <- Assays(obj)[[1]]
  message("requested assay not found; falling back to '", assay, "'")
}
DefaultAssay(obj) <- assay

# Robust counts fetch across Seurat 4/5.
counts <- tryCatch(
  GetAssayData(obj, assay = assay, layer = layer),
  error = function(e) GetAssayData(obj, assay = assay, slot = layer)
)
if (is.null(counts) || nrow(counts) == 0) {
  stop("empty counts for assay=", assay, " layer=", layer)
}
counts <- as(counts, "CsparseMatrix")
message("counts: ", nrow(counts), " genes x ", ncol(counts), " cells")

md <- obj@meta.data
md <- cbind(barcode = rownames(md), md)

# --- write bundle ----------------------------------------------------------
mtx_path <- file.path(outdir, "matrix.mtx")
Matrix::writeMM(counts, mtx_path)
system2("gzip", c("-f", shQuote(mtx_path)))

writeLines(rownames(counts), gzfile(file.path(outdir, "features.tsv.gz")))
writeLines(colnames(counts), gzfile(file.path(outdir, "barcodes.tsv.gz")))
write.csv(md, gzfile(file.path(outdir, "metadata.csv.gz")), row.names = FALSE)

# --- report so the user can choose the annotation column -------------------
message("\n=== metadata columns (candidate cell-type/state labels) ===")
for (cn in colnames(md)) {
  v <- md[[cn]]
  if (is.factor(v) || is.character(v)) {
    nu <- length(unique(v))
    if (nu >= 2 && nu <= 300) {
      message(sprintf("  %-32s %d levels", cn, nu))
    }
  }
}
message("\nwrote bundle to: ", normalizePath(outdir))
