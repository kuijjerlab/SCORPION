# Package the SERGIO simulation as a test fixture.
# Usage: Rscript data-raw/sergio/03_make_fixture.R <outDir>
# Writes: tests/testthat/fixtures/sergio_signed.rds

args <- commandArgs(trailingOnly = TRUE)
outDir <- if (length(args)) args[1] else "data-raw/sergio/out"

counts <- read.csv(file.path(outDir, "counts.csv"), row.names = 1, check.names = FALSE)
cells <- read.csv(file.path(outDir, "cells.csv"), stringsAsFactors = FALSE)
genes <- read.csv(file.path(outDir, "genes.csv"), stringsAsFactors = FALSE)
grn <- read.csv(file.path(outDir, "grn.csv"), stringsAsFactors = FALSE)
stopifnot(identical(colnames(counts), cells$cell_id), identical(rownames(counts), genes$gene))

load("data/scorpionTest.RData")
tfs <- genes$gene[genes$role != "target"]
prior <- scorpionTest$tf[scorpionTest$tf$source_genesymbol %in% tfs &
                           scorpionTest$tf$target_genesymbol %in% genes$gene, ]
prior$inGRN <- paste(prior$source_genesymbol, prior$target_genesymbol) %in% paste(grn$tf, grn$target)
rownames(prior) <- NULL
stopifnot(all(paste(grn$tf, grn$target) %in% paste(prior$source_genesymbol, prior$target_genesymbol)))

fixture <- list(
  gex = Matrix::Matrix(as.matrix(counts), sparse = TRUE),
  metadata = cells,
  grn = grn,
  genes = genes,
  prior = prior,
  info = "SERGIO steady-state simulation (data-raw/sergio); GRN sampled from scorpionTest$tf signed edges (mor)"
)
dir.create("tests/testthat/fixtures", showWarnings = FALSE, recursive = TRUE)
saveRDS(fixture, "tests/testthat/fixtures/sergio_signed.rds", compress = "xz")
cat(sprintf("Fixture: %d genes x %d cells | GRN edges %d | prior edges %d (%d in GRN) | %.0f KB\n",
            nrow(counts), ncol(counts), nrow(grn), nrow(prior), sum(prior$inGRN),
            file.size("tests/testthat/fixtures/sergio_signed.rds") / 1024))
