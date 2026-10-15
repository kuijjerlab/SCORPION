# Build a signed, acyclic GRN for SERGIO from the scorpionTest TF prior (mor column).
# Usage: Rscript data-raw/sergio/01_build_grn.R <outDir>
# Writes: targets.csv, regs.csv (SERGIO inputs), genes.csv, grn.csv (ground truth)

args <- commandArgs(trailingOnly = TRUE)
outDir <- if (length(args)) args[1] else "data-raw/sergio/out"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
set.seed(20261014)

load("data/scorpionTest.RData")
prior <- scorpionTest$tf
signed <- prior[prior$mor != 0 & prior$source_genesymbol != prior$target_genesymbol, ]

nMR <- 8 # master regulators (layer 0)
nL1 <- 4 # TFs regulated by master regulators (layer 1)
nTargetsPerTF <- 25 # signed targets sampled per TF
maxRegs <- 3 # maximum regulators per target
nBins <- 10 # cell states (bins 0-4: cell type A, bins 5-9: cell type B)

# Candidate TFs: enough activating and repressing signed targets
tfStats <- do.call(rbind, lapply(split(signed$mor, signed$source_genesymbol), function(m) {
  data.frame(nPos = sum(m > 0), nNeg = sum(m < 0))
}))
tfStats$tf <- rownames(tfStats)
cand <- tfStats$tf[tfStats$nPos >= 40 & tfStats$nNeg >= 10]

# Layer-1 TFs need a signed incoming edge from a master regulator candidate
tfEdges <- signed[signed$source_genesymbol %in% cand & signed$target_genesymbol %in% cand, ]
l1Pool <- unique(tfEdges$target_genesymbol)
L1 <- sample(l1Pool, nL1)
mrPool <- setdiff(unique(tfEdges$source_genesymbol[tfEdges$target_genesymbol %in% L1]), L1)
MR <- sample(mrPool, min(nMR, length(mrPool)))
if (length(MR) < nMR) MR <- c(MR, sample(setdiff(cand, c(MR, L1)), nMR - length(MR)))
TFs <- c(MR, L1)

# MR -> L1 edges (keep at most 2 per L1 TF); guarantee each L1 TF has a regulator
l1Edges <- do.call(rbind, lapply(L1, function(t) {
  e <- tfEdges[tfEdges$target_genesymbol == t & tfEdges$source_genesymbol %in% MR, ]
  e[seq_len(min(2, nrow(e))), ]
}))
stopifnot(all(L1 %in% l1Edges$target_genesymbol))

# Targets: sample signed targets per TF (non-TF genes), then cap regulators per target
tgtEdges <- signed[signed$source_genesymbol %in% TFs & !signed$target_genesymbol %in% TFs, ]
picked <- unique(unlist(lapply(TFs, function(t) {
  g <- unique(tgtEdges$target_genesymbol[tgtEdges$source_genesymbol == t])
  sample(g, min(nTargetsPerTF, length(g)))
})))
tgtEdges <- tgtEdges[tgtEdges$target_genesymbol %in% picked, ]
tgtEdges <- do.call(rbind, lapply(split(tgtEdges, tgtEdges$target_genesymbol), function(e) {
  e[sample(nrow(e), min(maxRegs, nrow(e))), ]
}))

grn <- rbind(l1Edges, tgtEdges)[, c("source_genesymbol", "target_genesymbol", "mor")]
names(grn) <- c("tf", "target", "mor")
grn$K <- grn$mor * round(runif(nrow(grn), 1, 5), 2)
rownames(grn) <- NULL

genes <- c(MR, L1, sort(unique(tgtEdges$target_genesymbol)))
idx <- setNames(seq_along(genes) - 1L, genes)

# SERGIO targets file: target, nRegs, regIdx..., K...
rows <- lapply(split(grn, grn$target), function(e) {
  c(idx[e$target[1]], nrow(e), idx[e$tf], e$K)
})
writeLines(vapply(rows, function(r) paste(r, collapse = ","), ""), file.path(outDir, "targets.csv"))

# SERGIO master regulator file: MR index and production rate per bin
rates <- matrix(round(runif(length(MR) * nBins, 0.2, 4), 3), nrow = length(MR))
writeLines(vapply(seq_along(MR), function(i) paste(c(idx[MR[i]], rates[i, ]), collapse = ","), ""),
           file.path(outDir, "regs.csv"))

write.csv(data.frame(index = idx, gene = genes, role = ifelse(genes %in% MR, "MR", ifelse(genes %in% L1, "L1", "target"))),
          file.path(outDir, "genes.csv"), row.names = FALSE)
write.csv(grn, file.path(outDir, "grn.csv"), row.names = FALSE)

cat(sprintf("GRN: %d genes (%d MR, %d L1, %d targets), %d edges (%d activating, %d repressing)\n",
            length(genes), length(MR), length(L1), length(genes) - length(TFs),
            nrow(grn), sum(grn$mor > 0), sum(grn$mor < 0)))
