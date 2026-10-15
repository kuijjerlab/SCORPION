# Shared helpers for inferEdgeSign() tests on the SERGIO ground-truth simulation
# (fixtures/sergio_signed.rds, built by data-raw/sergio/).

load_sergio <- function() {
  readRDS(test_path("fixtures", "sergio_signed.rds"))
}

sergio_prior <- function(fx) {
  fx$prior[, c("source_genesymbol", "target_genesymbol", "weight", "mor")]
}

run_sign <- function(fx, prior, gex = fx$gex, ...) {
  inferEdgeSign(gex, prior, fx$metadata, groupBy = "cell_type", sampleBy = "donor",
                gammaValue = 5, nBoot = 30, showProgress = FALSE, ...)$edges
}

# Edges whose sign was called from the data (with or against the prior)
data_called <- function(e) e$evidence %in% c("data", "data+prior", "data_vs_prior")

# Sign supported by the data: conflicting edges carry the opposite of the prior sign
data_sign <- function(e) {
  ifelse(e$label == "activating", 1,
         ifelse(e$label == "repressing", -1,
                ifelse(e$label == "conflicting", -e$priorSign, 0)))
}

true_sign <- function(fx, e) {
  fx$grn$mor[match(paste(e$tf, e$target), paste(fx$grn$tf, fx$grn$target))]
}

# Within-donor permutation of each row of a cells matrix (destroys gene-gene dependence)
permute_within_donor <- function(m, donor) {
  for (d in unique(donor)) {
    i <- which(donor == d)
    for (r in seq_len(nrow(m))) m[r, i] <- m[r, i[sample.int(length(i))]]
  }
  m
}
