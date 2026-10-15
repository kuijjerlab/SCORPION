#' @importFrom stats p.adjust pt model.matrix var
#' @importFrom cli cli_abort cli_alert_info cli_alert_success cli_alert_warning
#' @export inferEdgeSign
#' @title Infer the regulatory sign (activation or repression) of TF-target edges
#' @description
#' Estimates whether each transcription factor (TF) activates or represses each of
#' its prior targets from single-cell expression data. Signed, partially signed and
#' unsigned priors are all accepted and produce the same output; prior signs only
#' improve accuracy. TF activity is estimated internally from the expression of the
#' TF's targets, never from the TF's own mRNA (which, under SCORPION, reflects a future
#' rather than the current state of the cell). When a TF has enough signed targets, its
#' activity is their signed weighted sum; otherwise it is the leading axis of all its
#' targets, cross-fitted over 5 folds of the observations (each fold's activity comes
#' from axes learned on the other folds), with its direction set from the expected signs
#' of the TF's other targets (their prior sign when known, otherwise the prevailing sign
#' of the prior).
#' In every case the target being tested is left out of the activity it is tested
#' against, so neither its expression nor its own prior sign determines its estimate.
#' Within each group, every target is regressed jointly on the activities of all its
#' prior regulators after removing donor effects and optional covariates. Signs are
#' called only when they are significant after multiple-testing correction and stable
#' under bootstrap resampling; overlapping regulons that cannot be separated are
#' reported as unresolved. The original SCORPION networks are not modified.
#'
#' A sign call means the target varies consistently with activation (or repression)
#' by the TF in the data. It does not show that the TF binds the target directly.
#' @author Daniel Osorio <daniecos@uio.no>
#' @param gexMatrix An expression dataset with genes in the rows and barcodes (cells) in the columns.
#' @param tfMotifs A TF-target prior (data.frame or matrix) with columns TF, target,
#'   weight and, optionally, \code{mor} (mode of regulation: 1 activation, -1 repression,
#'   0 unknown). If no \code{mor} column is present, the fourth column is used; if there
#'   are only three columns, all signs are treated as unknown.
#' @param cellsMetadata A data.frame with cell-level metadata, one row per column of \code{gexMatrix}.
#' @param groupBy Character vector of column names in \code{cellsMetadata} defining the groups
#'   (e.g. cell types) analyzed separately. Group identifiers are joined with "--" as in
#'   \code{\link{runSCORPION}}.
#' @param sampleBy Optional column name in \code{cellsMetadata} identifying donors or samples.
#'   Super-cells are built within each sample, sample effects are removed with fixed effects,
#'   and the bootstrap resamples within (fewer than 5 samples) or across (5 or more) samples.
#' @param covariates Optional character vector of numeric columns in \code{cellsMetadata}
#'   (e.g. cell-cycle scores) to adjust for.
#' @param edges Optional data.frame with columns \code{tf} and \code{target} restricting
#'   which prior edges are tested. All prior regulators of a target are still used in its model.
#' @param scorpionEdges Optional output of \code{\link{runSCORPION}}. When a column matches a
#'   group identifier, its edge weights are returned in the \code{score} column. Edges selected
#'   on these scores should come from data independent of the data used here.
#' @param normalizeData Boolean to indicate log normalization of the expression data. Default TRUE.
#' @param gammaValue Graining level (cells per super-cell) used within each sample. Use 1 to
#'   analyze single cells. Default 10.
#' @param nPC Number of principal components used to build super-cells. Default 25.
#' @param minCells Minimum number of cells required for a group. Default 30.
#' @param minTargets Minimum number of targets (more than this) a TF needs in the data
#'   to be tested. Default 5.
#' @param vifThreshold Variance inflation factor above which a regulator's effect is considered
#'   not identifiable from co-regulators. Default 10.
#' @param nBoot Number of bootstrap resamples used to assess sign stability. Default 100.
#' @param stabilityThreshold Minimum fraction of bootstrap resamples that must agree with the
#'   estimated sign. Default 0.9.
#' @param alpha Significance threshold on adjusted p-values. Default 0.05.
#' @param padjustMethod Multiple-testing correction method, see \code{\link[stats]{p.adjust}}.
#'   Applied within each group. Default "BH".
#' @param seed Random seed for super-cell construction and the bootstrap. Default 1.
#' @param showProgress Boolean to indicate printing of progress messages. Default TRUE.
#' @return A list with:
#'   \itemize{
#'     \item{edges: data.frame with one row per tested prior edge and group: \code{group},
#'       \code{tf}, \code{target}, \code{score} (SCORPION edge weight, if supplied),
#'       \code{priorSign}, \code{estimate} (positive: activation), \code{SE},
#'       \code{tStatistic}, \code{df}, \code{pValue}, \code{pAdj}, \code{signStability},
#'       \code{vif}, \code{label} ("activating", "repressing", "conflicting", "unresolved"
#'       or "not_testable"), \code{evidence} ("data", "data+prior", "data_vs_prior",
#'       "prior_only" or "none")
#'       and \code{reason}, why the sign was not called from the data ("self_edge",
#'       "target_not_expressed", "tf_activity_unavailable", "insufficient_df", "collinear",
#'       "not_significant" or "unstable"; NA when it was called).}
#'     \item{activity: list with one TF x observation activity matrix per group (super-cells,
#'       or cells when \code{gammaValue = 1}, in columns; their sample in the "sample"
#'       attribute).}
#'     \item{params: list of the parameters used.}
#'   }
#' @seealso \code{\link{runSCORPION}}, \code{\link{maEdges}}
#' @examples
#' \donttest{
#' data(scorpionTest)
#' signs <- inferEdgeSign(
#'   gexMatrix = scorpionTest$gex,
#'   tfMotifs = scorpionTest$tf,
#'   cellsMetadata = scorpionTest$metadata,
#'   groupBy = "region",
#'   sampleBy = "donor",
#'   nBoot = 20
#' )
#' table(signs$edges$label)
#' }
inferEdgeSign <- function(gexMatrix,
                          tfMotifs,
                          cellsMetadata,
                          groupBy,
                          sampleBy = NULL,
                          covariates = NULL,
                          edges = NULL,
                          scorpionEdges = NULL,
                          normalizeData = TRUE,
                          gammaValue = 10,
                          nPC = 25,
                          minCells = 30,
                          minTargets = 5,
                          vifThreshold = 10,
                          nBoot = 100,
                          stabilityThreshold = 0.9,
                          alpha = 0.05,
                          padjustMethod = "BH",
                          seed = 1,
                          showProgress = TRUE) {
  # Input validation
  if (is.null(rownames(gexMatrix)) || is.null(colnames(gexMatrix))) {
    cli::cli_abort("gexMatrix must have gene names as rownames and cell barcodes as colnames")
  }
  cellsMetadata <- as.data.frame(cellsMetadata)
  if (nrow(cellsMetadata) != ncol(gexMatrix)) {
    cli::cli_abort("cellsMetadata must have one row per column of gexMatrix")
  }
  needed <- c(groupBy, sampleBy, covariates)
  if (!all(needed %in% colnames(cellsMetadata))) {
    cli::cli_abort("Columns not found in cellsMetadata: {paste(setdiff(needed, colnames(cellsMetadata)), collapse = ', ')}")
  }
  if (!is.null(covariates) && !all(vapply(cellsMetadata[covariates], is.numeric, logical(1)))) {
    cli::cli_abort("covariates must be numeric columns of cellsMetadata")
  }
  if (!is.numeric(gammaValue) || gammaValue < 1) {
    cli::cli_abort("gammaValue must be a number >= 1")
  }
  if (!is.numeric(stabilityThreshold) || stabilityThreshold < 0.5 || stabilityThreshold > 1) {
    cli::cli_abort("stabilityThreshold must be between 0.5 and 1")
  }
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) {
    cli::cli_abort("alpha must be between 0 and 1")
  }
  if (!is.numeric(nBoot) || nBoot < 1) {
    cli::cli_abort("nBoot must be a positive integer")
  }
  if (!is.null(edges) && !all(c("tf", "target") %in% colnames(edges))) {
    cli::cli_abort("edges must have columns 'tf' and 'target'")
  }

  prior <- .es_prior(tfMotifs)
  prior <- prior[prior$target %in% rownames(gexMatrix) | prior$tf %in% rownames(gexMatrix), ]
  if (nrow(prior) == 0) {
    cli::cli_abort("No prior targets or TFs are present in gexMatrix")
  }

  if (isTRUE(normalizeData)) {
    gexMatrix <- log_normalize_data(gexMatrix)
  }
  group <- do.call(paste, c(cellsMetadata[, groupBy, drop = FALSE], sep = "--"))
  donor <- if (is.null(sampleBy)) rep("all", ncol(gexMatrix)) else as.character(cellsMetadata[[sampleBy]])
  covMat <- if (is.null(covariates)) NULL else as.matrix(cellsMetadata[, covariates, drop = FALSE])

  groups <- names(which(table(group) >= minCells))
  if (length(groups) == 0) {
    cli::cli_abort("No groups have enough cells (>= {minCells})")
  }
  if (showProgress) {
    cli::cli_alert_info("Inferring edge signs in {length(groups)} group{?s}")
  }

  set.seed(seed)
  edgeOut <- list()
  activityOut <- list()
  for (g in groups) {
    idx <- which(group == g)
    agg <- .es_aggregate(gexMatrix[, idx, drop = FALSE], donor[idx],
                         if (is.null(covMat)) NULL else covMat[idx, , drop = FALSE],
                         gammaValue, nPC, seed)
    if (is.null(agg)) {
      cli::cli_alert_warning("Group {g}: not enough cells per sample to build super-cells; skipped")
      next
    }
    std <- .es_standardize(agg$Y, agg$donor, agg$cov)
    genes <- rownames(std$Z)
    st <- .es_structure(prior, genes, std$ok, minTargets)

    # Edges to report: prior edges (optionally restricted by `edges`)
    ge <- prior
    if (!is.null(edges)) {
      ge <- ge[paste(ge$tf, ge$target) %in% paste(edges$tf, edges$target), ]
    }
    testableTF <- ge$tf %in% st$tfs
    targetOK <- ge$target %in% genes[std$ok]
    selfEdge <- ge$tf == ge$target
    ge$reason <- ifelse(selfEdge, "self_edge",
      ifelse(!targetOK, "target_not_expressed",
        ifelse(!testableTF, "tf_activity_unavailable", NA_character_)
      )
    )

    # Model for each target: all testable prior regulators (not just the edges reported)
    tested <- ge[is.na(ge$reason), ]
    fit <- NULL
    if (nrow(tested) > 0) {
      allRegs <- prior[prior$tf %in% st$tfs & prior$target %in% unique(tested$target) & prior$tf != prior$target, ]
      regsByTarget <- lapply(split(allRegs$tf, allRegs$target), function(x) match(sort(unique(x)), st$tfs))
      key <- paste(tested$tf, tested$target)
      fit <- .es_estimate(agg$Y, agg$donor, agg$cov, st, regsByTarget, key)

      # Bootstrap sign stability
      agree <- numeric(nrow(fit))
      for (b in seq_len(nBoot)) {
        rs <- .es_resample(agg$donor)
        fb <- .es_estimate(agg$Y[, rs$idx, drop = FALSE], rs$donor,
                           if (is.null(agg$cov)) NULL else agg$cov[rs$idx, , drop = FALSE],
                           st, regsByTarget, key)
        eb <- fb$estimate
        agree <- agree + (!is.na(eb) & sign(eb) == sign(fit$estimate) & eb != 0)
      }
      fit$signStability <- agree / nBoot
      fit$tStatistic <- fit$estimate / fit$SE
      fit$pValue <- 2 * stats::pt(-abs(fit$tStatistic), df = pmax(fit$df, 1))
      fit$pValue[!is.finite(fit$pValue)] <- NA
    }

    res <- data.frame(group = g, ge[, c("tf", "target")], score = NA_real_,
                      priorSign = ge$mor, stringsAsFactors = FALSE)
    m <- if (is.null(fit)) rep(NA_integer_, nrow(res)) else match(paste(res$tf, res$target), paste(fit$tf, fit$target))
    for (col in c("estimate", "SE", "tStatistic", "df", "pValue", "signStability", "vif")) {
      res[[col]] <- if (is.null(fit)) NA_real_ else fit[[col]][m]
    }
    res$pAdj <- stats::p.adjust(res$pValue, method = padjustMethod)
    res$reason <- ge$reason

    # Labels
    called <- !is.na(res$pAdj) & res$pAdj < alpha & res$signStability >= stabilityThreshold &
      !is.na(res$vif) & res$vif <= vifThreshold
    dataSign <- ifelse(called, sign(res$estimate), 0)
    res$reason[is.na(res$reason) & !is.na(res$df) & res$df < 3] <- "insufficient_df"
    res$reason[is.na(res$reason) & (is.na(res$vif) | res$vif > vifThreshold)] <- "collinear"
    res$reason[is.na(res$reason) & (is.na(res$pAdj) | res$pAdj >= alpha)] <- "not_significant"
    res$reason[is.na(res$reason) & !called] <- "unstable"
    notTestable <- res$reason %in% c("self_edge", "target_not_expressed", "tf_activity_unavailable")
    res$label <- ifelse(dataSign > 0, "activating", ifelse(dataSign < 0, "repressing", "unresolved"))
    res$evidence <- ifelse(called, ifelse(res$priorSign == 0, "data", "data+prior"), "none")
    conflict <- called & res$priorSign != 0 & res$priorSign != dataSign
    res$label[conflict] <- "conflicting"
    res$evidence[conflict] <- "data_vs_prior"
    priorOnly <- !called & res$priorSign != 0
    res$label[priorOnly] <- ifelse(res$priorSign[priorOnly] > 0, "activating", "repressing")
    res$evidence[priorOnly] <- "prior_only"
    res$label[notTestable & !priorOnly] <- "not_testable"

    if (!is.null(scorpionEdges) && g %in% colnames(scorpionEdges)) {
      res$score <- scorpionEdges[[g]][match(paste(res$tf, res$target), paste(scorpionEdges$tf, scorpionEdges$target))]
    }
    edgeOut[[g]] <- res
    stFull <- .es_axes(st, std$Z, orient = TRUE)
    act <- .es_activity(std$Z, stFull) / sqrt(stFull$W)
    attr(act, "sample") <- agg$donor
    activityOut[[g]] <- act
    if (showProgress) {
      cli::cli_alert_success("{g}: {sum(called)} sign{?s} called from data, {sum(!is.na(res$estimate))} edge{?s} tested")
    }
  }

  edgesDF <- do.call(rbind, edgeOut)
  rownames(edgesDF) <- NULL
  list(
    edges = edgesDF,
    activity = activityOut,
    params = list(groupBy = groupBy, sampleBy = sampleBy, covariates = covariates,
                  gammaValue = gammaValue, minTargets = minTargets,
                  vifThreshold = vifThreshold, nBoot = nBoot, stabilityThreshold = stabilityThreshold,
                  alpha = alpha, padjustMethod = padjustMethod, seed = seed)
  )
}
