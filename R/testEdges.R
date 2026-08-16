#' @importFrom stats t.test p.adjust pf model.matrix
#' @importFrom cli cli_alert_info cli_alert_success cli_abort
#' @title Test edges from SCORPION networks
#' @description Performs statistical testing of network edges from runSCORPION output.
#' Supports single-sample tests (testing if edges differ from zero) and two-sample
#' tests (comparing edges between two groups).
#' @author Daniel Osorio <daniecos@uio.no>
#' @param networksDF A data.frame output from \code{\link{runSCORPION}} containing
#'   TF-target pairs as rows and network identifiers as columns.
#' @param testType Character specifying the test type. Options are:
#'   \itemize{
#'     \item{"single": Single-sample test (one-sample t-test against zero)}
#'     \item{"two.sample": Two-sample comparison (t-test between two groups)}
#'   }
#' @param group1 Character vector of column names in \code{networksDF} representing
#'   the first group (or the only group for single-sample tests).
#' @param group2 Character vector of column names in \code{networksDF} representing
#'   the second group. Required for two-sample tests, ignored for single-sample tests.
#' @param paired Logical indicating whether to perform a paired t-test. Default FALSE.
#'   When TRUE, group1 and group2 must have the same length and be in matched order
#'   (e.g., group1[1] is paired with group2[1]). Useful for comparing matched samples
#'   such as Tumor vs Normal from the same patient.
#' @param alternative Character specifying the alternative hypothesis. Options:
#'   "two.sided" (default), "greater", or "less".
#' @param padjustMethod Character specifying the p-value adjustment method for multiple
#'   testing correction. See \code{\link[stats]{p.adjust}} for options. Default "BH"
#'   (Benjamini-Hochberg FDR).
#' @param minLog2FC Numeric threshold for minimum absolute log2 fold change to
#'   include in testing. For two-sample and paired tests, edges with |log2FoldChange|
#'   below this threshold are excluded. Not applicable for single-sample tests.
#'   Default 0.
#' @param moderateVariance Logical indicating whether to apply SAM-style variance
#'   moderation. When TRUE, adds a fudge factor (s0, the median of all standard errors)
#'   to the denominator of the t-statistic. This prevents edges with very small
#'   variance from producing extreme t-statistics, resulting in volcano plots more
#'   similar to limma output. Default TRUE.
#' @param empiricalNull Logical indicating whether to estimate the null distribution
#'   empirically from the observed t-statistics. When TRUE, uses the median and MAD
#'   (median absolute deviation) of all t-statistics to recenter and rescale them,
#'   then computes p-values from the standard normal. This is Efron's empirical null
#'   correction (as in locfdr) and is essential when testing millions of correlated
#'   edges. Runs in O(n) time. Default TRUE.
#' @param nCores Integer specifying the number of parallel workers. Default 1
#'   (sequential processing). When greater than 1, edges are split into batches and
#'   processed in parallel using \code{furrr::future_map_dfr}. Requires the
#'   \pkg{furrr} and \pkg{future} packages to be installed.
#' @param batchSize Integer specifying the number of edges (rows) per batch for
#'   parallel processing. Default NULL, which auto-calculates as
#'   \code{ceiling(nrow(networksDF) / nCores)}. Only used when \code{nCores > 1}.
#'   Smaller batch sizes use less memory per worker but add communication overhead.
#' @return A data.frame containing:
#'   \itemize{
#'     \item{tf: Transcription factor}
#'     \item{target: Target gene}
#'     \item{meanEdge: Mean edge weight}
#'     \item{SE: Raw, unmoderated sampling standard error}
#'     \item{tStatistic: Test statistic}
#'     \item{pValue: Raw p-value}
#'     \item{pAdj: Adjusted p-value}
#'     \item{For two-sample tests: meanGroup1, meanGroup2, cohensD, log2FoldChange (Group1 - Group2)}
#'   }
#' @seealso \code{\link{runSCORPION}}, \code{\link{maEdges}}, \code{\link{circosEdges}}
#' @details
#' For single-sample tests, the function tests whether the mean edge weight across
#' replicates significantly differs from zero using a one-sample t-test.
#'
#' For two-sample tests, the function compares edge weights between two groups
#' using Welch's t-test (unequal variances assumed).
#'
#' For paired tests, the function calculates the difference between matched pairs
#' and performs a one-sample t-test on the differences (testing if mean difference
#' differs from zero). This is appropriate when samples are matched (e.g., Tumor
#' and Normal from the same patient).
#'
#' The returned \code{SE} is the raw sampling standard error before optional
#' SAM-style variance moderation. It is intended for downstream effect-size
#' meta-analysis. The moderated SE is used only internally for calculating
#' the test statistic and p-value.
#'
#' Edges are tested independently, and p-values are adjusted for multiple testing
#' using the specified method.
#'
#' The function uses fully vectorized computations for efficiency, making it suitable
#' for large-scale analyses with millions of edges. T-statistics and p-values are
#' calculated using matrix operations without iteration.
#' @examples
#' \dontrun{
#' # Load test data and build networks by donor and region
#' # Note: T = Tumor, N = Normal, B = Border regions
#' data(scorpionTest)
#' nets <- runSCORPION(
#'   gexMatrix = scorpionTest$gex,
#'   tfMotifs = scorpionTest$tf,
#'   ppiNet = scorpionTest$ppi,
#'   cellsMetadata = scorpionTest$metadata,
#'   groupBy = c("donor", "region")
#' )
#'
#' # Single-sample test: Test if edges in Tumor region differ from zero
#' tumor_nets <- grep("--T$", colnames(nets), value = TRUE)
#' results_single <- testEdges(
#'   networksDF = nets,
#'   testType = "single",
#'   group1 = tumor_nets
#' )
#'
#' # Two-sample test: Compare Tumor vs Border regions
#' tumor_nets <- grep("--T$", colnames(nets), value = TRUE)
#' border_nets <- grep("--B$", colnames(nets), value = TRUE)
#' results_tumor_vs_border <- testEdges(
#'   networksDF = nets,
#'   testType = "two.sample",
#'   group1 = tumor_nets,
#'   group2 = border_nets
#' )
#'
#' # View top differential edges (Tumor vs Border)
#' head(results_tumor_vs_border[order(results_tumor_vs_border$pAdj), ])
#'
#' # Compare Tumor vs Normal regions
#' normal_nets <- grep("--N$", colnames(nets), value = TRUE)
#' results_tumor_vs_normal <- testEdges(
#'   networksDF = nets,
#'   testType = "two.sample",
#'   group1 = tumor_nets,
#'   group2 = normal_nets
#' )
#'
#' # Filter by minimum log2 fold change for focused analysis
#' results_filtered <- testEdges(
#'   networksDF = nets,
#'   testType = "two.sample",
#'   group1 = tumor_nets,
#'   group2 = normal_nets,
#'   minLog2FC = 0.5
#' )
#'
#' # Paired t-test: Compare matched Tumor vs Normal samples
#' tumor_nets_ordered <- c("P31--T", "P32--T", "P33--T")
#' normal_nets_ordered <- c("P31--N", "P32--N", "P33--N")
#' results_paired <- testEdges(
#'   networksDF = nets,
#'   testType = "two.sample",
#'   group1 = tumor_nets_ordered,
#'   group2 = normal_nets_ordered,
#'   paired = TRUE
#' )
#' }
#' @export
#' @importFrom stats pt p.adjust var sd qnorm pnorm mad median
testEdges <- function(networksDF,
                      testType = c("single", "two.sample"),
                      group1,
                      group2 = NULL,
                      paired = FALSE,
                      alternative = c("two.sided", "greater", "less"),
                      padjustMethod = "BH",
                      minLog2FC = 0,
                      moderateVariance = TRUE,
                      empiricalNull = TRUE,
                      nCores = 1L,
                      batchSize = NULL) {
  
  testType <- match.arg(testType)
  alternative <- match.arg(alternative)
  
  if (missing(group1) || is.null(group1)) {
    cli::cli_abort("group1 must be specified")
  }
  
  if (!all(group1 %in% colnames(networksDF))) {
    missing_cols <- setdiff(group1, colnames(networksDF))
    cli::cli_abort(
      "Some group1 columns not found in networksDF: {paste(missing_cols, collapse = ', ')}"
    )
  }
  
  if (testType == "two.sample") {
    if (is.null(group2)) {
      cli::cli_abort("group2 must be specified for two.sample test")
    }
    if (!all(group2 %in% colnames(networksDF))) {
      missing_cols <- setdiff(group2, colnames(networksDF))
      cli::cli_abort(
        "Some group2 columns not found in networksDF: {paste(missing_cols, collapse = ', ')}"
      )
    }
    if (paired && length(group1) != length(group2)) {
      cli::cli_abort(
        "For paired tests, group1 and group2 must have the same length"
      )
    }
  }
  
  if (paired && testType == "single") {
    cli::cli_abort("Paired tests require testType = 'two.sample'")
  }
  
  nCores <- as.integer(nCores)
  if (length(nCores) != 1 || is.na(nCores) || nCores < 1L) {
    cli::cli_abort(
      "{.arg nCores} must be a positive integer, got {.val {nCores}}"
    )
  }
  
  if (!is.null(batchSize)) {
    batchSize <- as.integer(batchSize)
    if (length(batchSize) != 1 || is.na(batchSize) || batchSize < 1L) {
      cli::cli_abort(
        "{.arg batchSize} must be a positive integer or NULL, got {.val {batchSize}}"
      )
    }
  }
  
  if (nCores > 1L) {
    if (!requireNamespace("furrr", quietly = TRUE) ||
        !requireNamespace("future", quietly = TRUE)) {
      cli::cli_abort(c(
        "Packages {.pkg furrr} and {.pkg future} are required for parallel processing.",
        "i" = "Install them with {.code install.packages(c('furrr', 'future'))}",
        "i" = "Or set {.code nCores = 1} to run sequentially."
      ))
    }
  }
  
  tf_col <- which(colnames(networksDF) == "tf")
  target_col <- which(colnames(networksDF) == "target")
  
  if (length(tf_col) == 0 || length(target_col) == 0) {
    cli::cli_abort("networksDF must contain 'tf' and 'target' columns")
  }
  
  n_edges <- nrow(networksDF)
  
  cli::cli_alert_info(
    "Testing {n_edges} edges ({testType} test{if (paired) ', paired' else ''})"
  )
  
  s0 <- NULL
  if (moderateVariance && nCores > 1L) {
    s0 <- .computeGlobalS0(
      networksDF = networksDF,
      testType = testType,
      group1 = group1,
      group2 = group2,
      paired = paired
    )
  }
  
  if (testType == "single") {
    .helperFn <- testEdgesSingle
    .helperArgs <- list(
      group1 = group1,
      alternative = alternative,
      moderateVariance = moderateVariance,
      s0 = s0
    )
  } else if (paired) {
    .helperFn <- testEdgesPaired
    .helperArgs <- list(
      group1 = group1,
      group2 = group2,
      alternative = alternative,
      minLog2FC = minLog2FC,
      moderateVariance = moderateVariance,
      s0 = s0
    )
  } else {
    .helperFn <- testEdgesTwoSample
    .helperArgs <- list(
      group1 = group1,
      group2 = group2,
      alternative = alternative,
      minLog2FC = minLog2FC,
      moderateVariance = moderateVariance,
      s0 = s0
    )
  }
  
  if (nCores == 1L) {
    results <- do.call(
      .helperFn,
      c(list(networksDF = networksDF), .helperArgs)
    )
  } else {
    if (is.null(batchSize)) {
      batchSize <- ceiling(n_edges / nCores)
    }
    
    cli::cli_alert_info(
      "Using {nCores} workers, batch size {batchSize}"
    )
    
    batch_indices <- split(
      seq_len(n_edges),
      ceiling(seq_len(n_edges) / batchSize)
    )
    
    chunks <- lapply(
      batch_indices,
      function(idx) networksDF[idx, , drop = FALSE]
    )
    
    old_maxsize <- getOption("future.globals.maxSize")
    options(future.globals.maxSize = Inf)
    on.exit(
      options(future.globals.maxSize = old_maxsize),
      add = TRUE
    )
    
    old_plan <- future::plan(
      future::multisession,
      workers = nCores
    )
    
    on.exit({
      future::plan(future::sequential)
      future::plan(old_plan)
    }, add = TRUE)
    
    results <- furrr::future_map_dfr(
      chunks,
      function(chunk) {
        do.call(
          .helperFn,
          c(list(networksDF = chunk), .helperArgs)
        )
      },
      .options = furrr::furrr_options(seed = NULL)
    )
    
    future::plan(future::sequential)
    future::plan(old_plan)
    options(future.globals.maxSize = old_maxsize)
    on.exit()
  }
  
  # Apply empirical null correction
  if (empiricalNull) {
    t_stats <- results$tStatistic
    valid_idx <- !is.na(t_stats) & is.finite(t_stats)
    
    if (sum(valid_idx) > 100) {
      null_center <- median(t_stats[valid_idx])
      null_scale <- mad(
        t_stats[valid_idx],
        constant = 1.4826
      )
      
      if (null_scale > 0) {
        z_stats <- (t_stats - null_center) / null_scale
        
        results$pValue <- switch(
          alternative,
          "two.sided" = 2 * pnorm(
            abs(z_stats),
            lower.tail = FALSE
          ),
          "greater" = pnorm(
            z_stats,
            lower.tail = FALSE
          ),
          "less" = pnorm(
            z_stats,
            lower.tail = TRUE
          )
        )
      }
    }
  }
  
  results$pAdj <- p.adjust(
    results$pValue,
    method = padjustMethod
  )
  
  rownames(results) <- NULL
  
  cli::cli_alert_success(
    "Tested {nrow(results)} edges"
  )
  
  return(results)
}


# ============================================================
# Global SAM fudge factor
# ============================================================

.computeGlobalS0 <- function(
    networksDF,
    testType,
    group1,
    group2,
    paired) {
  
  if (testType == "single") {
    
    edge_matrix <- as.matrix(
      networksDF[, group1, drop = FALSE]
    )
    
    meanEdge <- rowMeans(
      edge_matrix,
      na.rm = TRUE
    )
    
    n_valid <- rowSums(
      !is.na(edge_matrix)
    )
    
    row_mean_sq <- rowMeans(
      edge_matrix^2,
      na.rm = TRUE
    )
    
    sd_edge <- sqrt(
      n_valid / (n_valid - 1) *
        (row_mean_sq - meanEdge^2)
    )
    
    se <- sd_edge / sqrt(n_valid)
    
  } else if (paired) {
    
    edge_matrix1 <- as.matrix(
      networksDF[, group1, drop = FALSE]
    )
    
    edge_matrix2 <- as.matrix(
      networksDF[, group2, drop = FALSE]
    )
    
    diff_matrix <- edge_matrix1 - edge_matrix2
    
    diffMean <- rowMeans(
      diff_matrix,
      na.rm = TRUE
    )
    
    valid_pairs <- !is.na(edge_matrix1) &
      !is.na(edge_matrix2)
    
    n_valid <- rowSums(valid_pairs)
    
    diff_mean_sq <- rowMeans(
      diff_matrix^2,
      na.rm = TRUE
    )
    
    sd_diff <- sqrt(
      n_valid / (n_valid - 1) *
        (diff_mean_sq - diffMean^2)
    )
    
    se <- sd_diff / sqrt(n_valid)
    
  } else {
    
    edge_matrix1 <- as.matrix(
      networksDF[, group1, drop = FALSE]
    )
    
    edge_matrix2 <- as.matrix(
      networksDF[, group2, drop = FALSE]
    )
    
    meanEdge1 <- rowMeans(
      edge_matrix1,
      na.rm = TRUE
    )
    
    meanEdge2 <- rowMeans(
      edge_matrix2,
      na.rm = TRUE
    )
    
    n1 <- rowSums(!is.na(edge_matrix1))
    n2 <- rowSums(!is.na(edge_matrix2))
    
    row_mean_sq1 <- rowMeans(
      edge_matrix1^2,
      na.rm = TRUE
    )
    
    row_mean_sq2 <- rowMeans(
      edge_matrix2^2,
      na.rm = TRUE
    )
    
    var1 <- n1 / (n1 - 1) *
      (row_mean_sq1 - meanEdge1^2)
    
    var2 <- n2 / (n2 - 1) *
      (row_mean_sq2 - meanEdge2^2)
    
    se <- sqrt(
      var1 / n1 +
        var2 / n2
    )
  }
  
  median(
    se,
    na.rm = TRUE
  )
}


# ============================================================
# Single-sample test
# ============================================================

testEdgesSingle <- function(
    networksDF,
    group1,
    alternative,
    moderateVariance = TRUE,
    s0 = NULL) {
  
  edge_data <- networksDF[
    ,
    group1,
    drop = FALSE
  ]
  
  edge_matrix <- as.matrix(edge_data)
  
  meanEdge <- rowMeans(
    edge_matrix,
    na.rm = TRUE
  )
  
  tf_target <- networksDF[
    ,
    c("tf", "target")
  ]
  
  n_samples <- ncol(edge_matrix)
  
  n_valid <- rowSums(
    !is.na(edge_matrix)
  )
  
  row_mean_sq <- rowMeans(
    edge_matrix^2,
    na.rm = TRUE
  )
  
  sd_edge <- sqrt(
    n_valid / (n_valid - 1) *
      (row_mean_sq - meanEdge^2)
  )
  
  # Raw/unmoderated SE for downstream meta-analysis
  rawSE <- sd_edge / sqrt(n_valid)
  
  # SE used for hypothesis testing
  se <- rawSE
  
  if (moderateVariance) {
    if (is.null(s0)) {
      s0 <- median(
        rawSE,
        na.rm = TRUE
      )
    }
    
    se <- rawSE + s0
  }
  
  test_stats <- meanEdge / se
  
  df <- n_valid - 1
  
  pvalues <- switch(
    alternative,
    "two.sided" = 2 * pt(
      abs(test_stats),
      df = df,
      lower.tail = FALSE
    ),
    "greater" = pt(
      test_stats,
      df = df,
      lower.tail = FALSE
    ),
    "less" = pt(
      test_stats,
      df = df,
      lower.tail = TRUE
    )
  )
  
  insufficient_data <-
    n_valid < 2 |
    is.na(sd_edge) |
    (!moderateVariance & sd_edge == 0)
  
  test_stats[insufficient_data] <- NA
  pvalues[insufficient_data] <- NA
  rawSE[insufficient_data] <- NA
  
  results <- data.frame(
    tf = tf_target$tf,
    target = tf_target$target,
    meanEdge = meanEdge,
    SE = rawSE,
    tStatistic = test_stats,
    pValue = pvalues,
    stringsAsFactors = FALSE
  )
  
  return(results)
}


# ============================================================
# Two-sample test
# ============================================================

testEdgesTwoSample <- function(
    networksDF,
    group1,
    group2,
    alternative,
    minLog2FC,
    moderateVariance = TRUE,
    s0 = NULL) {
  
  edge_data1 <- networksDF[
    ,
    group1,
    drop = FALSE
  ]
  
  edge_data2 <- networksDF[
    ,
    group2,
    drop = FALSE
  ]
  
  edge_matrix1 <- as.matrix(edge_data1)
  edge_matrix2 <- as.matrix(edge_data2)
  
  meanEdge1 <- rowMeans(
    edge_matrix1,
    na.rm = TRUE
  )
  
  meanEdge2 <- rowMeans(
    edge_matrix2,
    na.rm = TRUE
  )
  
  meanEdge <- (
    meanEdge1 + meanEdge2
  ) / 2
  
  diffMean <- meanEdge1 - meanEdge2
  
  log2FC <- meanEdge1 - meanEdge2
  
  keep_idx <- abs(log2FC) >= minLog2FC
  
  edge_matrix1 <- edge_matrix1[
    keep_idx,
    ,
    drop = FALSE
  ]
  
  edge_matrix2 <- edge_matrix2[
    keep_idx,
    ,
    drop = FALSE
  ]
  
  meanEdge1 <- meanEdge1[keep_idx]
  meanEdge2 <- meanEdge2[keep_idx]
  meanEdge <- meanEdge[keep_idx]
  diffMean <- diffMean[keep_idx]
  log2FC <- log2FC[keep_idx]
  
  tf_target <- networksDF[
    keep_idx,
    c("tf", "target")
  ]
  
  n1 <- rowSums(
    !is.na(edge_matrix1)
  )
  
  n2 <- rowSums(
    !is.na(edge_matrix2)
  )
  
  row_mean_sq1 <- rowMeans(
    edge_matrix1^2,
    na.rm = TRUE
  )
  
  row_mean_sq2 <- rowMeans(
    edge_matrix2^2,
    na.rm = TRUE
  )
  
  var1 <- n1 / (n1 - 1) *
    (row_mean_sq1 - meanEdge1^2)
  
  var2 <- n2 / (n2 - 1) *
    (row_mean_sq2 - meanEdge2^2)
  
  # Raw Welch SE
  rawSE <- sqrt(
    var1 / n1 +
      var2 / n2
  )
  
  # SE used for hypothesis testing
  se <- rawSE
  
  if (moderateVariance) {
    if (is.null(s0)) {
      s0 <- median(
        rawSE,
        na.rm = TRUE
      )
    }
    
    se <- rawSE + s0
  }
  
  test_stats <- diffMean / se
  
  df <- (
    var1 / n1 +
      var2 / n2
  )^2 / (
    (var1 / n1)^2 / (n1 - 1) +
      (var2 / n2)^2 / (n2 - 1)
  )
  
  pvalues <- switch(
    alternative,
    "two.sided" = 2 * pt(
      abs(test_stats),
      df = df,
      lower.tail = FALSE
    ),
    "greater" = pt(
      test_stats,
      df = df,
      lower.tail = FALSE
    ),
    "less" = pt(
      test_stats,
      df = df,
      lower.tail = TRUE
    )
  )
  
  se_before_mod <- rawSE
  
  insufficient_data <-
    n1 < 2 |
    n2 < 2 |
    is.na(var1) |
    is.na(var2) |
    (
      !moderateVariance &
        (
          var1 == 0 |
            var2 == 0 |
            se_before_mod == 0
        )
    ) |
    is.na(rawSE)
  
  test_stats[insufficient_data] <- NA
  pvalues[insufficient_data] <- NA
  df[insufficient_data] <- NA
  rawSE[insufficient_data] <- NA
  
  pooled_var <- (
    (n1 - 1) * var1 +
      (n2 - 1) * var2
  ) / (
    n1 + n2 - 2
  )
  
  pooled_sd <- sqrt(pooled_var)
  
  cohensD <- diffMean / pooled_sd
  
  cohensD[
    pooled_sd == 0 |
      is.na(pooled_sd)
  ] <- NA
  
  results <- data.frame(
    tf = tf_target$tf,
    target = tf_target$target,
    meanGroup1 = meanEdge1,
    meanGroup2 = meanEdge2,
    cohensD = cohensD,
    log2FoldChange = log2FC,
    meanEdge = meanEdge,
    SE = rawSE,
    tStatistic = test_stats,
    pValue = pvalues,
    stringsAsFactors = FALSE
  )
  
  return(results)
}


# ============================================================
# Paired two-sample test
# ============================================================

testEdgesPaired <- function(
    networksDF,
    group1,
    group2,
    alternative,
    minLog2FC,
    moderateVariance = TRUE,
    s0 = NULL) {
  
  edge_data1 <- networksDF[
    ,
    group1,
    drop = FALSE
  ]
  
  edge_data2 <- networksDF[
    ,
    group2,
    drop = FALSE
  ]
  
  edge_matrix1 <- as.matrix(edge_data1)
  edge_matrix2 <- as.matrix(edge_data2)
  
  meanEdge1 <- rowMeans(
    edge_matrix1,
    na.rm = TRUE
  )
  
  meanEdge2 <- rowMeans(
    edge_matrix2,
    na.rm = TRUE
  )
  
  meanEdge <- (
    meanEdge1 + meanEdge2
  ) / 2
  
  diff_matrix <- edge_matrix1 - edge_matrix2
  
  diffMean <- rowMeans(
    diff_matrix,
    na.rm = TRUE
  )
  
  log2FC <- meanEdge1 - meanEdge2
  
  keep_idx <- abs(log2FC) >= minLog2FC
  
  diff_matrix <- diff_matrix[
    keep_idx,
    ,
    drop = FALSE
  ]
  
  edge_matrix1 <- edge_matrix1[
    keep_idx,
    ,
    drop = FALSE
  ]
  
  edge_matrix2 <- edge_matrix2[
    keep_idx,
    ,
    drop = FALSE
  ]
  
  meanEdge1 <- meanEdge1[keep_idx]
  meanEdge2 <- meanEdge2[keep_idx]
  meanEdge <- meanEdge[keep_idx]
  diffMean <- diffMean[keep_idx]
  log2FC <- log2FC[keep_idx]
  
  tf_target <- networksDF[
    keep_idx,
    c("tf", "target")
  ]
  
  valid_pairs <- !is.na(edge_matrix1) &
    !is.na(edge_matrix2)
  
  n_valid <- rowSums(valid_pairs)
  
  diff_mean_sq <- rowMeans(
    diff_matrix^2,
    na.rm = TRUE
  )
  
  sd_diff <- sqrt(
    n_valid / (n_valid - 1) *
      (diff_mean_sq - diffMean^2)
  )
  
  # FIX:
  # Preserve the raw paired SE before variance moderation.
  rawSE <- sd_diff / sqrt(n_valid)
  
  # SE used for hypothesis testing
  se <- rawSE
  
  if (moderateVariance) {
    if (is.null(s0)) {
      s0 <- median(
        rawSE,
        na.rm = TRUE
      )
    }
    
    se <- rawSE + s0
  }
  
  test_stats <- diffMean / se
  
  df <- n_valid - 1
  
  pvalues <- switch(
    alternative,
    "two.sided" = 2 * pt(
      abs(test_stats),
      df = df,
      lower.tail = FALSE
    ),
    "greater" = pt(
      test_stats,
      df = df,
      lower.tail = FALSE
    ),
    "less" = pt(
      test_stats,
      df = df,
      lower.tail = TRUE
    )
  )
  
  insufficient_data <-
    n_valid < 2 |
    is.na(sd_diff) |
    (!moderateVariance & sd_diff == 0)
  
  test_stats[insufficient_data] <- NA
  pvalues[insufficient_data] <- NA
  rawSE[insufficient_data] <- NA
  
  cohensD <- diffMean / sd_diff
  
  cohensD[
    sd_diff == 0 |
      is.na(sd_diff)
  ] <- NA
  
  results <- data.frame(
    tf = tf_target$tf,
    target = tf_target$target,
    meanGroup1 = meanEdge1,
    meanGroup2 = meanEdge2,
    cohensD = cohensD,
    log2FoldChange = log2FC,
    meanEdge = meanEdge,
    SE = rawSE,
    tStatistic = test_stats,
    pValue = pvalues,
    stringsAsFactors = FALSE
  )
  
  return(results)
}
