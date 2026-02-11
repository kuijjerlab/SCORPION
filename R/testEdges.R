#' @importFrom stats t.test p.adjust pf model.matrix
#' @title Test edges from SCORPION networks
#' @description Performs statistical testing of network edges from runSCORPION output.
#' Supports single-sample tests (testing if edges differ from zero) and two-sample
#' tests (comparing edges between two groups).
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
#'   Default 1e-16.
#' @param moderateVariance Logical indicating whether to apply SAM-style variance
#'   moderation. When TRUE, adds a fudge factor (s0, the median of all standard errors)
#'   to the denominator of the t-statistic. This prevents edges with very small variance
#'   from producing extreme t-statistics, resulting in volcano plots more similar to
#'   limma output. Default TRUE.
#' @param empiricalNull Logical indicating whether to estimate the null distribution
#'   empirically from the observed t-statistics. When TRUE, uses the median and MAD
#'   (median absolute deviation) of all t-statistics to recenter and rescale them,
#'   then computes p-values from the standard normal. This is Efron's empirical null
#'   correction (as in locfdr) and is essential when testing millions of correlated
#'   edges. Runs in O(n) time. Default TRUE.
#' @return A data.frame containing:
#'   \itemize{
#'     \item{tf: Transcription factor}
#'     \item{target: Target gene}
#'     \item{meanEdge: Mean edge weight}
#'     \item{tStatistic: Test statistic}
#'     \item{pValue: Raw p-value}
#'     \item{pAdj: Adjusted p-value}
#'     \item{For two-sample tests: meanGroup1, meanGroup2, diffMean (Group1 - Group2), cohensD, log2FoldChange}
#'   }
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
#' tumor_nets <- grep("--T$", colnames(nets), value = TRUE)  # T = Tumor
#' results_single <- testEdges(
#'   networksDF = nets,
#'   testType = "single",
#'   group1 = tumor_nets
#' )
#'
#' # Two-sample test: Compare Tumor vs Border regions
#' tumor_nets <- grep("--T$", colnames(nets), value = TRUE)  # T = Tumor
#' border_nets <- grep("--B$", colnames(nets), value = TRUE)  # B = Border
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
#' normal_nets <- grep("--N$", colnames(nets), value = TRUE)  # N = Normal
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
#'   minLog2FC = 0.5  # Only test edges with |log2FoldChange| >= 0.5
#' )
#'
#' # Paired t-test: Compare matched Tumor vs Normal samples (same patient)
#' # Ensure columns are ordered by patient: P31--T with P31--N, P32--T with P32--N, etc.
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
                      minLog2FC = 1e-16,
                      moderateVariance = TRUE,
                      empiricalNull = TRUE) {
  
  # Input validation
  testType <- match.arg(testType)
  alternative <- match.arg(alternative)
  
  if (missing(group1) || is.null(group1)) {
    stop("group1 must be specified")
  }
  
  if (!all(group1 %in% colnames(networksDF))) {
    missing_cols <- setdiff(group1, colnames(networksDF))
    stop("Some group1 columns not found in networksDF: ", paste(missing_cols, collapse = ", "))
  }
  
  if (testType == "two.sample") {
    if (is.null(group2)) {
      stop("group2 must be specified for two.sample test")
    }
    if (!all(group2 %in% colnames(networksDF))) {
      missing_cols <- setdiff(group2, colnames(networksDF))
      stop("Some group2 columns not found in networksDF: ", paste(missing_cols, collapse = ", "))
    }
    if (paired && length(group1) != length(group2)) {
      stop("For paired tests, group1 and group2 must have the same length")
    }
  }
  
  if (paired && testType == "single") {
    stop("Paired tests require testType = 'two.sample'")
  }
  
  # Extract TF and target columns
  tf_col <- which(colnames(networksDF) == "tf")
  target_col <- which(colnames(networksDF) == "target")
  
  if (length(tf_col) == 0 || length(target_col) == 0) {
    stop("networksDF must contain 'tf' and 'target' columns")
  }
  
  # Perform tests based on testType
  if (testType == "single") {
    results <- testEdgesSingle(
      networksDF = networksDF,
      group1 = group1,
      alternative = alternative,
      moderateVariance = moderateVariance
    )
  } else if (paired) {
    results <- testEdgesPaired(
      networksDF = networksDF,
      group1 = group1,
      group2 = group2,
      alternative = alternative,
      minLog2FC = minLog2FC,
      moderateVariance = moderateVariance
    )
  } else {
    results <- testEdgesTwoSample(
      networksDF = networksDF,
      group1 = group1,
      group2 = group2,
      alternative = alternative,
      minLog2FC = minLog2FC,
      moderateVariance = moderateVariance
    )
  }
  
  # Apply empirical null correction (Efron's method) - O(n) time
  if (empiricalNull) {
    t_stats <- results$tStatistic
    valid_idx <- !is.na(t_stats) & is.finite(t_stats)
    
    if (sum(valid_idx) > 100) {
      # Estimate null distribution using median and MAD (robust to outliers)
      null_center <- median(t_stats[valid_idx])
      null_scale <- mad(t_stats[valid_idx], constant = 1.4826)  # scaled to match SD for normal
      
      if (null_scale > 0) {
        # Standardize t-statistics to empirical null
        z_stats <- (t_stats - null_center) / null_scale
        
        # Recompute p-values using standard normal
        results$pValue <- switch(alternative,
          "two.sided" = 2 * pnorm(abs(z_stats), lower.tail = FALSE),
          "greater" = pnorm(z_stats, lower.tail = FALSE),
          "less" = pnorm(z_stats, lower.tail = TRUE)
        )
      }
    }
  }
  
  # Adjust p-values
  results$pAdj <- p.adjust(results$pValue, method = padjustMethod)
  
  rownames(results) <- NULL
  
  return(results)
}

#' @keywords internal
testEdgesSingle <- function(networksDF, group1, alternative, 
                            moderateVariance = TRUE) {
  
  # Extract edge weights for group1
  edge_data <- networksDF[, group1, drop = FALSE]
  edge_matrix <- as.matrix(edge_data)
  
  # Calculate mean edge weight
  meanEdge <- rowMeans(edge_matrix, na.rm = TRUE)
  
  # No filtering for single-sample tests (minLog2FC not applicable)
  tf_target <- networksDF[, c("tf", "target")]
  
  # Vectorized one-sample t-test computation
  # Calculate statistics for all edges at once
  n_samples <- ncol(edge_matrix)
  
  # Count non-NA values per row
  n_valid <- rowSums(!is.na(edge_matrix))
  
  # Calculate standard deviation per row (using only non-NA values)
  sd_edge <- apply(edge_matrix, 1, sd, na.rm = TRUE)
  
  # Calculate standard error
  se <- sd_edge / sqrt(n_valid)
  
  # Apply SAM-style variance moderation if requested
  if (moderateVariance) {
    s0 <- median(se, na.rm = TRUE)
    se <- se + s0
  }
  
  # Calculate t-statistic: t = (mean - mu) / se
  # For one-sample test against mu = 0
  test_stats <- meanEdge / se
  
  # Calculate p-values from t-distribution
  # df = n - 1
  df <- n_valid - 1
  
  pvalues <- switch(alternative,
    "two.sided" = 2 * pt(abs(test_stats), df = df, lower.tail = FALSE),
    "greater" = pt(test_stats, df = df, lower.tail = FALSE),
    "less" = pt(test_stats, df = df, lower.tail = TRUE)
  )
  
  # Set NA for edges with insufficient data
  # When moderateVariance is TRUE, se=0 is still valid because s0 is added
  insufficient_data <- n_valid < 2 | is.na(sd_edge) | (!moderateVariance & sd_edge == 0)
  test_stats[insufficient_data] <- NA
  pvalues[insufficient_data] <- NA
  
  # Compile results
  results <- data.frame(
    tf = tf_target$tf,
    target = tf_target$target,
    meanEdge = meanEdge,
    tStatistic = test_stats,
    pValue = pvalues,
    stringsAsFactors = FALSE
  )
  
  return(results)
}

#' @keywords internal
testEdgesTwoSample <- function(networksDF, group1, group2, alternative, minLog2FC,
                               moderateVariance = TRUE) {
  
  # Extract edge weights for both groups
  edge_data1 <- networksDF[, group1, drop = FALSE]
  edge_data2 <- networksDF[, group2, drop = FALSE]
  edge_matrix1 <- as.matrix(edge_data1)
  edge_matrix2 <- as.matrix(edge_data2)
  
  # Calculate mean edge weights
  meanEdge1 <- rowMeans(edge_matrix1, na.rm = TRUE)
  meanEdge2 <- rowMeans(edge_matrix2, na.rm = TRUE)
  meanEdge <- (meanEdge1 + meanEdge2) / 2
  diffMean <- meanEdge1 - meanEdge2
  
  # Calculate log2 fold change from quantiles (needed for filtering)
  # Convert z-scores to quantiles using pnorm with log.p=TRUE for precision
  # log2(q1/q2) = (log(q1) - log(q2)) / log(2)
  log2FC <- (pnorm(meanEdge1, log.p = TRUE) - pnorm(meanEdge2, log.p = TRUE)) / log(2)
  
  # Filter by minimum log2 fold change
  keep_idx <- abs(log2FC) >= minLog2FC
  edge_matrix1 <- edge_matrix1[keep_idx, , drop = FALSE]
  edge_matrix2 <- edge_matrix2[keep_idx, , drop = FALSE]
  meanEdge1 <- meanEdge1[keep_idx]
  meanEdge2 <- meanEdge2[keep_idx]
  meanEdge <- meanEdge[keep_idx]
  diffMean <- diffMean[keep_idx]
  log2FC <- log2FC[keep_idx]
  tf_target <- networksDF[keep_idx, c("tf", "target")]
  
  # Vectorized two-sample t-test computation (Welch's t-test)
  # Count non-NA values per row for each group
  n1 <- rowSums(!is.na(edge_matrix1))
  n2 <- rowSums(!is.na(edge_matrix2))
  
  # Calculate variance per row (using only non-NA values)
  var1 <- apply(edge_matrix1, 1, var, na.rm = TRUE)
  var2 <- apply(edge_matrix2, 1, var, na.rm = TRUE)
  
  # Calculate Welch's t-statistic: t = (mean1 - mean2) / sqrt(var1/n1 + var2/n2)
  se <- sqrt(var1/n1 + var2/n2)
  
  # Apply SAM-style variance moderation if requested
  if (moderateVariance) {
    s0 <- median(se, na.rm = TRUE)
    se <- se + s0
  }
  
  test_stats <- diffMean / se
  
  # Calculate degrees of freedom (Welch-Satterthwaite equation)
  df <- (var1/n1 + var2/n2)^2 / ((var1/n1)^2/(n1-1) + (var2/n2)^2/(n2-1))
  
  # Calculate p-values from t-distribution
  pvalues <- switch(alternative,
    "two.sided" = 2 * pt(abs(test_stats), df = df, lower.tail = FALSE),
    "greater" = pt(test_stats, df = df, lower.tail = FALSE),
    "less" = pt(test_stats, df = df, lower.tail = TRUE)
  )
  
  # Set NA for edges with insufficient data
  # When moderateVariance is TRUE, var=0 is still valid because s0 is added
  se_before_mod <- sqrt(var1/n1 + var2/n2)
  insufficient_data <- n1 < 2 | n2 < 2 | is.na(var1) | is.na(var2) | 
                       (!moderateVariance & (var1 == 0 | var2 == 0 | se_before_mod == 0)) | is.na(se)
  test_stats[insufficient_data] <- NA
  pvalues[insufficient_data] <- NA
  df[insufficient_data] <- NA
  
  # log2FC was already calculated before filtering
  
  # Calculate Cohen's d (standardized effect size)
  # d = (mean1 - mean2) / pooled_sd
  # pooled_sd = sqrt(((n1-1)*var1 + (n2-1)*var2) / (n1+n2-2))
  pooled_var <- ((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2)
  pooled_sd <- sqrt(pooled_var)
  cohensD <- diffMean / pooled_sd
  cohensD[pooled_sd == 0 | is.na(pooled_sd)] <- NA
  
  # Compile results
  results <- data.frame(
    tf = tf_target$tf,
    target = tf_target$target,
    meanGroup1 = meanEdge1,
    meanGroup2 = meanEdge2,
    diffMean = diffMean,
    cohensD = cohensD,
    log2FoldChange = log2FC,
    meanEdge = meanEdge,
    tStatistic = test_stats,
    pValue = pvalues,
    stringsAsFactors = FALSE
  )
  
  return(results)
}

#' @keywords internal
testEdgesPaired <- function(networksDF, group1, group2, alternative, minLog2FC,
                            moderateVariance = TRUE) {
  
  # Extract edge weights for both groups
  edge_data1 <- networksDF[, group1, drop = FALSE]
  edge_data2 <- networksDF[, group2, drop = FALSE]
  edge_matrix1 <- as.matrix(edge_data1)
  edge_matrix2 <- as.matrix(edge_data2)
  
  # Calculate mean edge weights
  meanEdge1 <- rowMeans(edge_matrix1, na.rm = TRUE)
  meanEdge2 <- rowMeans(edge_matrix2, na.rm = TRUE)
  meanEdge <- (meanEdge1 + meanEdge2) / 2
  
  # Calculate paired differences
  diff_matrix <- edge_matrix1 - edge_matrix2
  diffMean <- rowMeans(diff_matrix, na.rm = TRUE)
  
  # Calculate log2 fold change from quantiles (needed for filtering)
  # Convert z-scores to quantiles using pnorm with log.p=TRUE for precision
  # log2(q1/q2) = (log(q1) - log(q2)) / log(2)
  log2FC <- (pnorm(meanEdge1, log.p = TRUE) - pnorm(meanEdge2, log.p = TRUE)) / log(2)
  
  # Filter by minimum log2 fold change
  keep_idx <- abs(log2FC) >= minLog2FC
  diff_matrix <- diff_matrix[keep_idx, , drop = FALSE]
  edge_matrix1 <- edge_matrix1[keep_idx, , drop = FALSE]
  edge_matrix2 <- edge_matrix2[keep_idx, , drop = FALSE]
  meanEdge1 <- meanEdge1[keep_idx]
  meanEdge2 <- meanEdge2[keep_idx]
  meanEdge <- meanEdge[keep_idx]
  diffMean <- diffMean[keep_idx]
  log2FC <- log2FC[keep_idx]
  tf_target <- networksDF[keep_idx, c("tf", "target")]
  
  # Vectorized paired t-test computation
  # For paired t-test: t = mean(d) / (sd(d) / sqrt(n))
  # where d = differences between pairs
  
  # Count valid pairs per row (both values must be non-NA)
  valid_pairs <- !is.na(edge_matrix1) & !is.na(edge_matrix2)
  n_valid <- rowSums(valid_pairs)
  
  # Calculate standard deviation of differences
  sd_diff <- apply(diff_matrix, 1, sd, na.rm = TRUE)
  
  # Calculate standard error
  se <- sd_diff / sqrt(n_valid)
  
  # Apply SAM-style variance moderation if requested
  if (moderateVariance) {
    s0 <- median(se, na.rm = TRUE)
    se <- se + s0
  }
  
  # Calculate t-statistic
  test_stats <- diffMean / se
  
  # Calculate degrees of freedom: df = n - 1
  df <- n_valid - 1
  
  # Calculate p-values from t-distribution
  pvalues <- switch(alternative,
    "two.sided" = 2 * pt(abs(test_stats), df = df, lower.tail = FALSE),
    "greater" = pt(test_stats, df = df, lower.tail = FALSE),
    "less" = pt(test_stats, df = df, lower.tail = TRUE)
  )
  
  # Set NA for edges with insufficient data
  # When moderateVariance is TRUE, sd_diff=0 is still valid because s0 is added
  insufficient_data <- n_valid < 2 | is.na(sd_diff) | (!moderateVariance & sd_diff == 0)
  test_stats[insufficient_data] <- NA
  pvalues[insufficient_data] <- NA
  
  # log2FC was already calculated before filtering
  
  # Calculate Cohen's d for paired data
  # d = mean(diff) / sd(diff)
  cohensD <- diffMean / sd_diff
  cohensD[sd_diff == 0 | is.na(sd_diff)] <- NA
  
  # Compile results
  results <- data.frame(
    tf = tf_target$tf,
    target = tf_target$target,
    meanGroup1 = meanEdge1,
    meanGroup2 = meanEdge2,
    diffMean = diffMean,
    cohensD = cohensD,
    log2FoldChange = log2FC,
    meanEdge = meanEdge,
    tStatistic = test_stats,
    pValue = pvalues,
    stringsAsFactors = FALSE
  )
  
  return(results)
}

#' @title Regression analysis of edges across ordered conditions
#' @description Performs linear regression on network edges from runSCORPION output
#' to identify edges that show significant trends across ordered conditions (e.g.,
#' disease progression: Normal -> Border -> Tumor).
#' @param networksDF A data.frame output from \code{\link{runSCORPION}} containing
#'   TF-target pairs as rows and network identifiers as columns.
#' @param orderedGroups A named list where each element is a character vector of
#'   column names in \code{networksDF}. Names represent ordered conditions
#'   (e.g., list(Normal = c("P31--N", "P32--N"), Border = c("P31--B", "P32--B"),
#'   Tumor = c("P31--T", "P32--T"))). The order of list elements defines the
#'   progression (first to last).
#' @param padjustMethod Character specifying the p-value adjustment method for multiple
#'   testing correction. See \code{\link[stats]{p.adjust}} for options. Default "BH"
#'   (Benjamini-Hochberg FDR).
#' @param minMeanEdge Numeric threshold for minimum mean absolute edge weight to
#'   include in testing. Edges with mean absolute weight below this threshold are
#'   excluded. Default 0 (no filtering).
#' @return A data.frame containing:
#'   \itemize{
#'     \item{tf: Transcription factor}
#'     \item{target: Target gene}
#'     \item{slope: Regression slope (change in edge weight per condition step)}
#'     \item{intercept: Regression intercept}
#'     \item{rSquared: R-squared value (proportion of variance explained)}
#'     \item{fStatistic: F-statistic for the regression}
#'     \item{pValue: Raw p-value for the slope}
#'     \item{pAdj: Adjusted p-value}
#'     \item{meanEdge: Overall mean edge weight across all conditions}
#'     \item{One column per condition showing mean edge weight in that condition}
#'   }
#' @details
#' This function performs simple linear regression for each edge, modeling edge weight
#' as a function of an ordered categorical variable (coded as 0, 1, 2, ... for each
#' condition level).
#'
#' The slope coefficient indicates the average change in edge weight per step along
#' the ordered progression. Positive slopes indicate increasing edge weights,
#' negative slopes indicate decreasing edge weights.
#'
#' The function uses vectorized computations for efficiency with large datasets.
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
#' # Define ordered progression: Normal -> Border -> Tumor
#' normal_nets <- grep("--N$", colnames(nets), value = TRUE)
#' border_nets <- grep("--B$", colnames(nets), value = TRUE)
#' tumor_nets <- grep("--T$", colnames(nets), value = TRUE)
#'
#' ordered_conditions <- list(
#'   Normal = normal_nets,
#'   Border = border_nets,
#'   Tumor = tumor_nets
#' )
#'
#' # Perform regression analysis
#' results_regression <- regressEdges(
#'   networksDF = nets,
#'   orderedGroups = ordered_conditions
#' )
#'
#' # View top edges with strongest trends
#' head(results_regression[order(results_regression$pAdj), ])
#'
#' # Edges with positive slopes (increasing from N to T)
#' increasing <- results_regression[results_regression$pAdj < 0.05 & 
#'                                   results_regression$slope > 0, ]
#' print(paste("Edges increasing along N->B->T:", nrow(increasing)))
#'
#' # Edges with negative slopes (decreasing from N to T)
#' decreasing <- results_regression[results_regression$pAdj < 0.05 & 
#'                                   results_regression$slope < 0, ]
#' print(paste("Edges decreasing along N->B->T:", nrow(decreasing)))
#'
#' # Filter by minimum edge weight and R-squared
#' strong_trends <- results_regression[results_regression$pAdj < 0.05 & 
#'                                      results_regression$rSquared > 0.7 &
#'                                      abs(results_regression$meanEdge) > 0.1, ]
#' }
#' @export
#' @importFrom stats lm coef summary.lm p.adjust
regressEdges <- function(networksDF,
                         orderedGroups,
                         padjustMethod = "BH",
                         minMeanEdge = 0) {
  
  # Input validation
  if (missing(orderedGroups) || is.null(orderedGroups)) {
    stop("orderedGroups must be specified")
  }
  
  if (!is.list(orderedGroups) || is.null(names(orderedGroups))) {
    stop("orderedGroups must be a named list")
  }
  
  if (length(orderedGroups) < 2) {
    stop("orderedGroups must contain at least 2 conditions")
  }
  
  # Validate all columns exist
  all_cols <- unlist(orderedGroups, use.names = FALSE)
  if (!all(all_cols %in% colnames(networksDF))) {
    missing_cols <- setdiff(all_cols, colnames(networksDF))
    stop("Some columns not found in networksDF: ", paste(missing_cols, collapse = ", "))
  }
  
  # Extract TF and target columns
  tf_col <- which(colnames(networksDF) == "tf")
  target_col <- which(colnames(networksDF) == "target")
  
  if (length(tf_col) == 0 || length(target_col) == 0) {
    stop("networksDF must contain 'tf' and 'target' columns")
  }
  
  # Prepare data for regression
  n_conditions <- length(orderedGroups)
  condition_names <- names(orderedGroups)
  
  # Create predictor variable (0, 1, 2, ... for ordered conditions)
  x <- numeric()
  edge_data_list <- list()
  
  for (i in seq_along(orderedGroups)) {
    cols <- orderedGroups[[i]]
    edge_data_list[[i]] <- as.matrix(networksDF[, cols, drop = FALSE])
    x <- c(x, rep(i - 1, length(cols)))  # 0-indexed for first condition
  }
  
  # Combine all edge data
  edge_matrix <- do.call(cbind, edge_data_list)
  
  # Calculate mean edge weight across all conditions
  meanEdge <- rowMeans(edge_matrix, na.rm = TRUE)
  
  # Calculate mean for each condition
  condition_means <- matrix(NA, nrow = nrow(edge_matrix), ncol = n_conditions)
  for (i in seq_along(orderedGroups)) {
    condition_means[, i] <- rowMeans(edge_data_list[[i]], na.rm = TRUE)
  }
  colnames(condition_means) <- paste0("mean", condition_names)
  
  # Filter by minimum mean edge weight
  keep_idx <- abs(meanEdge) >= minMeanEdge
  edge_matrix <- edge_matrix[keep_idx, , drop = FALSE]
  meanEdge <- meanEdge[keep_idx]
  condition_means <- condition_means[keep_idx, , drop = FALSE]
  tf_target <- networksDF[keep_idx, c("tf", "target")]
  
  # Vectorized linear regression
  n_edges <- nrow(edge_matrix)
  n_samples <- ncol(edge_matrix)
  
  # Pre-compute regression components
  x_mean <- mean(x)
  x_centered <- x - x_mean
  sxx <- sum(x_centered^2)
  
  # Initialize result vectors
  slopes <- numeric(n_edges)
  intercepts <- numeric(n_edges)
  r_squared <- numeric(n_edges)
  f_stats <- numeric(n_edges)
  pvalues <- numeric(n_edges)
  
  # Compute regression for each edge
  for (i in 1:n_edges) {
    y <- edge_matrix[i, ]
    
    # Remove NA values
    valid_idx <- !is.na(y)
    y_valid <- y[valid_idx]
    x_valid <- x[valid_idx]
    n_valid <- length(y_valid)
    
    if (n_valid < 3) {
      slopes[i] <- NA
      intercepts[i] <- NA
      r_squared[i] <- NA
      f_stats[i] <- NA
      pvalues[i] <- NA
      next
    }
    
    # Compute regression coefficients
    x_valid_mean <- mean(x_valid)
    y_valid_mean <- mean(y_valid)
    x_valid_centered <- x_valid - x_valid_mean
    y_valid_centered <- y_valid - y_valid_mean
    
    sxy <- sum(x_valid_centered * y_valid_centered)
    sxx_valid <- sum(x_valid_centered^2)
    syy <- sum(y_valid_centered^2)
    
    # Slope and intercept
    slopes[i] <- sxy / sxx_valid
    intercepts[i] <- y_valid_mean - slopes[i] * x_valid_mean
    
    # R-squared
    ss_res <- sum((y_valid - (intercepts[i] + slopes[i] * x_valid))^2)
    ss_tot <- syy
    r_squared[i] <- 1 - (ss_res / ss_tot)
    
    # F-statistic and p-value
    df_reg <- 1
    df_res <- n_valid - 2
    ms_reg <- (ss_tot - ss_res) / df_reg
    ms_res <- ss_res / df_res
    
    if (ms_res > 0) {
      f_stats[i] <- ms_reg / ms_res
      pvalues[i] <- pf(f_stats[i], df1 = df_reg, df2 = df_res, lower.tail = FALSE)
    } else {
      f_stats[i] <- NA
      pvalues[i] <- NA
    }
  }
  
  # Adjust p-values
  pAdj <- p.adjust(pvalues, method = padjustMethod)
  
  # Compile results
  results <- data.frame(
    tf = tf_target$tf,
    target = tf_target$target,
    slope = slopes,
    intercept = intercepts,
    rSquared = r_squared,
    fStatistic = f_stats,
    pValue = pvalues,
    pAdj = pAdj,
    meanEdge = meanEdge,
    condition_means,
    stringsAsFactors = FALSE
  )
  
  return(results)
}
