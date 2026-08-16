#' @title Regression analysis of edges across ordered conditions
#' @description Performs linear regression on network edges from runSCORPION output
#' to identify edges that show significant trends across ordered conditions (e.g.,
#' disease progression: Normal -> Border -> Tumor).
#' @author Daniel Osorio <daniecos@uio.no>
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
#' @seealso \code{\link{runSCORPION}}, \code{\link{testEdges}}
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
#' @importFrom stats pf p.adjust
regressEdges <- function(networksDF,
                         orderedGroups,
                         padjustMethod = "BH",
                         minMeanEdge = 0) {

  # Input validation
  if (missing(orderedGroups) || is.null(orderedGroups)) {
    cli::cli_abort("orderedGroups must be specified")
  }

  if (!is.list(orderedGroups) || is.null(names(orderedGroups))) {
    cli::cli_abort("orderedGroups must be a named list")
  }

  if (length(orderedGroups) < 2) {
    cli::cli_abort("orderedGroups must contain at least 2 conditions")
  }

  # Validate all columns exist
  all_cols <- unlist(orderedGroups, use.names = FALSE)
  if (!all(all_cols %in% colnames(networksDF))) {
    missing_cols <- setdiff(all_cols, colnames(networksDF))
    cli::cli_abort("Some columns not found in networksDF: {paste(missing_cols, collapse = ', ')}")
  }

  # Extract TF and target columns
  tf_col <- which(colnames(networksDF) == "tf")
  target_col <- which(colnames(networksDF) == "target")

  if (length(tf_col) == 0 || length(target_col) == 0) {
    cli::cli_abort("networksDF must contain 'tf' and 'target' columns")
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

  # Vectorized linear regression (no loop)
  n_edges <- nrow(edge_matrix)
  n_samples <- ncol(edge_matrix)

  # Build NA mask and zero-filled matrix for safe rowSums / matrix multiply
  mask <- !is.na(edge_matrix)
  edge_clean <- edge_matrix
  edge_clean[!mask] <- 0

  # Per-row valid counts
  n_valid <- rowSums(mask)

  # Per-row sums via matrix-vector products (mask converts NA positions to 0)
  sum_x  <- drop(mask %*% x)
  sum_x2 <- drop(mask %*% (x^2))
  sum_y  <- rowSums(edge_clean)
  sum_y2 <- rowSums(edge_clean^2)
  sum_xy <- drop(edge_clean %*% x)

  # Per-row means of x and y (over valid entries only)
  x_mean_r <- sum_x / n_valid
  y_mean_r <- sum_y / n_valid

  # Centered sums of squares and cross-products
  sxx <- sum_x2 - n_valid * x_mean_r^2
  sxy <- sum_xy - n_valid * x_mean_r * y_mean_r
  ss_tot <- sum_y2 - n_valid * y_mean_r^2  # = SYY

  # Regression coefficients
  slopes <- sxy / sxx
  intercepts <- y_mean_r - slopes * x_mean_r

  # SS_res via algebraic identity: SS_res = SS_tot - beta1^2 * Sxx
  ss_res <- ss_tot - slopes^2 * sxx

  # R-squared, F-statistic, p-value
  r_squared <- 1 - ss_res / ss_tot
  df_res <- n_valid - 2
  ms_res <- ss_res / df_res
  f_stats <- (ss_tot - ss_res) / ms_res
  pvalues <- pf(f_stats, df1 = 1, df2 = df_res, lower.tail = FALSE)

  # Mark edges with insufficient data or degenerate fits
  insufficient <- n_valid < 3 | sxx == 0 | ms_res <= 0
  slopes[insufficient] <- NA
  intercepts[insufficient] <- NA
  r_squared[insufficient] <- NA
  f_stats[insufficient] <- NA
  pvalues[insufficient] <- NA

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
