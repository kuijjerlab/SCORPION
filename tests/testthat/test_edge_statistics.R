# Tests for testEdges() and regressEdges() statistical accuracy
# Validates that vectorized implementations produce identical results to R's t.test and lm functions

# Create mock networksDF for testing
create_mock_network <- function(seed = 123) {
  set.seed(seed)
  n_edges <- 50
  n_samples_g1 <- 5
  n_samples_g2 <- 4
  
  # Create TF-target pairs
  tfs <- paste0("TF", 1:n_edges)
  targets <- paste0("Gene", 1:n_edges)
  
  # Create edge weights with some variation
  group1_cols <- paste0("G1_S", 1:n_samples_g1)
  group2_cols <- paste0("G2_S", 1:n_samples_g2)
  
  # Generate data with known differences
  g1_data <- matrix(rnorm(n_edges * n_samples_g1, mean = 0.5, sd = 0.3), 
                    nrow = n_edges, ncol = n_samples_g1)
  g2_data <- matrix(rnorm(n_edges * n_samples_g2, mean = 0.2, sd = 0.3), 
                    nrow = n_edges, ncol = n_samples_g2)
  
  colnames(g1_data) <- group1_cols
  colnames(g2_data) <- group2_cols
  
  df <- data.frame(
    tf = tfs,
    target = targets,
    g1_data,
    g2_data,
    stringsAsFactors = FALSE
  )
  
  list(
    df = df,
    group1 = group1_cols,
    group2 = group2_cols
  )
}

# Create mock data for regression testing (3 ordered conditions)
create_mock_network_regression <- function(seed = 123) {
  set.seed(seed)
  n_edges <- 30
  n_samples_per_condition <- 4
  
  tfs <- paste0("TF", 1:n_edges)
  targets <- paste0("Gene", 1:n_edges)
  
  # Create ordered conditions: Normal, Border, Tumor
  # Some edges increase, some decrease, some flat
  normal_cols <- paste0("N_S", 1:n_samples_per_condition)
  border_cols <- paste0("B_S", 1:n_samples_per_condition)
  tumor_cols <- paste0("T_S", 1:n_samples_per_condition)
  
  # Generate data with trends
  slopes <- seq(-0.5, 0.5, length.out = n_edges)  # Varying slopes
  
  normal_data <- matrix(NA, nrow = n_edges, ncol = n_samples_per_condition)
  border_data <- matrix(NA, nrow = n_edges, ncol = n_samples_per_condition)
  tumor_data <- matrix(NA, nrow = n_edges, ncol = n_samples_per_condition)
  
  for (i in 1:n_edges) {
    base <- rnorm(1, mean = 0, sd = 0.5)
    noise <- 0.1
    normal_data[i, ] <- base + 0 * slopes[i] + rnorm(n_samples_per_condition, 0, noise)
    border_data[i, ] <- base + 1 * slopes[i] + rnorm(n_samples_per_condition, 0, noise)
    tumor_data[i, ] <- base + 2 * slopes[i] + rnorm(n_samples_per_condition, 0, noise)
  }
  
  colnames(normal_data) <- normal_cols
  colnames(border_data) <- border_cols
  colnames(tumor_data) <- tumor_cols
  
  df <- data.frame(
    tf = tfs,
    target = targets,
    normal_data,
    border_data,
    tumor_data,
    stringsAsFactors = FALSE
  )
  
  list(
    df = df,
    orderedGroups = list(
      Normal = normal_cols,
      Border = border_cols,
      Tumor = tumor_cols
    )
  )
}

# =============================================================================
# Tests for single-sample t-test
# =============================================================================

test_that("testEdges single-sample t-statistic matches t.test()", {
  mock <- create_mock_network()
  
  # Run testEdges
  results <- testEdges(
    networksDF = mock$df,
    testType = "single",
    group1 = mock$group1
  )
  
  # Compare with t.test for first 10 edges
  for (i in 1:10) {
    edge_vals <- as.numeric(mock$df[i, mock$group1])
    
    # R's t.test
    ttest_result <- t.test(edge_vals, mu = 0, alternative = "two.sided")
    
    # Find corresponding result
    result_row <- results[results$tf == mock$df$tf[i] & 
                           results$target == mock$df$target[i], ]
    
    # Compare statistics (allow small numerical tolerance)
    expect_equal(result_row$tStatistic, as.numeric(ttest_result$statistic), 
                 tolerance = 1e-10,
                 label = paste("t-statistic for edge", i))
    expect_equal(result_row$pValue, ttest_result$p.value, 
                 tolerance = 1e-10,
                 label = paste("p-value for edge", i))
  }
})

test_that("testEdges single-sample one-sided tests match t.test()", {
  mock <- create_mock_network()
  
  # Test "greater" alternative
  results_greater <- testEdges(
    networksDF = mock$df,
    testType = "single",
    group1 = mock$group1,
    alternative = "greater"
  )
  
  # Test "less" alternative
  results_less <- testEdges(
    networksDF = mock$df,
    testType = "single",
    group1 = mock$group1,
    alternative = "less"
  )
  
  # Compare with t.test for first 5 edges
  for (i in 1:5) {
    edge_vals <- as.numeric(mock$df[i, mock$group1])
    
    ttest_greater <- t.test(edge_vals, mu = 0, alternative = "greater")
    ttest_less <- t.test(edge_vals, mu = 0, alternative = "less")
    
    result_greater <- results_greater[results_greater$tf == mock$df$tf[i] & 
                                       results_greater$target == mock$df$target[i], ]
    result_less <- results_less[results_less$tf == mock$df$tf[i] & 
                                 results_less$target == mock$df$target[i], ]
    
    expect_equal(result_greater$pValue, ttest_greater$p.value, 
                 tolerance = 1e-10,
                 label = paste("p-value (greater) for edge", i))
    expect_equal(result_less$pValue, ttest_less$p.value, 
                 tolerance = 1e-10,
                 label = paste("p-value (less) for edge", i))
  }
})

# =============================================================================
# Tests for two-sample t-test
# =============================================================================

test_that("testEdges two-sample t-statistic matches t.test() Welch", {
  mock <- create_mock_network()
  
  # Run testEdges
  results <- testEdges(
    networksDF = mock$df,
    testType = "two.sample",
    group1 = mock$group1,
    group2 = mock$group2
  )
  
  # Compare with t.test for first 10 edges
  for (i in 1:10) {
    edge_vals_g1 <- as.numeric(mock$df[i, mock$group1])
    edge_vals_g2 <- as.numeric(mock$df[i, mock$group2])
    
    # R's t.test (Welch's t-test, default)
    ttest_result <- t.test(edge_vals_g1, edge_vals_g2, alternative = "two.sided")
    
    # Find corresponding result
    result_row <- results[results$tf == mock$df$tf[i] & 
                           results$target == mock$df$target[i], ]
    
    # Compare statistics
    expect_equal(result_row$tStatistic, as.numeric(ttest_result$statistic), 
                 tolerance = 1e-10,
                 label = paste("t-statistic for edge", i))
    expect_equal(result_row$pValue, ttest_result$p.value, 
                 tolerance = 1e-10,
                 label = paste("p-value for edge", i))
  }
})

test_that("testEdges two-sample one-sided tests match t.test()", {
  mock <- create_mock_network()
  
  # Test "greater" alternative
  results_greater <- testEdges(
    networksDF = mock$df,
    testType = "two.sample",
    group1 = mock$group1,
    group2 = mock$group2,
    alternative = "greater"
  )
  
  # Test "less" alternative
  results_less <- testEdges(
    networksDF = mock$df,
    testType = "two.sample",
    group1 = mock$group1,
    group2 = mock$group2,
    alternative = "less"
  )
  
  # Compare with t.test for first 5 edges
  for (i in 1:5) {
    edge_vals_g1 <- as.numeric(mock$df[i, mock$group1])
    edge_vals_g2 <- as.numeric(mock$df[i, mock$group2])
    
    ttest_greater <- t.test(edge_vals_g1, edge_vals_g2, alternative = "greater")
    ttest_less <- t.test(edge_vals_g1, edge_vals_g2, alternative = "less")
    
    result_greater <- results_greater[results_greater$tf == mock$df$tf[i] & 
                                       results_greater$target == mock$df$target[i], ]
    result_less <- results_less[results_less$tf == mock$df$tf[i] & 
                                 results_less$target == mock$df$target[i], ]
    
    expect_equal(result_greater$pValue, ttest_greater$p.value, 
                 tolerance = 1e-10,
                 label = paste("p-value (greater) for edge", i))
    expect_equal(result_less$pValue, ttest_less$p.value, 
                 tolerance = 1e-10,
                 label = paste("p-value (less) for edge", i))
  }
})

test_that("testEdges two-sample mean difference is correct", {
  mock <- create_mock_network()
  
  results <- testEdges(
    networksDF = mock$df,
    testType = "two.sample",
    group1 = mock$group1,
    group2 = mock$group2
  )
  
  # Check mean differences for first 10 edges
  for (i in 1:10) {
    edge_vals_g1 <- as.numeric(mock$df[i, mock$group1])
    edge_vals_g2 <- as.numeric(mock$df[i, mock$group2])
    
    expected_diff <- mean(edge_vals_g1) - mean(edge_vals_g2)
    
    result_row <- results[results$tf == mock$df$tf[i] & 
                           results$target == mock$df$target[i], ]
    
    expect_equal(result_row$diffMean, expected_diff, 
                 tolerance = 1e-10,
                 label = paste("mean difference for edge", i))
    expect_equal(result_row$meanGroup1, mean(edge_vals_g1), 
                 tolerance = 1e-10,
                 label = paste("meanGroup1 for edge", i))
    expect_equal(result_row$meanGroup2, mean(edge_vals_g2), 
                 tolerance = 1e-10,
                 label = paste("meanGroup2 for edge", i))
  }
})

# =============================================================================
# Tests for paired t-test
# =============================================================================

# Create mock data for paired testing (equal sample sizes in matched order)
create_mock_network_paired <- function(seed = 456) {
  set.seed(seed)
  n_edges <- 50
  n_pairs <- 6  # Same number of samples in each group
  
  tfs <- paste0("TF", 1:n_edges)
  targets <- paste0("Gene", 1:n_edges)
  
  # Create paired samples (e.g., same patient, different conditions)
  # P1--T paired with P1--N, P2--T paired with P2--N, etc.
  group1_cols <- paste0("P", 1:n_pairs, "--T")
  group2_cols <- paste0("P", 1:n_pairs, "--N")
  
  # Generate correlated paired data
  g1_data <- matrix(NA, nrow = n_edges, ncol = n_pairs)
  g2_data <- matrix(NA, nrow = n_edges, ncol = n_pairs)
  
  for (i in 1:n_edges) {
    # Base value for each patient (creates correlation between pairs)
    patient_base <- rnorm(n_pairs, mean = 0, sd = 0.5)
    
    # Tumor adds a consistent effect (varying by edge)
    tumor_effect <- rnorm(1, mean = 0.3, sd = 0.2)
    
    g1_data[i, ] <- patient_base + tumor_effect + rnorm(n_pairs, 0, 0.1)
    g2_data[i, ] <- patient_base + rnorm(n_pairs, 0, 0.1)
  }
  
  colnames(g1_data) <- group1_cols
  colnames(g2_data) <- group2_cols
  
  df <- data.frame(
    tf = tfs,
    target = targets,
    g1_data,
    g2_data,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  list(
    df = df,
    group1 = group1_cols,
    group2 = group2_cols
  )
}

test_that("testEdges paired t-statistic matches t.test(paired = TRUE)", {
  mock <- create_mock_network_paired()
  
  # Run testEdges with paired = TRUE
  results <- testEdges(
    networksDF = mock$df,
    testType = "two.sample",
    group1 = mock$group1,
    group2 = mock$group2,
    paired = TRUE
  )
  
  # Compare with t.test(paired = TRUE) for first 15 edges
  for (i in 1:15) {
    edge_vals_g1 <- as.numeric(mock$df[i, mock$group1])
    edge_vals_g2 <- as.numeric(mock$df[i, mock$group2])
    
    # R's t.test with paired = TRUE
    ttest_result <- t.test(edge_vals_g1, edge_vals_g2, paired = TRUE, 
                           alternative = "two.sided")
    
    # Find corresponding result
    result_row <- results[results$tf == mock$df$tf[i] & 
                           results$target == mock$df$target[i], ]
    
    # Compare t-statistic
    expect_equal(result_row$tStatistic, as.numeric(ttest_result$statistic), 
                 tolerance = 1e-10,
                 label = paste("paired t-statistic for edge", i))
    
    # Compare p-value
    expect_equal(result_row$pValue, ttest_result$p.value, 
                 tolerance = 1e-10,
                 label = paste("paired p-value for edge", i))
  }
})

test_that("testEdges paired one-sided tests match t.test(paired = TRUE)", {
  mock <- create_mock_network_paired()
  
  # Test "greater" alternative
  results_greater <- testEdges(
    networksDF = mock$df,
    testType = "two.sample",
    group1 = mock$group1,
    group2 = mock$group2,
    paired = TRUE,
    alternative = "greater"
  )
  
  # Test "less" alternative
  results_less <- testEdges(
    networksDF = mock$df,
    testType = "two.sample",
    group1 = mock$group1,
    group2 = mock$group2,
    paired = TRUE,
    alternative = "less"
  )
  
  # Compare with t.test for first 10 edges
  for (i in 1:10) {
    edge_vals_g1 <- as.numeric(mock$df[i, mock$group1])
    edge_vals_g2 <- as.numeric(mock$df[i, mock$group2])
    
    ttest_greater <- t.test(edge_vals_g1, edge_vals_g2, paired = TRUE, 
                            alternative = "greater")
    ttest_less <- t.test(edge_vals_g1, edge_vals_g2, paired = TRUE, 
                         alternative = "less")
    
    result_greater <- results_greater[results_greater$tf == mock$df$tf[i] & 
                                       results_greater$target == mock$df$target[i], ]
    result_less <- results_less[results_less$tf == mock$df$tf[i] & 
                                 results_less$target == mock$df$target[i], ]
    
    expect_equal(result_greater$pValue, ttest_greater$p.value, 
                 tolerance = 1e-10,
                 label = paste("paired p-value (greater) for edge", i))
    expect_equal(result_less$pValue, ttest_less$p.value, 
                 tolerance = 1e-10,
                 label = paste("paired p-value (less) for edge", i))
  }
})

test_that("testEdges paired mean difference is correct", {
  mock <- create_mock_network_paired()
  
  results <- testEdges(
    networksDF = mock$df,
    testType = "two.sample",
    group1 = mock$group1,
    group2 = mock$group2,
    paired = TRUE
  )
  
  # Check that diffMean equals mean of differences (not difference of means)
  for (i in 1:10) {
    edge_vals_g1 <- as.numeric(mock$df[i, mock$group1])
    edge_vals_g2 <- as.numeric(mock$df[i, mock$group2])
    
    # For paired test, diffMean should be mean(g1 - g2)
    expected_diff <- mean(edge_vals_g1 - edge_vals_g2)
    
    result_row <- results[results$tf == mock$df$tf[i] & 
                           results$target == mock$df$target[i], ]
    
    expect_equal(result_row$diffMean, expected_diff, 
                 tolerance = 1e-10,
                 label = paste("paired mean difference for edge", i))
  }
})

test_that("testEdges paired vs unpaired gives different results", {
  mock <- create_mock_network_paired()
  
  results_paired <- testEdges(
    networksDF = mock$df,
    testType = "two.sample",
    group1 = mock$group1,
    group2 = mock$group2,
    paired = TRUE
  )
  
  results_unpaired <- testEdges(
    networksDF = mock$df,
    testType = "two.sample",
    group1 = mock$group1,
    group2 = mock$group2,
    paired = FALSE
  )
  
  # Results should be different (paired test accounts for within-subject correlation)
  expect_false(all(results_paired$tStatistic == results_unpaired$tStatistic, na.rm = TRUE))
  expect_false(all(results_paired$pValue == results_unpaired$pValue, na.rm = TRUE))
})

test_that("testEdges paired validation: requires equal length groups", {
  mock <- create_mock_network()
  
  # This should fail because group1 has 5 samples and group2 has 4 samples
  expect_error(
    testEdges(
      networksDF = mock$df,
      testType = "two.sample",
      group1 = mock$group1,
      group2 = mock$group2,
      paired = TRUE
    ),
    "group1 and group2 must have the same length"
  )
})

test_that("testEdges paired validation: paired requires two.sample test", {
  mock <- create_mock_network_paired()
  
  expect_error(
    testEdges(
      networksDF = mock$df,
      testType = "single",
      group1 = mock$group1,
      paired = TRUE
    ),
    "Paired tests require testType = 'two.sample'"
  )
})

test_that("testEdges paired handles NA values correctly", {
  mock <- create_mock_network_paired()
  
  # Introduce NA in both groups (same pair)
  mock$df[1, mock$group1[1]] <- NA
  mock$df[1, mock$group2[1]] <- NA
  
  results <- testEdges(
    networksDF = mock$df,
    testType = "two.sample",
    group1 = mock$group1,
    group2 = mock$group2,
    paired = TRUE
  )
  
  # Edge 1 should have n_pairs - 1 valid pairs
  edge_vals_g1 <- as.numeric(mock$df[1, mock$group1])
  edge_vals_g2 <- as.numeric(mock$df[1, mock$group2])
  
  # R's t.test handles NAs in paired test
  ttest_result <- t.test(edge_vals_g1, edge_vals_g2, paired = TRUE, 
                         alternative = "two.sided")
  
  result_row <- results[results$tf == mock$df$tf[1] & 
                         results$target == mock$df$target[1], ]
  
  expect_equal(result_row$tStatistic, as.numeric(ttest_result$statistic), 
               tolerance = 1e-10)
  expect_equal(result_row$pValue, ttest_result$p.value, 
               tolerance = 1e-10)
})

# =============================================================================
# Tests for regression (regressEdges)
# =============================================================================

test_that("regressEdges slope matches lm()", {
  mock <- create_mock_network_regression()
  
  # Run regressEdges
  results <- regressEdges(
    networksDF = mock$df,
    orderedGroups = mock$orderedGroups
  )
  
  # Compare with lm for first 10 edges
  for (i in 1:10) {
    # Get edge values for each condition
    normal_vals <- as.numeric(mock$df[i, mock$orderedGroups$Normal])
    border_vals <- as.numeric(mock$df[i, mock$orderedGroups$Border])
    tumor_vals <- as.numeric(mock$df[i, mock$orderedGroups$Tumor])
    
    # Create data for lm
    y <- c(normal_vals, border_vals, tumor_vals)
    x <- c(rep(0, length(normal_vals)), 
           rep(1, length(border_vals)), 
           rep(2, length(tumor_vals)))
    
    # R's lm
    lm_result <- lm(y ~ x)
    lm_summary <- summary(lm_result)
    
    # Find corresponding result
    result_row <- results[results$tf == mock$df$tf[i] & 
                           results$target == mock$df$target[i], ]
    
    # Compare slope
    expect_equal(result_row$slope, as.numeric(coef(lm_result)["x"]), 
                 tolerance = 1e-10,
                 label = paste("slope for edge", i))
    
    # Compare intercept
    expect_equal(result_row$intercept, as.numeric(coef(lm_result)["(Intercept)"]), 
                 tolerance = 1e-10,
                 label = paste("intercept for edge", i))
    
    # Compare R-squared
    expect_equal(result_row$rSquared, lm_summary$r.squared, 
                 tolerance = 1e-10,
                 label = paste("R-squared for edge", i))
    
    # Compare F-statistic
    f_stat <- lm_summary$fstatistic[1]
    expect_equal(result_row$fStatistic, as.numeric(f_stat), 
                 tolerance = 1e-10,
                 label = paste("F-statistic for edge", i))
    
    # Compare p-value (for the F-test)
    lm_pvalue <- pf(f_stat, lm_summary$fstatistic[2], lm_summary$fstatistic[3], 
                    lower.tail = FALSE)
    expect_equal(result_row$pValue, as.numeric(lm_pvalue), 
                 tolerance = 1e-10,
                 label = paste("p-value for edge", i))
  }
})

test_that("regressEdges condition means are correct", {
  mock <- create_mock_network_regression()
  
  results <- regressEdges(
    networksDF = mock$df,
    orderedGroups = mock$orderedGroups
  )
  
  # Check condition means for first 10 edges
  for (i in 1:10) {
    normal_vals <- as.numeric(mock$df[i, mock$orderedGroups$Normal])
    border_vals <- as.numeric(mock$df[i, mock$orderedGroups$Border])
    tumor_vals <- as.numeric(mock$df[i, mock$orderedGroups$Tumor])
    
    result_row <- results[results$tf == mock$df$tf[i] & 
                           results$target == mock$df$target[i], ]
    
    expect_equal(result_row$meanNormal, mean(normal_vals), 
                 tolerance = 1e-10,
                 label = paste("meanNormal for edge", i))
    expect_equal(result_row$meanBorder, mean(border_vals), 
                 tolerance = 1e-10,
                 label = paste("meanBorder for edge", i))
    expect_equal(result_row$meanTumor, mean(tumor_vals), 
                 tolerance = 1e-10,
                 label = paste("meanTumor for edge", i))
  }
})

test_that("regressEdges with 2 conditions matches lm()", {
  mock <- create_mock_network_regression()
  
  # Test with only 2 conditions
  ordered_2 <- list(
    Normal = mock$orderedGroups$Normal,
    Tumor = mock$orderedGroups$Tumor
  )
  
  results <- regressEdges(
    networksDF = mock$df,
    orderedGroups = ordered_2
  )
  
  # Compare with lm for first 5 edges
  for (i in 1:5) {
    normal_vals <- as.numeric(mock$df[i, ordered_2$Normal])
    tumor_vals <- as.numeric(mock$df[i, ordered_2$Tumor])
    
    y <- c(normal_vals, tumor_vals)
    x <- c(rep(0, length(normal_vals)), rep(1, length(tumor_vals)))
    
    lm_result <- lm(y ~ x)
    lm_summary <- summary(lm_result)
    
    result_row <- results[results$tf == mock$df$tf[i] & 
                           results$target == mock$df$target[i], ]
    
    expect_equal(result_row$slope, as.numeric(coef(lm_result)["x"]), 
                 tolerance = 1e-10,
                 label = paste("slope (2 conditions) for edge", i))
    expect_equal(result_row$rSquared, lm_summary$r.squared, 
                 tolerance = 1e-10,
                 label = paste("R-squared (2 conditions) for edge", i))
  }
})

# =============================================================================
# Tests for edge cases
# =============================================================================

test_that("testEdges handles NA values correctly", {
  mock <- create_mock_network()
  
  # Introduce some NAs
  mock$df[1, mock$group1[1]] <- NA
  mock$df[2, mock$group1[c(1, 2)]] <- NA
  
  results <- testEdges(
    networksDF = mock$df,
    testType = "single",
    group1 = mock$group1
  )
  
  # Edge 1: should have n-1 valid samples
  edge1_vals <- as.numeric(mock$df[1, mock$group1])
  edge1_vals <- edge1_vals[!is.na(edge1_vals)]
  ttest1 <- t.test(edge1_vals, mu = 0)
  
  result1 <- results[results$tf == mock$df$tf[1] & 
                      results$target == mock$df$target[1], ]
  
  expect_equal(result1$tStatistic, as.numeric(ttest1$statistic), tolerance = 1e-10)
  expect_equal(result1$pValue, ttest1$p.value, tolerance = 1e-10)
})

test_that("testEdges handles constant values (zero variance)", {
  mock <- create_mock_network()
  
  # Set one edge to constant value
  mock$df[1, mock$group1] <- 0.5
  
  results <- testEdges(
    networksDF = mock$df,
    testType = "single",
    group1 = mock$group1
  )
  
  result1 <- results[results$tf == mock$df$tf[1] & 
                      results$target == mock$df$target[1], ]
  
  # Should have NA for statistic and p-value

  expect_true(is.na(result1$tStatistic))
  expect_true(is.na(result1$pValue))
})

test_that("p-value adjustment is applied correctly", {
  mock <- create_mock_network()
  
  results <- testEdges(
    networksDF = mock$df,
    testType = "single",
    group1 = mock$group1,
    padjustMethod = "BH"
  )
  
  # Manually calculate adjusted p-values
  expected_padj <- p.adjust(results$pValue, method = "BH")
  
  expect_equal(results$pAdj, expected_padj, tolerance = 1e-10)
})

test_that("minMeanEdge filter works correctly", {
  mock <- create_mock_network()
  
  results_all <- testEdges(
    networksDF = mock$df,
    testType = "single",
    group1 = mock$group1,
    minMeanEdge = 0
  )
  
  results_filtered <- testEdges(
    networksDF = mock$df,
    testType = "single",
    group1 = mock$group1,
    minMeanEdge = 0.3
  )
  
  # Filtered results should have fewer rows
  expect_true(nrow(results_filtered) <= nrow(results_all))
  
  # All filtered edges should have |meanEdge| >= 0.3
  expect_true(all(abs(results_filtered$meanEdge) >= 0.3))
})
