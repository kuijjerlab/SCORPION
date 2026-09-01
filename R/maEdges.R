#' @title Meta-analysis of TF-target edges across studies
#' @description
#' Performs a meta-analysis of TF-target edges across multiple studies using
#' either a fixed-effect or DerSimonian-Laird random-effects model. Missing or
#' non-finite effect sizes and standard errors are excluded from the corresponding
#' study. A TF-target pair is only counted as contributing to a study when both
#' its effect size and SE are valid.
#'
#' @author Daniel Osorio <daniecos@uio.no>
#' @param edgesList A list of data.frames, one per study, typically produced by
#'   \code{\link{testEdges}}. Each data.frame must contain the columns
#'   \code{tf}, \code{target}, \code{log2FoldChange} and \code{SE}.
#' @param method Meta-analysis model. Either \code{"random"} (DerSimonian-Laird
#'   random-effects) or \code{"fixed"} (inverse-variance fixed-effect). Default
#'   is \code{"random"}.
#' @param minStudies Minimum number of studies with valid numeric information
#'   required for a TF-target pair to be included. Default 2.
#' @param padjustMethod Character specifying the p-value adjustment method for multiple
#'   testing correction. See \code{\link[stats]{p.adjust}} for options. Default "BH"
#'   (Benjamini-Hochberg FDR).
#' @param moderateVariance Logical indicating whether to apply SAM-style
#'   variance moderation to the meta-analysis SE. Default TRUE.
#' @param s0 Optional variance-moderation fudge factor. If NULL and
#'   \code{moderateVariance = TRUE}, the median of all valid meta-analysis SEs
#'   is used.
#' @return A data.frame containing:
#'   \itemize{
#'     \item{tf: Transcription factor}
#'     \item{target: Target gene}
#'     \item{k: Number of studies contributing to the meta-analysis}
#'     \item{log2FoldChange: Meta-analytic effect size}
#'     \item{SE: Meta-analysis standard error}
#'     \item{ciLow: Lower bound of the 95\% confidence interval}
#'     \item{ciHigh: Upper bound of the 95\% confidence interval}
#'     \item{zStatistic: Test statistic}
#'     \item{pValue: Raw p-value}
#'     \item{pAdj: Adjusted p-value}
#'     \item{Q: Cochran's Q heterogeneity statistic}
#'     \item{iSquared: I-squared heterogeneity (percentage)}
#'     \item{tauSquared: DerSimonian-Laird between-study variance}
#'   }
#' @seealso \code{\link{runSCORPION}}, \code{\link{testEdges}}, \code{\link{circosEdges}}
#' @export
#' @importFrom stats p.adjust qnorm pnorm median
#' @importFrom cli cli_alert_info cli_progress_bar cli_progress_update cli_progress_done
maEdges <- function(edgesList,
                    method = c("random", "fixed"),
                    minStudies = 2L,
                    padjustMethod = "BH",
                    moderateVariance = TRUE,
                    s0 = NULL) {
  
  method <- match.arg(method)
  n_studies <- length(edgesList)
  
  stopifnot(is.list(edgesList), n_studies >= 2L)
  
  # ============================================================
  # Build the UNION of TFs across all studies.
  # ============================================================
  
  tfs <- unique(unlist(lapply(edgesList, function(x) {
    unique(x$tf[!is.na(x$tf)])
  }), use.names = FALSE))
  
  n_tf <- length(tfs)
  
  cli_alert_info("TFs to process: {format(n_tf, big.mark = ',')}")
  
  # ============================================================
  # Output list
  # ============================================================
  
  result_list <- vector("list", n_tf)
  
  # ============================================================
  # Progress
  # ============================================================
  
  pb <- cli_progress_bar("TF-level meta-analysis",
                         total = n_tf,
                         clear = FALSE)
  
  # ============================================================
  # Process one TF at a time
  # ============================================================
  
  for (tf_i in seq_along(tfs)) {
    
    tf <- tfs[tf_i]
    
    # --------------------------------------------------------
    # Extract this TF from each study.
    #
    # SE is taken directly from testEdges().
    # --------------------------------------------------------
    
    study_data <- vector("list", n_studies)
    
    for (s in seq_len(n_studies)) {
      
      x <- edgesList[[s]]
      
      idx_tf <-
        !is.na(x$tf) &
        x$tf == tf
      
      if (!any(idx_tf)) {
        next
      }
      
      study_data[[s]] <- data.frame(
        target =
          x$target[idx_tf],
        effect =
          x$log2FoldChange[idx_tf],
        SE =
          x$SE[idx_tf],
        stringsAsFactors = FALSE
      )
    }
    
    # --------------------------------------------------------
    # Get target universe for this TF.
    # --------------------------------------------------------
    
    target_union <- unique(unlist(lapply(study_data, function(x) {
      
      if (is.null(x)) {
        NULL
      } else {
        x$target
      }
      
    }), use.names = FALSE))
    
    target_union <- target_union[!is.na(target_union)]
    
    if (!length(target_union)) {
      
      cli_progress_update(
        id = pb,
        status = sprintf("TF %d/%d", tf_i, n_tf)
      )
      
      next
    }
    
    n_targets <- length(target_union)
    
    # ========================================================
    # PASS 1 ACCUMULATORS
    # ========================================================
    
    k <- integer(n_targets)
    
    sum_w <- numeric(n_targets)
    sum_wy <- numeric(n_targets)
    sum_wy2 <- numeric(n_targets)
    sum_w2 <- numeric(n_targets)
    
    # ========================================================
    # PASS 1
    #
    # IMPORTANT:
    # Use the SE returned by testEdges() directly.
    # Do NOT reconstruct SE from p-values.
    # ========================================================
    
    for (s in seq_len(n_studies)) {
      
      x <- study_data[[s]]
      
      if (is.null(x)) {
        next
      }
      
      idx <- match(
        x$target,
        target_union
      )
      
      ok <- !is.na(idx)
      
      if (!any(ok)) {
        next
      }
      
      idx <- idx[ok]
      
      effect <- x$effect[ok]
      SE <- x$SE[ok]
      
      # ----------------------------------------------------
      # A study contributes only if BOTH effect and SE
      # are valid.
      # ----------------------------------------------------
      
      valid_input <-
        is.finite(effect) &
        is.finite(SE) &
        SE > 0
      
      if (!any(valid_input)) {
        next
      }
      
      idx <- idx[valid_input]
      effect <- effect[valid_input]
      SE <- SE[valid_input]
      
      # ----------------------------------------------------
      # Inverse-variance weight
      # ----------------------------------------------------
      
      w <- 1 / SE^2
      
      valid_w <-
        is.finite(w) &
        w > 0
      
      if (!any(valid_w)) {
        next
      }
      
      idx <- idx[valid_w]
      effect <- effect[valid_w]
      w <- w[valid_w]
      
      # ----------------------------------------------------
      # Accumulate
      # ----------------------------------------------------
      
      k[idx] <-
        k[idx] + 1L
      
      sum_w[idx] <-
        sum_w[idx] + w
      
      sum_wy[idx] <-
        sum_wy[idx] +
        w * effect
      
      sum_wy2[idx] <-
        sum_wy2[idx] +
        w * effect^2
      
      sum_w2[idx] <-
        sum_w2[idx] +
        w^2
    }
    
    # ========================================================
    # Keep targets with enough VALID studies
    # ========================================================
    
    valid_meta <-
      k >= minStudies &
      is.finite(sum_w) &
      sum_w > 0
    
    if (!any(valid_meta)) {
      
      cli_progress_update(
        id = pb,
        status = sprintf("TF %d/%d", tf_i, n_tf)
      )
      
      next
    }
    
    # --------------------------------------------------------
    # Restrict all vectors to meta-analyzed targets
    # --------------------------------------------------------
    
    target_meta <-
      target_union[valid_meta]
    
    k <-
      k[valid_meta]
    
    sum_w <-
      sum_w[valid_meta]
    
    sum_wy <-
      sum_wy[valid_meta]
    
    sum_wy2 <-
      sum_wy2[valid_meta]
    
    sum_w2 <-
      sum_w2[valid_meta]
    
    n_meta <-
      length(target_meta)
    
    # ========================================================
    # Fixed effect
    # ========================================================
    
    fixed_log2FC <-
      sum_wy / sum_w
    
    fixed_SE <-
      sqrt(1 / sum_w)
    
    # ========================================================
    # Cochran Q
    # ========================================================
    
    Q <-
      pmax(
        0,
        sum_wy2 -
          sum_wy^2 / sum_w
      )
    
    df <-
      k - 1L
    
    # ========================================================
    # DerSimonian-Laird tau²
    # ========================================================
    
    C <-
      sum_w -
      sum_w2 / sum_w
    
    tau2 <-
      numeric(n_meta)
    
    tau_ok <-
      k > 1L &
      is.finite(C) &
      C > 0
    
    tau2[tau_ok] <-
      pmax(
        0,
        (Q[tau_ok] - df[tau_ok]) /
          C[tau_ok]
      )
    
    # ========================================================
    # PASS 2: random effects
    # ========================================================
    
    if (method == "random") {
      
      sum_wr <-
        numeric(n_meta)
      
      sum_wry <-
        numeric(n_meta)
      
      for (s in seq_len(n_studies)) {
        
        x <- study_data[[s]]
        
        if (is.null(x)) {
          next
        }
        
        idx <- match(
          x$target,
          target_meta
        )
        
        ok <- !is.na(idx)
        
        if (!any(ok)) {
          next
        }
        
        idx <- idx[ok]
        
        effect <- x$effect[ok]
        SE <- x$SE[ok]
        
        # ----------------------------------------------------
        # Same validity criteria as PASS 1.
        # ----------------------------------------------------
        
        valid_input <-
          is.finite(effect) &
          is.finite(SE) &
          SE > 0
        
        if (!any(valid_input)) {
          next
        }
        
        idx <- idx[valid_input]
        effect <- effect[valid_input]
        SE <- SE[valid_input]
        
        # ----------------------------------------------------
        # Random-effects weight
        #
        # Raw study-level SE is used here.
        # ----------------------------------------------------
        
        wr <-
          1 /
          (SE^2 +
             tau2[idx])
        
        valid_w <-
          is.finite(wr) &
          wr > 0
        
        if (!any(valid_w)) {
          next
        }
        
        idx <- idx[valid_w]
        effect <- effect[valid_w]
        wr <- wr[valid_w]
        
        sum_wr[idx] <-
          sum_wr[idx] + wr
        
        sum_wry[idx] <-
          sum_wry[idx] +
          wr * effect
      }
      
      # ------------------------------------------------------
      # Random-effects estimate
      # ------------------------------------------------------
      
      meta_log2FC <-
        numeric(n_meta)
      
      meta_SE <-
        numeric(n_meta)
      
      random_ok <-
        is.finite(sum_wr) &
        sum_wr > 0
      
      meta_log2FC[random_ok] <-
        sum_wry[random_ok] /
        sum_wr[random_ok]
      
      meta_SE[random_ok] <-
        sqrt(
          1 /
            sum_wr[random_ok]
        )
      
    } else {
      
      meta_log2FC <-
        fixed_log2FC
      
      meta_SE <-
        fixed_SE
    }
    
    # ========================================================
    # Meta-analysis statistics
    #
    # NOTE:
    # The actual variance moderation is intentionally done
    # AFTER all TFs have been processed so that s0 is GLOBAL.
    #
    # Here we retain the raw meta-analysis SE.
    # ========================================================
    
    valid_effect <-
      is.finite(meta_log2FC) &
      is.finite(meta_SE) &
      meta_SE > 0
    
    if (!any(valid_effect)) {
      
      cli_progress_update(
        id = pb,
        status = sprintf("TF %d/%d", tf_i, n_tf)
      )
      
      next
    }
    
    # --------------------------------------------------------
    # Raw z-statistic
    #
    # This is retained temporarily. It will be replaced below
    # by the variance-moderated statistic after global s0 is
    # calculated.
    # --------------------------------------------------------
    
    meta_z <-
      numeric(n_meta)
    
    meta_z[valid_effect] <-
      meta_log2FC[valid_effect] /
      meta_SE[valid_effect]
    
    pValue <-
      numeric(n_meta)
    
    pValue[valid_effect] <-
      2 *
      pnorm(
        -abs(meta_z[valid_effect])
      )
    
    # --------------------------------------------------------
    # Confidence intervals
    # --------------------------------------------------------
    
    CI_low <-
      numeric(n_meta)
    
    CI_high <-
      numeric(n_meta)
    
    CI_low[valid_effect] <-
      meta_log2FC[valid_effect] -
      qnorm(.975) *
      meta_SE[valid_effect]
    
    CI_high[valid_effect] <-
      meta_log2FC[valid_effect] +
      qnorm(.975) *
      meta_SE[valid_effect]
    
    # --------------------------------------------------------
    # I²
    # --------------------------------------------------------
    
    I2 <-
      numeric(n_meta)
    
    i2_ok <-
      is.finite(Q) &
      Q > 0 &
      df > 0
    
    I2[i2_ok] <-
      pmax(
        0,
        (Q[i2_ok] - df[i2_ok]) /
          Q[i2_ok]
      ) * 100
    
    # ========================================================
    # Output
    # ========================================================
    
    out <- data.frame(
      tf =
        rep(
          tf,
          sum(valid_effect)
        ),
      
      target =
        target_meta[valid_effect],
      
      k =
        k[valid_effect],
      
      log2FoldChange =
        meta_log2FC[valid_effect],
      
      # Raw meta-analysis SE
      SE =
        meta_SE[valid_effect],
      
      ciLow =
        CI_low[valid_effect],
      
      ciHigh =
        CI_high[valid_effect],
      
      zStatistic =
        meta_z[valid_effect],
      
      pValue =
        pValue[valid_effect],
      
      Q =
        Q[valid_effect],
      
      iSquared =
        I2[valid_effect],
      
      tauSquared =
        tau2[valid_effect],
      
      stringsAsFactors = FALSE
    )
    
    # --------------------------------------------------------
    # Defensive final filter
    # --------------------------------------------------------
    
    out <-
      out[
        is.finite(out$log2FoldChange) &
          is.finite(out$SE) &
          is.finite(out$pValue),
        ,
        drop = FALSE
      ]
    
    if (nrow(out) > 0L) {
      result_list[[tf_i]] <- out
    }
    
    # --------------------------------------------------------
    # Cleanup
    # --------------------------------------------------------
    
    rm(
      study_data,
      target_union,
      target_meta,
      k,
      sum_w,
      sum_wy,
      sum_wy2,
      sum_w2,
      fixed_log2FC,
      fixed_SE,
      Q,
      df,
      C,
      tau2,
      meta_log2FC,
      meta_SE,
      meta_z,
      pValue,
      CI_low,
      CI_high,
      I2,
      out
    )
    
    gc(FALSE)
    
    # ========================================================
    # EXISTING PROGRESS BAR -- UNCHANGED
    # ========================================================
    
    cli_progress_update(
      id = pb,
      status = sprintf(
        "TF %d/%d",
        tf_i,
        n_tf
      )
    )
  }
  
  cli_progress_done(id = pb)
  
  # ============================================================
  # Combine
  # ============================================================
  
  cli_alert_info("Combining results")
  
  valid_results <-
    result_list[
      vapply(
        result_list,
        function(x) {
          !is.null(x) &&
            nrow(x) > 0L
        },
        logical(1)
      )
    ]
  
  if (!length(valid_results)) {
    
    return(
      data.frame(
        tf = character(),
        target = character(),
        k = integer(),
        log2FoldChange = numeric(),
        SE = numeric(),
        ciLow = numeric(),
        ciHigh = numeric(),
        zStatistic = numeric(),
        pValue = numeric(),
        pAdj = numeric(),
        Q = numeric(),
        iSquared = numeric(),
        tauSquared = numeric(),
        stringsAsFactors = FALSE
      )
    )
  }
  
  result <-
    do.call(
      rbind,
      valid_results
    )
  
  rownames(result) <- NULL
  
  # ============================================================
  # GLOBAL SAM-STYLE VARIANCE MODERATION
  #
  # This follows testEdges():
  #
  #     s0 <- median(SE)
  #     SE <- SE + s0
  #
  # Crucially, s0 is calculated ONCE over all meta-analysis
  # edges, rather than separately for each TF.
  # ============================================================
  
  valid_se <-
    is.finite(result$SE) &
    result$SE > 0
  
  if (moderateVariance && any(valid_se)) {
    
    if (is.null(s0)) {
      
      s0 <-
        median(
          result$SE[valid_se],
          na.rm = TRUE
        )
    }
    
    moderatedSE <-
      result$SE +
      s0
    
  } else {
    
    if (is.null(s0)) {
      s0 <- 0
    }
    
    moderatedSE <-
      result$SE
  }
  
  # ============================================================
  # Recalculate meta-analysis statistics using moderated SE
  # ============================================================
  
  valid_inference <-
    is.finite(result$log2FoldChange) &
    is.finite(moderatedSE) &
    moderatedSE > 0
  
  result$zStatistic <-
    NA_real_
  
  result$zStatistic[valid_inference] <-
    result$log2FoldChange[valid_inference] /
    moderatedSE[valid_inference]
  
  # ============================================================
  # P-values
  # ============================================================
  
  result$pValue <-
    NA_real_
  
  result$pValue[valid_inference] <-
    2 *
    pnorm(
      -abs(
        result$zStatistic[valid_inference]
      )
    )
  
  # ============================================================
  # Confidence intervals
  # ============================================================
  
  result$ciLow <-
    NA_real_
  
  result$ciHigh <-
    NA_real_
  
  result$ciLow[valid_inference] <-
    result$log2FoldChange[valid_inference] -
    qnorm(.975) *
    moderatedSE[valid_inference]
  
  result$ciHigh[valid_inference] <-
    result$log2FoldChange[valid_inference] +
    qnorm(.975) *
    moderatedSE[valid_inference]
  
  # ============================================================
  # BH correction
  # ============================================================
  
  result$pAdj <-
    p.adjust(
      result$pValue,
      method = padjustMethod
    )
  
  # ============================================================
  # Place pAdj next to pValue
  # ============================================================
  
  result <-
    result[
      ,
      c(
        "tf",
        "target",
        "k",
        "log2FoldChange",
        "SE",
        "ciLow",
        "ciHigh",
        "zStatistic",
        "pValue",
        "pAdj",
        "Q",
        "iSquared",
        "tauSquared"
      ),
      drop = FALSE
    ]
  
  # ============================================================
  # Sort safely
  # ============================================================
  
  result <-
    result[
      order(
        result$tf,
        result$target,
        na.last = TRUE
      ),
      ,
      drop = FALSE
    ]
  
  rownames(result) <- NULL
  
  gc()
  
  return(result)
}
