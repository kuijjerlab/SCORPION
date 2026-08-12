#' @title Gene set enrichment analysis of TF-target edges
#' @description
#' Performs gene set enrichment analysis separately for each transcription
#' factor (TF) using the edge-level values supplied in \code{numericValue}.
#' Enrichment is performed with the multilevel implementation of \pkg{fgsea}.
#' Calculations for individual TFs are performed in parallel.
#' @author Daniel Osorio <daniecos@uio.no>
#' @param edgesDF A data.frame of TF-target edges, typically produced by
#'   \code{\link{testEdges}}. Must contain a \code{tf} column, a \code{target}
#'   column, and the numeric column named by \code{numericValue}.
#' @param geneSets A named list of gene sets. The names of the list elements are
#'   used as gene set identifiers.
#' @param numericValue Character string naming the column in \code{edgesDF} used
#'   as the ranking statistic for enrichment analysis.
#' @param nCores Integer specifying the number of parallel workers to use.
#'   Default 3.
#' @param seed Integer specifying the random seed used by the parallel
#'   enrichment calculations. Default 1.
#' @return A data.frame of enrichment results with one row per TF-gene set pair:
#'   \itemize{
#'     \item{tf: Transcription factor}
#'     \item{geneSet: Gene set identifier}
#'     \item{pValue: Raw enrichment p-value}
#'     \item{pAdj: Benjamini-Hochberg adjusted p-value}
#'     \item{log2Err: Expected log2 error of the p-value estimate}
#'     \item{ES: Enrichment score}
#'     \item{NES: Normalized enrichment score}
#'     \item{geneSetSize: Number of genes from the set found among the targets}
#'   }
#' @details
#' For each TF, the values in \code{numericValue} are used as ranked statistics
#' for its target genes. Edges with missing targets or missing or non-finite
#' values in the selected numeric column are excluded before enrichment
#' analysis. If a target occurs more than once for a TF, only the observation
#' with the largest absolute value of the selected ranking statistic is
#' retained.
#'
#' Gene set enrichment is performed using \code{fgsea::fgseaMultilevel}; the
#' \pkg{fgsea} package (Bioconductor) is required. The \code{leadingEdge} column
#' returned by \code{fgseaMultilevel} is not included in the output. P-values
#' are adjusted across all TF-gene set enrichment tests using the
#' Benjamini-Hochberg procedure.
#' @seealso \code{\link{testEdges}}, \code{\link{maEdges}}
#' @examples
#' \dontrun{
#' data(scorpionTest)
#' nets <- runSCORPION(
#'   gexMatrix = scorpionTest$gex,
#'   tfMotifs = scorpionTest$tf,
#'   ppiNet = scorpionTest$ppi,
#'   cellsMetadata = scorpionTest$metadata,
#'   groupBy = c("donor", "region")
#' )
#' res <- testEdges(
#'   networksDF = nets,
#'   testType = "two.sample",
#'   group1 = grep("--T$", colnames(nets), value = TRUE),
#'   group2 = grep("--N$", colnames(nets), value = TRUE)
#' )
#'
#' geneSets <- list(SetA = c("ACKR1", "ACTA2"), SetB = c("ACTG2", "ADAMDEC1"))
#' enr <- enrichEdges(
#'   edgesDF = res,
#'   geneSets = geneSets,
#'   numericValue = "log2FoldChange"
#' )
#' }
#' @export
#' @importFrom stats p.adjust
enrichEdges <- function(edgesDF,
                        geneSets,
                        numericValue,
                        nCores = 3,
                        seed = 1) {
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop(
      "Package 'fgsea' (Bioconductor) is required for enrichEdges(). ",
      "Install it with BiocManager::install('fgsea').",
      call. = FALSE
    )
  }
  # ============================================================
  # Check input
  # ============================================================

  stopifnot(
    is.data.frame(edgesDF),
    is.list(geneSets),
    length(geneSets) > 0L,
    is.character(numericValue),
    length(numericValue) == 1L,
    is.numeric(nCores),
    length(nCores) == 1L,
    is.finite(nCores),
    nCores >= 1,
    nCores == as.integer(nCores),
    is.numeric(seed),
    length(seed) == 1L,
    is.finite(seed),
    seed == as.integer(seed)
  )

  required_columns <- c("tf", "target", numericValue)

  missing_columns <- setdiff(required_columns, colnames(edgesDF))

  if (length(missing_columns)) {
    stop(
      "Missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.numeric(edgesDF[[numericValue]])) {
    stop("`numericValue` must specify a numeric column.", call. = FALSE)
  }

  if (is.null(names(geneSets)) ||
      any(!nzchar(names(geneSets)))) {
    stop("`geneSets` must be a named list.", call. = FALSE)
  }

  # ============================================================
  # TFs
  # ============================================================

  tf_list <- unique(edgesDF$tf[!is.na(edgesDF$tf)])

  # ============================================================
  # Enrichment for one TF
  # ============================================================

  enrich_tf_edges <- function(selected_tf) {
    tf_edges <- edgesDF[
      !is.na(edgesDF$tf) & edgesDF$tf == selected_tf,
      ,
      drop = FALSE
    ]

    # ----------------------------------------------------------
    # Keep only edges with valid ranking values and targets.
    # ----------------------------------------------------------

    valid <- (
      !is.na(tf_edges$target) &
        is.finite(tf_edges[[numericValue]])
    )

    tf_edges <- tf_edges[valid, , drop = FALSE]

    if (!nrow(tf_edges)) {
      return(NULL)
    }

    # ----------------------------------------------------------
    # Construct ranked statistics.
    #
    # If a target occurs more than once for a TF, retain the
    # observation with the largest absolute statistic.
    # ----------------------------------------------------------

    tf_edges <- tf_edges[
      order(
        abs(tf_edges[[numericValue]]),
        decreasing = TRUE
      ),
      ,
      drop = FALSE
    ]

    tf_edges <- tf_edges[
      !duplicated(tf_edges$target),
      ,
      drop = FALSE
    ]

    weights <- tf_edges[[numericValue]]

    names(weights) <- tf_edges$target

    # ----------------------------------------------------------
    # fgsea
    # ----------------------------------------------------------

    enrichment <- suppressWarnings(
      fgsea::fgseaMultilevel(
        pathways = geneSets,
        stats = weights
      )
    )

    if (!nrow(enrichment)) {
      return(NULL)
    }

    enrichment <- as.data.frame(enrichment)

    # ----------------------------------------------------------
    # Keep the relevant output columns.
    # fgsea returns padj, which is replaced below by the global
    # adjustment across all TF-gene set tests.
    # ----------------------------------------------------------

    enrichment <- enrichment[
      setdiff(colnames(enrichment), "leadingEdge"),
      drop = FALSE
    ]

    enrichment <- enrichment[
      ,
      c(
        "pathway",
        "pval",
        "log2err",
        "ES",
        "NES",
        "size"
      ),
      drop = FALSE
    ]

    colnames(enrichment) <- c(
      "geneSet",
      "pValue",
      "log2Err",
      "ES",
      "NES",
      "geneSetSize"
    )

    enrichment$tf <- selected_tf

    enrichment <- enrichment[
      ,
      c(
        "tf",
        "geneSet",
        "pValue",
        "log2Err",
        "ES",
        "NES",
        "geneSetSize"
      ),
      drop = FALSE
    ]

    enrichment
  }

  # ============================================================
  # Parallel enrichment
  # ============================================================

  old_plan <- future::plan()

  on.exit(
    future::plan(old_plan),
    add = TRUE
  )

  future::plan(
    future::multisession,
    workers = nCores
  )

  enrichment <- furrr::future_map_dfr(
    tf_list,
    enrich_tf_edges,
    .progress = TRUE,
    .options = furrr::furrr_options(seed = seed)
  )

  # ============================================================
  # Handle case where no TF has valid information
  # ============================================================

  if (!nrow(enrichment)) {
    return(
      data.frame(
        tf = character(),
        geneSet = character(),
        pValue = numeric(),
        pAdj = numeric(),
        log2Err = numeric(),
        ES = numeric(),
        NES = numeric(),
        geneSetSize = integer(),
        stringsAsFactors = FALSE
      )
    )
  }

  # ============================================================
  # BH correction across all TF-gene set tests
  # ============================================================

  enrichment$pAdj <- p.adjust(
    enrichment$pValue,
    method = "BH"
  )

  # ============================================================
  # Final column order
  # ============================================================

  enrichment <- enrichment[
    ,
    c(
      "tf",
      "geneSet",
      "pValue",
      "pAdj",
      "log2Err",
      "ES",
      "NES",
      "geneSetSize"
    ),
    drop = FALSE
  ]

  rownames(enrichment) <- NULL

  enrichment
}
