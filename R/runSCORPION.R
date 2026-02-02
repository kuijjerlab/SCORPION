#' @title Run SCORPION across cell groups and return combined networks
#' @description Builds per-group regulatory networks by running \code{\link{scorpion}} on subsets of cells defined by \code{cellsMetadata} and combining the resulting networks into a wide-format data frame where each column corresponds to a network.
#' @param gexMatrix An expression dataset with genes in the rows and barcodes (cells) in the columns.
#' @param tfMotifs A motif dataset, a data.frame or a matrix containing 3 columns. Each row describes a motif associated with a transcription factor (column 1) a gene (column 2) and a score (column 3).
#' @param ppiNet A Protein-Protein-Interaction dataset, a data.frame or matrix containing 3 columns. Each row describes a protein-protein interaction between transcription factor 1 (column 1), transcription factor 2 (column 2) and a score (column 3).
#' @param cellsMetadata A data.frame with cell-level metadata; must contain columns specified in \code{groupBy}.
#' @param groupBy Character vector of one or more column names in \code{cellsMetadata} to use for grouping cells into networks.
#' @param normalizeData Boolean to indicate normalization of expression data. Default TRUE performs log normalization.
#' @param removeBatchEffect Boolean to indicate batch effect correction. Default FALSE.
#' @param batch Factor or vector giving batch assignment for each cell; required if \code{removeBatchEffect = TRUE}.
#' @param minCells Minimum number of cells per group required to build a network. Default is 30.
#' @param computingEngine Either 'cpu' or 'gpu'. Passed to \code{\link{scorpion}}.
#' @param nCores Number of processors to be used if BLAS or MPI is active.
#' @param gammaValue Graining level of data (proportion of number of single cells to super-cells). Default 10.
#' @param nPC Number of principal components to use for kNN network construction. Default 25.
#' @param assocMethod Association method. Must be one of 'pearson', 'spearman' or 'pcNet'. Default 'pearson'.
#' @param alphaValue Value to be used for update variable in PANDA. Default 0.1.
#' @param hammingValue Value at which to terminate the process based on Hamming distance. Default 0.001.
#' @param nIter Sets the maximum number of iterations PANDA can run before exiting. Default Inf.
#' @param outNet Character specifying which network to extract. Options include "regNet", "coregNet", "coopNet". Default "regNet".
#' @param zScaling Boolean to indicate use of Z-Scores in output. FALSE will use [0,1] scale. Default TRUE.
#' @param showProgress Boolean to indicate printing of output for algorithm progress. Default TRUE.
#' @param randomizationMethod Method by which to randomize gene expression matrix. Default "None". Must be one of "None", "within.gene", "by.gene".
#' @param scaleByPresent Boolean to indicate scaling of correlations by percentage of positive samples. Default FALSE.
#' @param filterExpr Boolean to indicate whether or not to remove genes with 0 expression across all cells. Default FALSE.
#' @return A data.frame in wide format where rows represent TF-target pairs (union across all networks) and columns represent network identifiers. Cell values are edge weights from the corresponding network.
#' @details
#' This function is a wrapper around \code{\link{scorpion}} that groups cells according to metadata columns, filters out groups with insufficient cells, runs network inference on each remaining group independently, and finally combines all resulting networks into a single wide-format data frame.
#' @examples
#' \dontrun{
#' # Load test data
#' data(scorpionTest)
#'
#' # Example 1: Group by single column (region)
#' nets_by_region <- runSCORPION(
#'   gexMatrix = scorpionTest$gex,
#'   tfMotifs = scorpionTest$tf,
#'   ppiNet = scorpionTest$ppi,
#'   cellsMetadata = scorpionTest$metadata,
#'   groupBy = "region"
#' )
#'
#' # -- SCORPION ----------------------------------------------------------------
#' # + Normalizing data (log scale)
#' # i 3 networks requested
#' # + 3 networks meet the minimum cell requirement (30)
#' # i Computing networks
#' # + Networks successfully constructed
#' # + Networks successfully combined
#'
#' # head(nets_by_region)
#' #                           tf target           T          B           N
#' # 1                       AATF  ACKR1 -0.31433856 -0.3569918 -0.33734920
#' # 2                       ABL1  ACKR1 -0.32915008 -0.3648895 -0.34437341
#' # 3                      ACSS2  ACKR1 -0.31418599 -0.3557854 -0.33663144
#' # 4                       ADNP  ACKR1  0.04105895  0.1109288  0.09910822
#' # 5                      AEBP2  ACKR1 -0.18964574 -0.2202269 -0.17558140
#' # 6 AEBP2_EED_EZH2_RBBP4_SUZ12  ACKR1 -0.31024700 -0.3508320 -0.33054519
#'
#' # Example 2: Group by single column (donor)
#' nets_by_donor <- runSCORPION(
#'   gexMatrix = scorpionTest$gex,
#'   tfMotifs = scorpionTest$tf,
#'   ppiNet = scorpionTest$ppi,
#'   cellsMetadata = scorpionTest$metadata,
#'   groupBy = "donor"
#' )
#'
#' # -- SCORPION ----------------------------------------------------------------
#' # + Normalizing data (log scale)
#' # i 3 networks requested
#' # + 3 networks meet the minimum cell requirement (30)
#' # i Computing networks
#' # + Networks successfully constructed
#' # + Networks successfully combined
#' # head(nets_by_donor)
#' #                           tf target         P31        P32         P33
#' # 1                       AATF  ACKR1 -0.34869366 -0.3557884 -0.35010835
#' # 2                       ABL1  ACKR1 -0.33724323 -0.3575331 -0.32875974
#' # 3                      ACSS2  ACKR1 -0.34569954 -0.3573108 -0.34980657
#' # 4                       ADNP  ACKR1  0.09933951  0.1045316  0.06046914
#' # 5                      AEBP2  ACKR1 -0.25111137 -0.2245655 -0.23157035
#' # 6 AEBP2_EED_EZH2_RBBP4_SUZ12  ACKR1 -0.34148264 -0.3518686 -0.34398594
#'
#' # Example 3: Group by two columns (donor and region)
#' nets_by_donor_region <- runSCORPION(
#'   gexMatrix = scorpionTest$gex,
#'   tfMotifs = scorpionTest$tf,
#'   ppiNet = scorpionTest$ppi,
#'   cellsMetadata = scorpionTest$metadata,
#'   groupBy = c("donor", "region")
#' )
#'
#' # -- SCORPION ----------------------------------------------------------------
#' # + Normalizing data (log scale)
#' # i 9 networks requested
#' # + 9 networks meet the minimum cell requirement (30)
#' # i Computing networks
#' # + Networks successfully constructed
#' # + Networks successfully combined
#' # head(nets_by_donor_region)
#' #                           tf target      P31--T      P31--B     P31--N
#' # 1                       AATF  ACKR1 -0.32634975 -0.33717677 -0.3442886
#' # 2                       ABL1  ACKR1 -0.34048759 -0.33890429 -0.3509986
#' # 3                      ACSS2  ACKR1 -0.32570697 -0.33600811 -0.3436603
#' # 4                       ADNP  ACKR1  0.07975735  0.05354279  0.1048301
#' # 5                      AEBP2  ACKR1 -0.21472437 -0.20545660 -0.1815737
#' # 6 AEBP2_EED_EZH2_RBBP4_SUZ12  ACKR1 -0.31861592 -0.32809314 -0.3375652
#'
#' # Example 4: Group by three columns (donor, region, and cell_type)
#' nets_by_donor_region_cell_type <- runSCORPION(
#'   gexMatrix = scorpionTest$gex,
#'   tfMotifs = scorpionTest$tf,
#'   ppiNet = scorpionTest$ppi,
#'   cellsMetadata = scorpionTest$metadata,
#'   groupBy = c("donor", "region", "cell_type")
#' )
#'
#' # -- SCORPION ----------------------------------------------------------------
#' # + Normalizing data (log scale)
#' # i 9 networks requested
#' # + 9 networks meet the minimum cell requirement (30)
#' # i Computing networks
#' # + Networks successfully constructed
#' # + Networks successfully combined
#' # head(nets_by_donor_region_cell_type)
#' #                           tf target P31--T--Epithelial P31--B--Epithelial
#' # 1                       AATF  ACKR1        -0.32634975        -0.33717677
#' # 2                       ABL1  ACKR1        -0.34048759        -0.33890429
#' # 3                      ACSS2  ACKR1        -0.32570697        -0.33600811
#' # 4                       ADNP  ACKR1         0.07975735         0.05354279
#' # 5                      AEBP2  ACKR1        -0.21472437        -0.20545660
#' # 6 AEBP2_EED_EZH2_RBBP4_SUZ12  ACKR1        -0.31861592        -0.32809314
#'
#' # Example 5: Using GPU computing engine (if available)
#' nets_gpu <- runSCORPION(
#'   gexMatrix = scorpionTest$gex,
#'   tfMotifs = scorpionTest$tf,
#'   ppiNet = scorpionTest$ppi,
#'   cellsMetadata = scorpionTest$metadata,
#'   groupBy = "region",
#'   computingEngine = "gpu"
#' )
#'
#' # -- SCORPION ----------------------------------------------------------------
#' # + Normalizing data (log scale)
#' # i 3 networks requested
#' # + 3 networks meet the minimum cell requirement (30)
#' # i Computing networks
#' # + Networks successfully constructed
#' # + Networks successfully combined
#' # head(nets_gpu)
#' #                           tf target           T          B           N
#' # 1                       AATF  ACKR1 -0.31433821 -0.3569913 -0.33734894
#' # 2                       ABL1  ACKR1 -0.32915005 -0.3648892 -0.34437302
#' # 3                      ACSS2  ACKR1 -0.31418574 -0.3557851 -0.33663106
#' # 4                       ADNP  ACKR1  0.04105883  0.1109285  0.09910798
#' # 5                      AEBP2  ACKR1 -0.18964562 -0.2202267 -0.17558131
#' # 6 AEBP2_EED_EZH2_RBBP4_SUZ12  ACKR1 -0.31024694 -0.3508317 -0.33054504
#'
#' # Example 6: Removing batch effect using donor as batch
#' nets_batch_corrected <- runSCORPION(
#'   gexMatrix = scorpionTest$gex,
#'   tfMotifs = scorpionTest$tf,
#'   ppiNet = scorpionTest$ppi,
#'   cellsMetadata = scorpionTest$metadata,
#'   groupBy = "region",
#'   removeBatchEffect = TRUE,
#'   batch = scorpionTest$metadata$donor
#' )
#'
#' # -- SCORPION ----------------------------------------------------------------
#' # + Normalizing data (log scale)
#' # + Correcting for batch effects
#' # i 3 networks requested
#' # + 3 networks meet the minimum cell requirement (30)
#' # i Computing networks
#' # + Networks successfully constructed
#' # + Networks successfully combined
#' # head(nets_batch_corrected)
#' #                           tf target          T           B           N
#' # 1                       AATF  ACKR1 -0.3337298 -0.34885471 -0.13011777
#' # 2                       ABL1  ACKR1 -0.3408020 -0.35409813 -0.17694266
#' # 3                      ACSS2  ACKR1 -0.3325270 -0.35115311 -0.12661518
#' # 4                       ADNP  ACKR1  0.1117504  0.08691481  0.01608898
#' # 5                      AEBP2  ACKR1 -0.2334648 -0.22113011  0.12519312
#' # 6 AEBP2_EED_EZH2_RBBP4_SUZ12  ACKR1 -0.3274770 -0.34475499 -0.12449908
#' }
#' @export
#' @importFrom dplyr %>% mutate group_by filter bind_rows .data
#' @importFrom stats reshape model.matrix
#' @importFrom cli cli_h1 cli_alert_success cli_alert_info cli_abort cli_progress_along
runSCORPION <- function(gexMatrix,
                        tfMotifs,
                        ppiNet,
                        cellsMetadata,
                        groupBy,
                        normalizeData = TRUE,
                        removeBatchEffect = FALSE,
                        batch = NULL,
                        minCells = 30,
                        computingEngine = "cpu",
                        nCores = 1,
                        gammaValue = 10,
                        nPC = 25,
                        assocMethod = "pearson",
                        alphaValue = 0.1,
                        hammingValue = 0.001,
                        nIter = Inf,
                        outNet = "regNet",
                        zScaling = TRUE,
                        showProgress = TRUE,
                        randomizationMethod = "None",
                        scaleByPresent = FALSE,
                        filterExpr = FALSE
) {
  # Input validation
  if (ncol(gexMatrix) != nrow(cellsMetadata)) {
    cli::cli_abort("gexMatrix must have the same number of columns as cellsMetadata has rows")
  }
  if (!all(groupBy %in% colnames(cellsMetadata))) {
    cli::cli_abort("groupBy columns not found in cellsMetadata: {paste(setdiff(groupBy, colnames(cellsMetadata)), collapse=', ')}")
  }
  if (!is.numeric(minCells) || minCells < 1) {
    cli::cli_abort("minCells must be a positive integer")
  }
  if (removeBatchEffect && is.null(batch)) {
    cli::cli_abort("batch must be provided when removeBatchEffect = TRUE")
  }
  if (!is.null(batch) && length(batch) != ncol(gexMatrix)) {
    cli::cli_abort("batch must have the same length as gexMatrix columns (cells)")
  }
  if (alphaValue < 0 || alphaValue > 1) {
    cli::cli_abort("alphaValue must be a numeric value between 0 and 1")
  }

  if (showProgress) {
    cli::cli_h1("SCORPION")
  }

  # Normalizing data
  if (normalizeData) {
    if (showProgress) {
      cli::cli_alert_success("Normalizing data (log scale)")
    }
    gexMatrix <- log_normalize_data(gexMatrix)
  }

  # Removing batch effect
  if (removeBatchEffect) {
    if (showProgress) {
      cli::cli_alert_success("Correcting for batch effects")
    }
    if (is.null(batch)) {
      cli::cli_abort('batch is needed for batch effect correction')
    }
    mean_expr <- apply(gexMatrix, 1, median)
    gexMatrix <- remove_batch(X = gexMatrix, batch = batch)
    gexMatrix <- gexMatrix + mean_expr
  }

  # Setting min number of cells to construct network
  min_cells <- max(minCells, 30)

  if (all(groupBy %in% colnames(cellsMetadata))) {
    collapsedGroup = apply(cellsMetadata[, groupBy, drop = FALSE], 1, function(X) { paste0(X, collapse = '--') })
    metadata <- data.frame(cell_id = colnames(gexMatrix), network_id = collapsedGroup)
    metadata <- metadata %>%
      group_by(.data$network_id) %>%
      mutate(n_cells = length(.data$cell_id))
  } else {
    cli::cli_abort('groupBy must match cellsMetadata column name')
  }

  total_net <- length(unique(metadata$network_id))
  metadata <- metadata %>% filter(.data$n_cells >= min_cells)
  filtered_net <- length(unique(metadata$network_id))

  if (showProgress) {
    cli::cli_alert_info(paste0(total_net, " networks requested"))
    cli::cli_alert_success(paste0(filtered_net, " networks meet the minimum cell requirement (", min_cells, ")"))
  }

  compute_network <- function(selected_network) {
    selected_network <- network_ids[selected_network]

    selected_cells <- metadata %>%
      filter(.data$network_id %in% selected_network)
    selected_cells <- gexMatrix[, selected_cells$cell_id]

    network <- scorpion(
      gexMatrix = selected_cells,
      tfMotifs = as.data.frame(tfMotifs),
      ppiNet = as.data.frame(ppiNet),
      computingEngine = computingEngine,
      nCores = nCores,
      gammaValue = gammaValue,
      nPC = nPC,
      assocMethod = assocMethod,
      alphaValue = alphaValue,
      hammingValue = hammingValue,
      nIter = nIter,
      outNet = outNet,
      zScaling = zScaling,
      showProgress = FALSE,
      randomizationMethod = randomizationMethod,
      scaleByPresent = scaleByPresent,
      filterExpr = filterExpr
    )[[outNet]]
    network <- as.data.frame(as.table(network))
    colnames(network) <- c("tf", "target", "weights")
    network$network_id <- selected_network
    return(network)
  }

  network_ids <- unique(metadata$network_id)

  if (showProgress) {
    cli::cli_alert_info("Computing networks")
    networks <- lapply(cli_progress_along(network_ids), compute_network)
    cli::cli_alert_success("Networks successfully constructed")
  } else {
    networks <- lapply(seq_along(network_ids), compute_network)
  }

  networks <- bind_rows(networks)
  networks <- reshape(
    networks,
    idvar = c("tf", "target"),
    timevar = "network_id",
    direction = "wide"
  )
  names(networks) <- sub("^weights\\.", "", names(networks))
  if (showProgress) {
    cli::cli_alert_success("Networks successfully combined")
  }
  return(networks)
}
