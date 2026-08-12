#' @title Circos plot of differential network edges
#' @description Draws a circular (Circos) plot of transcription factor to target
#'   links from a \code{\link{testEdges}} two-sample result. Genes are placed on
#'   their genomic coordinates, links are coloured continuously by significance,
#'   flagged as known or novel against an optional a priori network, and genes
#'   belonging to supplied gene sets can be labelled around the circle.
#' @author Daniel Osorio <daniecos@uio.no>
#' @param edgesDF A data.frame produced by \code{\link{testEdges}} (two-sample or
#'   paired). Must contain the columns \code{tf}, \code{target},
#'   \code{log2FoldChange} and \code{pAdj} (\code{pValue} is used as a fallback
#'   when \code{pAdj} is absent).
#' @param species Ensembl dataset name passed to \pkg{biomaRt} when
#'   \code{geneCoords} is \code{NULL}, e.g. \code{"hsapiens_gene_ensembl"} for
#'   human or \code{"mmusculus_gene_ensembl"} for mouse. See
#'   \code{biomaRt::listDatasets()} for the full multi-species list.
#' @param geneCoords Optional data.frame supplying gene coordinates from any
#'   source (overrides the \pkg{biomaRt} download). Must have columns
#'   \code{gene}, \code{chr}, \code{start} and \code{end}.
#' @param priorNet Optional a priori TF-target network whose first two columns
#'   are the TF and target. Links present here are labelled \code{"known"}, all
#'   others \code{"novel"}. Accepts a data.frame or a matrix.
#' @param geneSets Optional gene-set annotation used to label genes around the
#'   circle: either a path to a GMT file or a named list of character vectors.
#' @param pAdjThreshold Numeric significance cutoff applied to \code{pAdj}
#'   (or \code{pValue} when \code{pAdj} is missing). Default \code{0.05}.
#' @param log2FCThreshold Numeric minimum absolute \code{log2FoldChange} required
#'   to draw a link. Default \code{0}.
#' @param maxEdges Integer cap on the number of links drawn; when exceeded, the
#'   most significant edges are kept. Default \code{500}.
#' @param colorBy Name of the \code{edgesDF} column mapped to the link colour
#'   ramp. Default \code{"log2FoldChange"}, giving a continuous diverging colour
#'   scale centred at zero. When the default is used but \code{log2FoldChange} is
#'   absent (e.g. single-sample \code{\link{testEdges}} output), it falls back to
#'   \code{meanEdge}.
#' @param linkColors Length-3 vector of colours for the low, mid and high ends of
#'   \code{colorBy}. Default blue-white-red; a diverging ramp is used when
#'   \code{colorBy} has negative values, otherwise a sequential low-to-high ramp.
#' @param lwdRange Length-2 numeric giving the minimum and maximum link line
#'   width; each link's thickness is scaled linearly within this range by its
#'   \code{-log10} adjusted p-value. Default \code{c(0.5, 4)}.
#' @param nmaxTF,nmaxTarget Integers giving how many TFs and targets to label,
#'   selected by the largest absolute out-degree and in-degree respectively.
#'   Use \code{NULL} or \code{Inf} to label all. Defaults \code{20}.
#' @param knownColor,novelColor Border colours distinguishing known from novel
#'   links. Defaults grey and orange.
#' @param geneSetColors Optional named vector mapping gene-set names to colours.
#'   When \code{NULL}, colours are generated automatically.
#' @param chromosomes Optional character vector restricting and ordering the
#'   chromosomes shown. When \code{NULL}, all chromosomes present are used.
#' @param mainChromosomesOnly Logical; when \code{TRUE} (the default) and
#'   \code{chromosomes} is \code{NULL}, only the main chromosomes (numbered,
#'   plus X, Y and MT) are kept and unplaced scaffolds/contigs are dropped.
#' @param ensemblMirror \pkg{biomaRt} mirror to query: one of \code{"www"},
#'   \code{"useast"} or \code{"asia"}. Default \code{"www"}.
#' @param transparency Numeric link transparency in \code{[0, 1]} (0 is opaque).
#'   Default \code{0.5}.
#' @param hRatio Numeric in \code{[0, 1]} controlling how far link ribbons bend
#'   toward the circle centre; smaller values give flatter, less tangled links.
#'   Default \code{0.6}.
#' @param fontFamily Font family used for all plot text, e.g. \code{"sans"}
#'   (Helvetica/Arial, the default) for a publication look.
#' @param legend Logical; whether to draw legends for effect size, novelty and
#'   gene sets. Default \code{TRUE}.
#' @return Invisibly, a list with \code{edges} (the plotted links annotated with
#'   coordinates and novelty) and \code{coords} (the gene coordinate table with
#'   \code{outDegree}, \code{inDegree} and total \code{degree} columns, each the
#'   sum of \code{log2FoldChange} over a gene's outgoing / incoming links). The
#'   function is called for the side effect of drawing the plot.
#' @seealso \code{\link{testEdges}}, \code{\link{runSCORPION}}
#' @details Requires the \pkg{circlize} package, and \pkg{biomaRt} when gene
#'   coordinates are downloaded automatically (\code{geneCoords = NULL}). Genes
#'   without coordinates, and links whose TF or target lacks coordinates, are
#'   dropped with a message.
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
#' # Human coordinates auto-downloaded from Ensembl, known/novel vs a prior net
#' circosEdges(
#'   edgesDF = res,
#'   species = "hsapiens_gene_ensembl",
#'   priorNet = scorpionTest$tf,
#'   geneSets = "hallmark.gmt"
#' )
#' }
#' @export
#' @importFrom cli cli_abort cli_alert_info cli_alert_success cli_alert_warning
#' @importFrom grDevices adjustcolor hcl.colors
#' @importFrom graphics legend
circosEdges <- function(edgesDF,
                        species = "hsapiens_gene_ensembl",
                        geneCoords = NULL,
                        priorNet = NULL,
                        geneSets = NULL,
                        pAdjThreshold = 0.05,
                        log2FCThreshold = 0,
                        maxEdges = 500L,
                        colorBy = "log2FoldChange",
                        linkColors = c("#2166AC", "#F7F7F7", "#B2182B"),
                        lwdRange = c(0.5, 4),
                        nmaxTF = 20L,
                        nmaxTarget = 20L,
                        knownColor = "grey60",
                        novelColor = "#D95F02",
                        geneSetColors = NULL,
                        chromosomes = NULL,
                        mainChromosomesOnly = TRUE,
                        ensemblMirror = "www",
                        transparency = 0.5,
                        hRatio = 0.6,
                        fontFamily = "sans",
                        legend = TRUE) {

  if (!requireNamespace("circlize", quietly = TRUE)) {
    cli::cli_abort(c(
      "Package {.pkg circlize} is required to draw Circos plots.",
      "i" = "Install it with {.code install.packages('circlize')}."
    ))
  }

  # Validate input columns
  if (!is.data.frame(edgesDF)) {
    cli::cli_abort("{.arg edgesDF} must be a data.frame from {.fn testEdges}.")
  }
  # Single-sample testEdges output has no log2FoldChange; fall back to meanEdge
  if (missing(colorBy) && !("log2FoldChange" %in% colnames(edgesDF)) &&
      "meanEdge" %in% colnames(edgesDF)) {
    colorBy <- "meanEdge"
    cli::cli_alert_info("No {.field log2FoldChange} column; colouring by {.field meanEdge}.")
  }
  required_cols <- c("tf", "target")
  if (!identical(colorBy, "log10Padj")) required_cols <- c(required_cols, colorBy)
  missing_cols <- setdiff(required_cols, colnames(edgesDF))
  if (length(missing_cols)) {
    cli::cli_abort("{.arg edgesDF} is missing required column{?s}: {.val {missing_cols}}")
  }
  if (length(linkColors) != 3) {
    cli::cli_abort("{.arg linkColors} must be a vector of exactly 3 colours.")
  }
  if (transparency < 0 || transparency > 1) {
    cli::cli_abort("{.arg transparency} must be between 0 and 1.")
  }

  # Significance column (pAdj preferred, pValue fallback)
  sig_col <- if ("pAdj" %in% colnames(edgesDF)) "pAdj" else if ("pValue" %in% colnames(edgesDF)) "pValue" else NA_character_

  edges <- edgesDF
  edges$tf <- as.character(edges$tf)
  edges$target <- as.character(edges$target)

  # Filter by significance and effect size
  if (!is.na(sig_col)) {
    keep <- !is.na(edges[[sig_col]]) & edges[[sig_col]] < pAdjThreshold
    edges <- edges[keep, , drop = FALSE]
  }
  # Keep every significant edge for degree computation (before fold-change
  # filtering and the maxEdges cap, which only affect which links are drawn)
  edges_sig <- edges
  if (log2FCThreshold > 0 && "log2FoldChange" %in% colnames(edges)) {
    keep <- !is.na(edges$log2FoldChange) & abs(edges$log2FoldChange) >= log2FCThreshold
    edges <- edges[keep, , drop = FALSE]
  }
  if (identical(colorBy, "log10Padj")) {
    if (is.na(sig_col)) {
      cli::cli_abort("{.arg colorBy = 'log10Padj'} requires a {.field pAdj} or {.field pValue} column.")
    }
    edges$log10Padj <- -log10(pmax(edges[[sig_col]], .Machine$double.xmin))
  }
  edges <- edges[!is.na(edges[[colorBy]]), , drop = FALSE]
  if (nrow(edges) == 0) {
    cli::cli_abort("No edges pass the significance / fold-change thresholds.")
  }

  # Cap to the most significant edges
  if (!is.na(sig_col) && nrow(edges) > maxEdges) {
    edges <- edges[order(edges[[sig_col]]), , drop = FALSE][seq_len(maxEdges), , drop = FALSE]
    cli::cli_alert_info("Showing top {maxEdges} edges by {sig_col}.")
  }

  genes <- unique(c(edges$tf, edges$target))

  # Resolve gene coordinates
  if (is.null(geneCoords)) {
    coords <- .fetchGeneCoords(genes, species, ensemblMirror)
  } else {
    coords <- .validateGeneCoords(geneCoords)
  }
  coords <- coords[coords$gene %in% genes, , drop = FALSE]
  coords <- coords[!duplicated(coords$gene), , drop = FALSE]
  if (nrow(coords) == 0) {
    cli::cli_abort("None of the edge genes could be mapped to coordinates.")
  }

  # Restrict to main chromosomes (drop unplaced scaffolds/contigs)
  if (isTRUE(mainChromosomesOnly) && is.null(chromosomes)) {
    bare <- sub("^chr", "", coords$chr, ignore.case = TRUE)
    is_main <- grepl("^[0-9]+$", bare) | toupper(bare) %in% c("X", "Y", "MT", "M")
    dropped <- unique(coords$chr[!is_main])
    coords <- coords[is_main, , drop = FALSE]
    if (length(dropped)) {
      cli::cli_alert_info("Excluding {length(dropped)} non-main sequence{?s} (scaffolds/contigs).")
    }
  }

  # Optionally restrict chromosomes
  if (!is.null(chromosomes)) {
    coords <- coords[coords$chr %in% as.character(chromosomes), , drop = FALSE]
  }
  if (nrow(coords) == 0) {
    cli::cli_abort("No genes remain after chromosome filtering.")
  }
  coord_idx <- match(genes, coords$gene)
  names(coord_idx) <- genes

  # Keep only links whose endpoints both have coordinates
  tf_i <- match(edges$tf, coords$gene)
  tg_i <- match(edges$target, coords$gene)
  mappable <- !is.na(tf_i) & !is.na(tg_i)
  if (!all(mappable)) {
    cli::cli_alert_warning("Dropping {sum(!mappable)} link{?s} with unmapped TF or target.")
    edges <- edges[mappable, , drop = FALSE]
    tf_i <- tf_i[mappable]
    tg_i <- tg_i[mappable]
  }
  if (nrow(edges) == 0) {
    cli::cli_abort("No links remain after coordinate mapping.")
  }

  # Flag known vs novel against the a priori network
  novelty <- rep("novel", nrow(edges))
  if (!is.null(priorNet)) {
    prior <- as.data.frame(priorNet, stringsAsFactors = FALSE)
    if (ncol(prior) < 2) {
      cli::cli_abort("{.arg priorNet} must have at least two columns (TF, target).")
    }
    prior_key <- paste(as.character(prior[[1]]), as.character(prior[[2]]), sep = "\r")
    edge_key <- paste(edges$tf, edges$target, sep = "\r")
    novelty[edge_key %in% prior_key] <- "known"
  }
  edges$novelty <- novelty

  # Per-gene degree: sum of log2FoldChange over outgoing (TF) and incoming
  # (target) links, using every edge under pAdjThreshold (not just drawn links)
  deg_col <- if ("log2FoldChange" %in% colnames(edges_sig)) "log2FoldChange" else colorBy
  out_deg <- tapply(edges_sig[[deg_col]], edges_sig$tf, sum, na.rm = TRUE)
  in_deg <- tapply(edges_sig[[deg_col]], edges_sig$target, sum, na.rm = TRUE)
  coords$outDegree <- as.numeric(out_deg[coords$gene])
  coords$outDegree[is.na(coords$outDegree)] <- 0
  coords$inDegree <- as.numeric(in_deg[coords$gene])
  coords$inDegree[is.na(coords$inDegree)] <- 0
  coords$degree <- coords$outDegree + coords$inDegree

  # Build chromosome layout (BED-like ranges spanning each chromosome)
  chr_order <- if (!is.null(chromosomes)) as.character(chromosomes) else .orderChr(coords$chr)
  chr_order <- chr_order[chr_order %in% coords$chr]
  chr_layout <- do.call(rbind, lapply(chr_order, function(ch) {
    s <- coords[coords$chr == ch, , drop = FALSE]
    rng <- range(c(s$start, s$end), na.rm = TRUE)
    data.frame(chr = ch, start = rng[1], end = rng[2], stringsAsFactors = FALSE)
  }))

  # Link fill colours: diverging (signed values) or sequential (all non-negative)
  vals <- edges[[colorBy]]
  diverging <- any(vals < 0, na.rm = TRUE)
  if (diverging) {
    maxabs <- max(abs(vals), na.rm = TRUE)
    if (!is.finite(maxabs) || maxabs == 0) maxabs <- 1
    value_range <- c(-maxabs, maxabs)
    col_fun <- circlize::colorRamp2(c(-maxabs, 0, maxabs), linkColors)
  } else {
    rng <- range(vals, na.rm = TRUE)
    if (!all(is.finite(rng)) || diff(rng) == 0) rng <- c(0, max(1, rng[2]))
    value_range <- rng
    col_fun <- circlize::colorRamp2(c(rng[1], mean(rng), rng[2]), linkColors)
  }
  fill_cols <- grDevices::adjustcolor(col_fun(vals), alpha.f = 1 - transparency)
  border_cols <- ifelse(edges$novelty == "known", knownColor, novelColor)

  # Link thickness scaled linearly by -log10 adjusted p-value
  if (!is.na(sig_col)) {
    sig_vals <- -log10(pmax(edges[[sig_col]], .Machine$double.xmin))
    srng <- range(sig_vals, na.rm = TRUE)
    if (!all(is.finite(srng)) || diff(srng) == 0) {
      link_lwd <- rep(mean(lwdRange), nrow(edges))
    } else {
      link_lwd <- lwdRange[1] + (sig_vals - srng[1]) / (srng[2] - srng[1]) * (lwdRange[2] - lwdRange[1])
    }
  } else {
    link_lwd <- rep(lwdRange[1], nrow(edges))
  }

  region1 <- data.frame(
    chr = coords$chr[tf_i], start = coords$start[tf_i], end = coords$start[tf_i],
    stringsAsFactors = FALSE
  )
  region2 <- data.frame(
    chr = coords$chr[tg_i], start = coords$start[tg_i], end = coords$start[tg_i],
    stringsAsFactors = FALSE
  )

  # Gene-set colour assignment (optional)
  set_colors <- NULL
  gene_set_list <- NULL
  if (!is.null(geneSets)) {
    sets <- if (is.character(geneSets) && length(geneSets) == 1) .parseGMT(geneSets) else geneSets
    if (!is.list(sets) || is.null(names(sets))) {
      cli::cli_abort("{.arg geneSets} must be a GMT file path or a named list of gene vectors.")
    }
    gene_set_list <- sets
    set_names <- names(sets)
    if (is.null(geneSetColors)) {
      set_colors <- stats::setNames(grDevices::hcl.colors(length(set_names), "Dark3"), set_names)
    } else {
      set_colors <- geneSetColors
    }
  }

  # Select the top TFs and targets to label by |out-degree| / |in-degree|,
  # styling TFs in bold (font 2) and targets in italics (font 3).
  tf_genes <- unique(edges$tf)
  target_genes <- unique(edges$target)
  pick_top <- function(cands, score_col, n) {
    if (!length(cands)) return(character(0))
    if (is.null(n) || !is.finite(n)) return(cands)
    sc <- abs(coords[[score_col]][match(cands, coords$gene)])
    cands[order(-sc)][seq_len(min(as.integer(n), length(cands)))]
  }
  label_genes <- union(pick_top(tf_genes, "outDegree", nmaxTF),
                       pick_top(target_genes, "inDegree", nmaxTarget))
  lab_idx <- which(coords$gene %in% label_genes)
  label_font <- ifelse(coords$gene[lab_idx] %in% tf_genes, 2L, 3L)
  # Colour TF labels by out-degree sign: blue if negative, red if positive
  # (fully saturated endpoints, no fade). Targets are black.
  # Cap the ring colour scale at the 90th percentile of non-zero degrees so
  # typical genes show saturated colour instead of washing out near the max.
  deg_vals <- c(coords$outDegree, coords$inDegree)
  nz <- abs(deg_vals[deg_vals != 0 & is.finite(deg_vals)])
  deg_absmax <- if (length(nz)) as.numeric(stats::quantile(nz, 0.9, names = FALSE)) else 0
  if (!is.finite(deg_absmax) || deg_absmax == 0) {
    deg_absmax <- max(abs(deg_vals), na.rm = TRUE)
  }
  if (!is.finite(deg_absmax) || deg_absmax == 0) deg_absmax <- 1
  deg_col_fun <- circlize::colorRamp2(c(-deg_absmax, 0, deg_absmax), linkColors)
  is_tf_lab <- coords$gene[lab_idx] %in% tf_genes
  label_col <- rep("black", length(lab_idx))
  label_col[is_tf_lab] <- ifelse(
    coords$outDegree[lab_idx][is_tf_lab] < 0, linkColors[1], linkColors[3]
  )
  label_bed <- data.frame(
    chr = coords$chr[lab_idx], start = coords$start[lab_idx], end = coords$start[lab_idx],
    gene = coords$gene[lab_idx], font = label_font, col = label_col,
    stringsAsFactors = FALSE
  )

  # Draw
  op <- graphics::par(family = fontFamily, xpd = NA, mar = c(0, 0, 0, 0))
  on.exit(graphics::par(op), add = TRUE)
  circlize::circos.clear()
  on.exit(circlize::circos.clear(), add = TRUE)
  n_sectors <- nrow(chr_layout)
  gaps <- rep(1.5, n_sectors)
  # Seam gap widens with the longest track label so names always fit
  seam_labels <- c("in-degree", "out-degree", names(gene_set_list))
  gaps[n_sectors] <- min(90, max(15, max(nchar(seam_labels)) * 2.1))
  circlize::circos.par(
    cell.padding = c(0, 0, 0, 0), points.overflow.warning = FALSE,
    gap.after = gaps, start.degree = 90, track.margin = c(0.006, 0.004),
    canvas.xlim = c(-1.2, 1.2), canvas.ylim = c(-1.2, 1.2)
  )
  circlize::circos.genomicInitialize(chr_layout, plotType = NULL)

  # Outermost ring: gene / TF labels
  if (!is.null(label_bed) && nrow(label_bed) > 0) {
    circlize::circos.genomicLabels(
      label_bed[, c("chr", "start", "end", "gene")],
      labels.column = 4,
      col = label_bed$col,
      line_col = label_bed$col,
      font = label_bed$font,
      side = "outside",
      cex = 0.5
    )
  }

  # Shared tile geometry for the degree rings (deg_col_fun defined above)
  span_max <- max(chr_layout$end - chr_layout$start, na.rm = TRUE)
  tile_w <- span_max / 300

  ring_bed <- data.frame(
    chr = coords$chr,
    start = pmax(0, coords$start - tile_w), end = coords$start + tile_w,
    stringsAsFactors = FALSE
  )
  seam_label <- function(name, col = "grey30") {
    if (circlize::CELL_META$sector.numeric.index == 1) {
      circlize::circos.text(
        circlize::CELL_META$cell.xlim[1], 0.5, labels = name,
        facing = "downward", adj = c(1, 0.5), cex = 0.5, col = col
      )
    }
  }

  # Thin heatmap ring of a per-gene value on the diverging degree scale
  circlize::circos.par(track.margin = c(0.006, 0.004))
  deg_ring <- function(v, name) {
    circlize::circos.genomicTrack(
      cbind(ring_bed, value = v), ylim = c(0, 1), track.height = 0.03,
      bg.border = "grey85", bg.lwd = 0.5,
      panel.fun = function(region, value, ...) {
        # Draw weakest first so the strongest tiles sit on top when they overlap
        ord <- order(abs(value[[1]]), decreasing = FALSE, na.last = TRUE)
        circlize::circos.genomicRect(
          region[ord, , drop = FALSE], value[ord, , drop = FALSE],
          ytop = 1, ybottom = 0, border = NA,
          col = deg_col_fun(value[[1]][ord])
        )
        seam_label(name)
      }
    )
  }

  # In-degree ring: always the outermost data track, just inside the labels
  deg_ring(coords$inDegree, "in-degree")

  # One thin presence ring per gene set (member genes coloured by set),
  # drawn just inside the in-degree ring
  if (!is.null(gene_set_list) && length(gene_set_list)) {
    for (sn in names(gene_set_list)) {
      scol <- if (!is.null(set_colors)) unname(set_colors[[sn]]) else "grey40"
      member <- as.integer(coords$gene %in% gene_set_list[[sn]])
      local({
        col_s <- scol
        nm <- sn
        circlize::circos.genomicTrack(
          cbind(ring_bed, value = member), ylim = c(0, 1), track.height = 0.022,
          bg.border = "grey85", bg.lwd = 0.5,
          panel.fun = function(region, value, ...) {
            keep <- value[[1]] == 1
            if (any(keep)) {
              circlize::circos.genomicRect(
                region[keep, , drop = FALSE], value[keep, , drop = FALSE],
                ytop = 1, ybottom = 0, col = col_s, border = NA
              )
            }
            seam_label(nm)
          }
        )
      })
    }
  }

  # Chromosome name band, then the out-degree ring innermost
  circlize::circos.track(
    ylim = c(0, 1), track.height = 0.045, bg.col = "grey92", bg.border = NA,
    panel.fun = function(x, y) {
      circlize::circos.text(
        circlize::CELL_META$xcenter, circlize::CELL_META$ycenter,
        labels = circlize::CELL_META$sector.index,
        cex = 0.7, font = 2, facing = "inside", niceFacing = TRUE
      )
    }
  )
  deg_ring(coords$outDegree, "out-degree")

  # Draw strongest links last so the signal sits on top of faint edges
  draw_order <- order(abs(vals), decreasing = FALSE, na.last = TRUE)
  circlize::circos.genomicLink(
    region1[draw_order, , drop = FALSE], region2[draw_order, , drop = FALSE],
    col = fill_cols[draw_order], border = border_cols[draw_order],
    lwd = link_lwd[draw_order], h.ratio = hRatio
  )

  if (isTRUE(legend)) {
    present_sets <- names(gene_set_list)
    color_label <- if (identical(colorBy, "log10Padj")) "-log10(padj)" else colorBy
    .drawCircosLegends(
      col_fun = col_fun, valueRange = value_range, colorBy = color_label,
      hasPrior = !is.null(priorNet), knownColor = knownColor, novelColor = novelColor,
      set_colors = if (length(present_sets)) set_colors[present_sets] else NULL,
      lwdRange = lwdRange, widthShown = !is.na(sig_col),
      degCol_fun = deg_col_fun, degRange = c(-deg_absmax, deg_absmax),
      degColBy = deg_col
    )
  }

  cli::cli_alert_success("Drew {nrow(edges)} link{?s} across {nrow(chr_layout)} chromosome{?s}.")
  invisible(list(edges = edges, coords = coords))
}

#' Parse a GMT gene-set file into a named list.
#' @keywords internal
.parseGMT <- function(path) {
  if (!file.exists(path)) {
    cli::cli_abort("GMT file not found: {.path {path}}")
  }
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  sets <- lapply(lines, function(l) {
    parts <- strsplit(l, "\t", fixed = TRUE)[[1]]
    genes <- if (length(parts) > 2) parts[-(1:2)] else character(0)
    genes[nzchar(genes)]
  })
  names(sets) <- vapply(lines, function(l) strsplit(l, "\t", fixed = TRUE)[[1]][1], character(1))
  sets
}

#' Validate a user-supplied gene coordinate table.
#' @keywords internal
.validateGeneCoords <- function(geneCoords) {
  if (!is.data.frame(geneCoords)) {
    cli::cli_abort("{.arg geneCoords} must be a data.frame.")
  }
  needed <- c("gene", "chr", "start", "end")
  missing_cols <- setdiff(needed, colnames(geneCoords))
  if (length(missing_cols)) {
    cli::cli_abort("{.arg geneCoords} is missing column{?s}: {.val {missing_cols}}")
  }
  data.frame(
    gene = as.character(geneCoords$gene),
    chr = as.character(geneCoords$chr),
    start = as.numeric(geneCoords$start),
    end = as.numeric(geneCoords$end),
    stringsAsFactors = FALSE
  )
}

#' Download gene coordinates from Ensembl via biomaRt.
#' @keywords internal
.fetchGeneCoords <- function(genes, species, mirror) {
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    cli::cli_abort(c(
      "Package {.pkg biomaRt} is required to download gene coordinates.",
      "i" = "Install it with {.code BiocManager::install('biomaRt')}, or pass {.arg geneCoords}."
    ))
  }
  cli::cli_alert_info("Querying Ensembl ({species}) for {length(genes)} genes via biomaRt.")
  mart <- biomaRt::useEnsembl(biomart = "genes", dataset = species, mirror = mirror)
  attrs <- c("external_gene_name", "chromosome_name", "start_position", "end_position")
  bm <- biomaRt::getBM(
    attributes = attrs, filters = "external_gene_name",
    values = genes, mart = mart
  )
  bm <- bm[!is.na(bm$start_position) & nzchar(bm$chromosome_name), , drop = FALSE]
  data.frame(
    gene = as.character(bm$external_gene_name),
    chr = as.character(bm$chromosome_name),
    start = as.numeric(bm$start_position),
    end = as.numeric(bm$end_position),
    stringsAsFactors = FALSE
  )
}

#' Natural-sort chromosome names (1..N, then X, Y, MT, then the rest).
#' @keywords internal
.orderChr <- function(chr) {
  u <- unique(as.character(chr))
  bare <- sub("^chr", "", u, ignore.case = TRUE)
  num <- suppressWarnings(as.numeric(bare))
  special <- match(toupper(bare), c("X", "Y", "MT", "M"))
  # numeric first (by value), then X/Y/MT, then everything else alphabetically
  key1 <- ifelse(!is.na(num), 1L, ifelse(!is.na(special), 2L, 3L))
  key2 <- ifelse(!is.na(num), num, ifelse(!is.na(special), special, 0))
  u[order(key1, key2, u)]
}

#' Draw effect-size, novelty and gene-set legends on a Circos plot.
#' @keywords internal
.drawCircosLegends <- function(col_fun, valueRange, colorBy, hasPrior,
                               knownColor, novelColor, set_colors,
                               lwdRange = NULL, widthShown = FALSE,
                               degCol_fun = NULL, degRange = NULL,
                               degColBy = "log2FoldChange") {
  # Render known measures as tidy plotmath; fall back to the raw string
  nice_label <- function(txt) {
    if (!is.character(txt)) return(txt)
    switch(txt,
      "log2FoldChange" = expression(log[2] ~ "fold change"),
      "-log10(padj)"   = expression(-log[10] * "(" * p[adj] * ")"),
      txt
    )
  }

  y_top <- -1.05     # common top edge shared by every legend title
  gap   <- 0.08      # horizontal gap between blocks
  barw  <- 0.45      # colour-bar length
  bw    <- 0.03      # colour-bar thickness (thin)
  # Measure text so colour bars line up with the legend swatch rows
  th     <- graphics::strheight("Ag", cex = 0.68, font = 2)  # title height
  kh     <- graphics::strheight("Ag", cex = 0.60)            # key-row height
  bar_cy <- y_top - th - 0.5 * kh - 0.015                    # first-key centre

  # A continuous colour-bar block, self-centred within its measured width
  bar_block <- function(cfun, rng, title) {
    tw <- graphics::strwidth(title, cex = 0.68, font = 2)
    w  <- max(barw, tw)
    list(w = w, draw = function(xl) {
      cx  <- xl + w / 2
      bxl <- cx - barw / 2
      bxr <- cx + barw / 2
      ytt <- bar_cy + bw / 2
      ybb <- bar_cy - bw / 2
      n <- 128L
      xs <- seq(bxl, bxr, length.out = n + 1L)
      gvals <- seq(rng[1], rng[2], length.out = n)
      graphics::rect(xs[-(n + 1L)], ybb, xs[-1L], ytt,
                     col = cfun(gvals), border = NA, xpd = NA)
      graphics::rect(bxl, ybb, bxr, ytt, border = "grey40", lwd = 0.6, xpd = NA)
      breaks <- pretty(rng, n = 4)
      breaks <- breaks[breaks >= rng[1] & breaks <= rng[2]]
      tick_x <- bxl + (breaks - rng[1]) / diff(rng) * (bxr - bxl)
      graphics::segments(tick_x, ybb, tick_x, ybb - 0.018, lwd = 0.6, xpd = NA)
      graphics::text(tick_x, ybb - 0.04, labels = format(breaks, digits = 2),
                     adj = c(0.5, 1), cex = 0.5, xpd = NA)
      graphics::text(cx, y_top, labels = title, adj = c(0.5, 1),
                     cex = 0.68, font = 2, xpd = NA)
    })
  }

  # A discrete graphics::legend block: measure its width, then draw at an anchor
  legend_block <- function(...) {
    m <- graphics::legend(x = 0, y = 0, ..., plot = FALSE)
    list(w = m$rect$w, draw = function(xl)
      graphics::legend(x = xl, y = y_top, ..., xjust = 0, yjust = 1, xpd = NA))
  }

  # Assemble only the blocks that apply to this plot. The gene-set ring seam
  # labels already name each set, so no separate gene-set legend is drawn.
  blocks <- list(bar_block(col_fun, valueRange, nice_label(colorBy)))
  if (!is.null(degCol_fun) && !is.null(degRange)) {
    deg_title <- if (identical(degColBy, "log2FoldChange")) {
      expression("degree" ~ (Sigma ~ log[2] * "FC"))
    } else {
      bquote("degree" ~ (Sigma ~ .(degColBy)))
    }
    blocks <- c(blocks, list(bar_block(degCol_fun, degRange, deg_title)))
  }
  if (isTRUE(hasPrior)) {
    blocks <- c(blocks, list(legend_block(
      title = "Link", legend = c("known", "novel"),
      col = c(knownColor, novelColor), lwd = 3, bty = "n", cex = 0.6,
      y.intersp = 0.5)))
  }
  if (isTRUE(widthShown) && !is.null(lwdRange)) {
    blocks <- c(blocks, list(legend_block(
      title = expression(-log[10] * "(" * p[adj] * ")"),
      legend = c("low", "high"), lwd = lwdRange, col = "grey30",
      bty = "n", cex = 0.6, y.intersp = 0.5)))
  }

  # Pack the blocks left-to-right and centre the whole row under the circle
  widths <- vapply(blocks, function(b) b$w, numeric(1))
  x <- -(sum(widths) + gap * (length(blocks) - 1)) / 2
  for (b in blocks) {
    b$draw(x)
    x <- x + b$w + gap
  }
}
