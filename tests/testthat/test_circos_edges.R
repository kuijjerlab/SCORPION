test_that("circosEdges draws with user-supplied coordinates (offline)", {
  skip_on_cran()
  skip_if_not_installed("circlize")

  set.seed(1)
  n <- 30
  edges <- data.frame(
    tf = paste0("TF", sample(1:5, n, replace = TRUE)),
    target = paste0("G", sample(1:20, n, replace = TRUE)),
    log2FoldChange = rnorm(n),
    tStatistic = rnorm(n),
    pValue = runif(n, 0, 0.01),
    pAdj = runif(n, 0, 0.04),
    stringsAsFactors = FALSE
  )
  edges <- edges[edges$tf != edges$target, ]

  all_genes <- unique(c(edges$tf, edges$target))
  coords <- data.frame(
    gene = all_genes,
    chr = sample(c("1", "2", "X"), length(all_genes), replace = TRUE),
    start = sample(1:1e6, length(all_genes)),
    end = sample(1:1e6, length(all_genes)) + 1000,
    stringsAsFactors = FALSE
  )

  prior <- data.frame(tf = edges$tf[1:5], target = edges$target[1:5])
  gene_sets <- list(SetA = all_genes[1:5], SetB = all_genes[6:10])

  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp)
  on.exit({
    grDevices::dev.off()
    unlink(tmp)
  }, add = TRUE)

  out <- circosEdges(
    edgesDF = edges,
    geneCoords = coords,
    priorNet = prior,
    geneSets = gene_sets,
    pAdjThreshold = 0.05,
    nmaxTF = 3,
    nmaxTarget = 3
  )

  expect_type(out, "list")
  expect_true(all(c("edges", "coords") %in% names(out)))
  expect_true(all(c("outDegree", "inDegree", "degree") %in% names(out$coords)))
  expect_equal(out$coords$degree, out$coords$outDegree + out$coords$inDegree)
  expect_true(all(out$edges$novelty %in% c("known", "novel")))
  expect_true(any(out$edges$novelty == "known"))
  expect_true(all(out$edges[["pAdj"]] < 0.05))
})

test_that("circosEdges parses a GMT file", {
  skip_on_cran()
  gmt <- tempfile(fileext = ".gmt")
  writeLines(c(
    "SET1\tdescription\tGENEA\tGENEB\tGENEC",
    "SET2\thttp://example\tGENEC\tGENED"
  ), gmt)
  on.exit(unlink(gmt), add = TRUE)

  parseGMT <- getFromNamespace(".parseGMT", "SCORPION")
  sets <- parseGMT(gmt)
  expect_named(sets, c("SET1", "SET2"))
  expect_equal(sets$SET1, c("GENEA", "GENEB", "GENEC"))
  expect_equal(sets$SET2, c("GENEC", "GENED"))
})

test_that("circosEdges errors on missing columns and empty selections", {
  skip_on_cran()
  skip_if_not_installed("circlize")

  bad <- data.frame(tf = "A", target = "B")
  expect_error(circosEdges(bad, colorBy = "log2FoldChange"), "log2FoldChange")

  edges <- data.frame(
    tf = "TF1", target = "G1", log2FoldChange = 1,
    pValue = 0.5, pAdj = 0.9, stringsAsFactors = FALSE
  )
  coords <- data.frame(
    gene = c("TF1", "G1"), chr = c("1", "1"),
    start = c(1, 100), end = c(10, 110), stringsAsFactors = FALSE
  )
  expect_error(
    circosEdges(edges, geneCoords = coords, pAdjThreshold = 0.05),
    "thresholds"
  )
})

test_that(".orderChr sorts chromosomes naturally", {
  orderChr <- getFromNamespace(".orderChr", "SCORPION")
  expect_equal(
    orderChr(c("2", "X", "1", "MT", "10", "Y")),
    c("1", "2", "10", "X", "Y", "MT")
  )
})
