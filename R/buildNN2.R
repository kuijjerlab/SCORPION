#' @importFrom RANN nn2
#' @importFrom igraph graph_from_adj_list simplify
buildNN2 <-
  function(X,
           k = min(5, ncol(X)),
           mode = "all") {
    nn2.res <- RANN::nn2(data = X, k = k)
    nn2.res <- nn2.res$nn.idx

    adj.knn <-
      split(nn2.res, rep(1:nrow(nn2.res), times = ncol(nn2.res))) # get adj list

    graph.knn <-
      igraph::graph_from_adj_list(adj.knn, duplicate = F, mode = mode)

    graph.knn <- igraph::simplify(graph.knn, remove.multiple = T)
    igraph::E(graph.knn)$weight <- 1

    return(res <- list(graph.knn = graph.knn))
  }
