#' @importFrom methods is
buildKNND <- function(D, k = 5, return_neighbors_order = T, mode = 'all'){

  ##print("Start knn_graph_from_dist")
  if(!is(D, 'matrix') | !is(D, 'dist')){
    stop("D (X) mast be a matrix or dist!")
  }

  if(!is(D, 'dist')){
    D <- as.dist(D)
  }



  N        <- (1 + sqrt(1+8*length(D)))/2 # number of cells

  if (k >= N)
    stop("Not enought neighbors in data set!")
  if (k < 1)
    stop("Invalid number of nearest neighbors, k must be >= 1!")

  row <- function(i, N){
    return(c(if(i>1) D[(i-1)+c(0:(i-2))*(N - 1 - c(1:(i-1))/2)],
             NA,
             if(i < N) D[((i-1)*(N-1) - ((i-1)*(i-2)/2) + 1) : (((i-1)*(N-1) - ((i-1)*(i-2)/2) + 1) + N-i-1)]))
  }

  neighbors <- t(sapply(1:N, function(i) {order(row(i,N))[1:k]}))

  adj.knn <- split(neighbors, rep(1:nrow(neighbors), times = ncol(neighbors)))


  graph.knn     <- igraph::graph_from_adj_list(adj.knn,  duplicate = F, mode = mode)
  graph.knn     <- igraph::simplify(graph.knn, remove.multiple = T)
  igraph::E(graph.knn)$weight <- 1

  if(return_neighbors_order){
    res <- list(graph.knn = graph.knn,
                order = neighbors)
  } else {res <- list(graph.knn = graph.knn)}

  return(res)
}
