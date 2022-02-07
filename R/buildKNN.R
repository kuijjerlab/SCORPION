buildKNN <- function(X,
                            k = 5,
                            from = c("dist", "coordinates"),
                            use.nn2 = TRUE,
                            return_neighbors_order = F,
                            dist_method = "euclidean",
                            cor_method = "pearson",
                            p = 2,
                            directed = FALSE)
{
  av.methods <- c("dist", "coordinates")
  method <-  pmatch(from[1], av.methods)
  if (is.na(method)) {
    stop(paste(
      "Unlnown method",
      from,
      "Available methods are",
      paste(av.methods, collapse = ", ")
    ))
  }


  if (method == 2) {
    # from coordinates
    if (use.nn2) {
      if (dist_method != "euclidean") {
        stop(
          paste0(
            "Fast nn2 function from RANN package is used, so",
            dist_method,
            "distnce is not acceptable.
            To use nn2 method, please, choose eucleadian distance.
            If you want to use",
            dist_method,
            "distance, please set parameter use.nn2 to FALSE"
          )
        )
      }
      mode <- ifelse(directed, 'out', 'all')
      return(buildNN2(X = X, k = k, mode = mode))
    } else {
      av.dist      <-
        c("cor",
          "euclidean",
          "maximum",
          "manhattan",
          "canberra",
          "binary",
          "minkowski")
      dist_method_ <-  pmatch(dist_method, av.dist)
      if (is.na(dist_method_)) {
        stop(
          paste(
            "Unknown distance method:",
            dist_method,
            "Available dist methods are",
            paste(av.dist, collapse = ", ")
          )
        )
      }
      if (dist_method_ == 1) {
        #print("cor")
        #print(cor_method)
        av.cor_methods <- c("pearson", "kendall", "spearman")
        cor_method_    <- pmatch(cor_method, av.cor_methods)
        if (is.na(cor_method_)) {
          stop(
            paste(
              "Unknown cor method:",
              cor_method,
              "Available cor methods are",
              paste(av.cor_methods)
            )
          )
        }
        X <- as.dist(as.matrix(1 - cor(t(X), method = cor_method)))
      } else {
        X <- dist(X, method = dist_method)
      }
    }
  } else {
    if (use.nn2 == TRUE) {
      stop(
        "Method nn2 cannot be applied to distance, to use fast nn2 method, please provide coordinates rather than distance
        and set parameter from to coordinates"
      )
    }
    return(knn_graph_from_dist(
      D = X,
      k = k,
      return_neighbors_order = return_neighbors_order
    ))
    }


  ### now X is distance in any case
  return(knn_graph_from_dist(
    D = X,
    k = k,
    return_neighbors_order = return_neighbors_order
  ))

  }
