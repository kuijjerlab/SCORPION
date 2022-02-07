spreadNet <- function(df) {
  df[, 3] <- as.numeric(df[, 3])
  row_names <- unique(df[, 1])
  col_names <- unique(df[, 2])
  spread.df <-
    data.frame(matrix(0, nrow = length(row_names), ncol = length(col_names)), row.names =
                 row_names)
  colnames(spread.df) <- col_names
  for (i in 1:nrow(df)) {
    spread.df[as.character(df[i, 1]), as.character(df[i, 2])] <-
      df[i, 3]
  }
  spread.df
}
