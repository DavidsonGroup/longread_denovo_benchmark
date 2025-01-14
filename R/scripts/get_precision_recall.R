# x is the true set, y is the predicted set
get_precision_recall <- function(x,y) {
  tp <- length(intersect(x,y))
  fn <- length(setdiff(x,y))
  fp <- length(setdiff(y,x))
  Precision = tp/(tp+fp)
  Recall = tp/(tp+fn)
  F1 = 2*Precision*Recall/(Precision + Recall)
  res <- c(tp, fn, fp, Precision, Recall, F1)
  names(res) <- c('tp','fn','fp', 'Precision', 'Recall', 'F1')
  return(res)
}