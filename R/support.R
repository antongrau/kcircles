#' Title
#'
#' @param incidence
#'
#' @return
#' @export
#'
#' @examples
largest.component.in.incidence <- function(incidence){

  g <- graph_from_incidence_matrix(incidence)
  g <- largest.component(g)
  incidence <- incidence[rownames(incidence) %in% V(g)$name, colnames(incidence) %in% V(g)$name]
  incidence
}

level.summary <- function(l.inc){
  l.inc  <- compact(l.inc)
  l.g    <- map(l.inc, ~graph_from_incidence_matrix(incidence = .x))
  l.cl   <- map(l.g, clusters)

  map_dbl(l.cl, "no")
  l.g %>% map(~bipartite.projection(.x)[[1]]) %>% map(degree) %>% map_dbl(mean)

}

print.k.circles <- function(x){
  x
}

#' Title
#'
#' @param x
#'
#' @return
#' @export cross_sum
#'
#' @examples
#'
cross_sum <- function(x) {
  M1 <- x %*% Matrix::t((x != 0)) + Matrix::t(x %*% Matrix::t(x != 0))
  Matrix::diag(M1) <- Matrix::diag(M1)/2
  M1
}
