#' Title
#'
#' @param k
#' @param base
#' @param sigma
#'
#' @return
#' @export
#'
#' @examples
weighted.k.circle.projection.members <- function(k, base = 10, sigma = 14){

  seq <- seq_along(k)
  m <- level.membership(k, "two-mode", levels = seq) %>% filter(type == 2)

  i <- k[[1]]
  ind <- k[[1]]

  stopifnot(all.equal(colnames(i), m$Name))
  il <- t(i) * (m$Level - 1)
  il <- t(il)
  rs <- rowSums(il)

  # Occasions weight
  Y                       <- slam::as.simple_triplet_matrix(drop0(Matrix::t(ind)))
  col.max                 <- as(slam::rollup(Y, 2, FUN = max), "sparseVector")
  col.max                 <- as.numeric(col.max)
  ind                    <- Matrix::t(Matrix::t(ind) * (1 / col.max))

  # Size weight
  affil.members           <- Matrix::colSums(ind)
  affil.memberships       <- Matrix::rowSums(ind)
  affil.weight                    <- sqrt((sigma/affil.members))
  affil.weight[affil.weight > 1]  <- 1

  ind              <- Matrix::t(ind) * affil.weight
  ind              <- Matrix::t(ind)
  weighted.memberships <- rowSums(ind)

  ind              <- t(ind) * (m$Level - 1)
  ind              <- t(ind)

  adj.ind                <- Matrix::crossprod(t(sqrt(ind)))
  affil.diag             <- Matrix::diag(adj.ind) - weighted.memberships

  a <- tibble(NAME = rownames(ind) ,"Number of memberships" = affil.memberships, "Weighted memberships" = weighted.memberships, "Weighted k-sum" = affil.diag |> round(3), "k-sum" = rs)
  a <- a |> mutate(rank.dif = rank(`Weighted k-sum`) - rank(`k-sum`) )

  adj.ind@x <- 1/adj.ind@x
  g <- graph_from_adjacency_matrix(adj.ind, weighted = TRUE, mode = "undirected")
  g <- delete.edges(g, which(E(g)$weight == Inf))

  list(scores = a, graph = g, ind = ind)
}

#' Title
#'
#' @param k
#' @param base
#' @param sigma
#'
#' @return
#' @export
#'
#' @examples
weighted.k.circle.projection.affiliations <- function(k, base = 10, sigma = 14){

  seq <- seq_along(k)
  m <- level.membership(k, "two-mode", levels = seq) %>% filter(type == 1)
  i <- k[[1]]
  ind <- i

  stopifnot(all.equal(rownames(i), m$Name))

  il <- i * (m$Level - 1)

  # Occasions weight
  Y                       <- slam::as.simple_triplet_matrix(drop0(Matrix::t(ind)))
  col.max                 <- as(slam::rollup(Y, 2, FUN = max), "sparseVector")
  col.max                 <- as.numeric(col.max)
  ind                    <- Matrix::t(Matrix::t(ind) * (1 / col.max))

  # Size weight
  affil.members           <- Matrix::colSums(ind)
  affil.weight                    <- sqrt((sigma/affil.members))
  affil.weight[affil.weight > 1]  <- 1

  ind              <- Matrix::t(ind) * affil.weight
  ind              <- Matrix::t(ind)
  weighted.members <- colSums(ind)
  ind              <- ind * (m$Level - 1)

  adj.affil              <- Matrix::crossprod(sqrt(ind))
  affil.diag             <- Matrix::diag(adj.affil) - weighted.members

  cs <- colSums(il) - colSums(i)

  a <- tibble(AFFILIATION = colnames(ind) ,"Number of members" = affil.members, "Affiliation weight" = affil.weight, "Weighted k-sum" = affil.diag |> round(3), "k-sum" = cs)
  a <- a |> mutate(rank.dif = rank(`Weighted k-sum`) - rank(`k-sum`) )

  adj.affil@x <- 1/adj.affil@x
  g <- graph_from_adjacency_matrix(adj.affil, weighted = TRUE, mode = "directed", diag = FALSE)
  g <- delete.edges(g, which(E(g)$weight == Inf))


  list(scores = a, graph = g, ind = ind)
}

#' Title
#'
#' @param ind
#' @param base
#' @param sigma
#'
#' @return
#' @export
#'
#' @examples
size.weight.incidence <- function(ind, base = 10, sigma = 14){

  # Occasions weight
  Y                       <- slam::as.simple_triplet_matrix(drop0(Matrix::t(ind)))
  col.max                 <- as(slam::rollup(Y, 2, FUN = max), "sparseVector")
  col.max                 <- as.numeric(col.max)
  ind                    <- Matrix::t(Matrix::t(ind) * (1 / col.max))

  # Size weight
  affil.members           <- Matrix::colSums(ind)
  affil.weight                    <- sqrt((sigma/affil.members))
  affil.weight[affil.weight > 1]  <- 1

  ind              <- Matrix::t(ind) * affil.weight
  ind              <- Matrix::t(ind)
  ind
}
