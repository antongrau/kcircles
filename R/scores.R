#' Calculate the k-circle scores for clusters
#'
#' @param k a list of incidence matrices as produced by \link{k.circles}
#'
#' @return a list of data.frames
#' @export circle.esteem.for.clusters
#'
#' @examples
circle.esteem.for.clusters <- function(k){
  warning("esteem value not adjusted correctly")
  seq <- seq_along(k)

  m <- level.membership(k, "two-mode", levels = seq) %>% filter(type == 1)
  i <- k[[1]]

  stopifnot(all.equal(rownames(i), m$Name))
  il <- i * (m$Level - 1)

  l.il <- Matrix::mat2triplet(il) %>% bind_cols()
  l.il$NAME        <- rownames(il)[l.il$i]
  l.il$AFFILIATION <- colnames(il)[l.il$j]

  #mc   <- attr(last(k), "col.circle.merge.cluster")
  mck  <- map(k, ~attr(.x, "col.circle.merge.cluster")) %>% bind_rows(.id = "level" )
  mck  <- mck %>% group_by(AFFILIATION) %>% filter(level == max(level))
  mc   <- distinct(mck, AFFILIATION, overlap.cluster, .keep_all = TRUE) %>% filter(AFFILIATION != overlap.cluster) %>% group_by(AFFILIATION) %>% mutate(n = n())

  if(any(mc$n > 1)) warning("Some affiliations are in multiple overlap.clusters - this is probably an error in merge.perfect.overlap()")

  l.il <- left_join(l.il, mc, by = c("AFFILIATION"), relationship = "many-to-many")
  s    <- is.na(l.il$overlap.cluster)
  l.il$overlap.cluster[s] <- l.il$AFFILIATION[s]

  m          <- level.membership(k, "two-mode", levels = seq) %>% filter(type == 2)
  by.affil   <- l.il %>% group_by(AFFILIATION) %>% summarise(circle.esteem = sum(x) - n(), n.linkers = sum(x != 0), n.members = n())
  by.affil   <- left_join(by.affil, mc %>% select(-level, -n), by = "AFFILIATION")
  by.affil   <- by.affil |> left_join(m |> select(-type), by = c("AFFILIATION" = "Name"))

  l.il       <- l.il |> left_join(m |> select(-type), by = c("AFFILIATION" = "Name"))
  by.cluster <- l.il %>% group_by(NAME, overlap.cluster) %>% summarise(x = max(x), Level = max(Level)) %>% group_by(AFFILIATION = overlap.cluster) %>% summarise(circle.esteem = sum(x) - n(), n.linkers = sum(x != 0), n.members = n(), Level = max(Level))

  by.cluster$is.cluster <- by.cluster$AFFILIATION %in% mc$overlap.cluster

  o <- list(affil = by.affil, cluster = by.cluster)

  # By individual ----
  m <- level.membership(k, "two-mode", levels = seq) %>% filter(type == 2)
  il <- t(t(i) * (m$Level - 1))

  l.il <- Matrix::mat2triplet(il) %>% bind_cols()
  l.il$NAME        <- rownames(il)[l.il$i]
  l.il$AFFILIATION <- colnames(il)[l.il$j]

  #mc   <- attr(last(k), "col.circle.merge.cluster")
  mck  <- map(k, ~attr(.x, "col.circle.merge.cluster")) %>% bind_rows(.id = "level" )
  mck  <- mck %>% group_by(AFFILIATION) %>% filter(level == max(level))
  mc   <- distinct(mck, AFFILIATION, overlap.cluster, .keep_all = TRUE) %>% filter(AFFILIATION != overlap.cluster) %>% group_by(AFFILIATION) %>% mutate(n = n())

  if(any(mc$n > 1)) warning("Some affiliations are in multiple overlap.clusters - this is probably an error in merge.perfect.overlap()")

  l.il <- left_join(l.il, mc, by = c("AFFILIATION"), relationship = "many-to-many")
  s    <- is.na(l.il$overlap.cluster)
  l.il$overlap.cluster[s] <- l.il$AFFILIATION[s]

  by.name    <- l.il %>% group_by(NAME) %>% summarise(circle.esteem = sum(x) - n(), n.linking.positions = sum(x != 0), n.positions = n())
  #by.affil   <- left_join(by.affil, mc %>% select(-level, -n), by = "AFFILIATION")
  by.cluster.name <- l.il %>% group_by(NAME, overlap.cluster) %>% summarise(x = max(x)) %>% group_by(NAME) %>% summarise(circle.esteem.by.cluster = sum(x) - n(), n.linking.positions = sum(x != 0), n.positions = n())
  m               <- level.membership(k, "two-mode", levels = seq) %>% filter(type == 1)
  by.cluster.name <- left_join(by.cluster.name, by.name %>% select(NAME, circle.esteem), by = "NAME") |> left_join(m |> select(-type), by = c("NAME" = "Name"))

  o$name <- by.cluster.name
  o
}



#' Create graphs and adjacency matrixes weighted by circle esteem
#'
#' @param k a list of incidence matrices as produced by \link{k.circles}
#'
#' @return list of data.frames
#' @export esteem.graphs
#'
#' @examples
esteem.graphs <- function(k){
  seq <- seq_along(k)
  mi <- level.membership(k, "two-mode", levels = seq) %>% filter(type == 1)
  mo <- level.membership(k, "two-mode", levels = seq) %>% filter(type == 2)
  # mi$Level[mi$Level <= 1] <- -1
  # mo$Level[mo$Level <= 1] <- -1

  mi$Level[mi$Level <= 1] <- 0
  mo$Level[mo$Level <= 1] <- 0


  bi_adj <- k[[1]] # the bi_adjacency matrix of individuals and affiliations

  # iio
  ioo <- t(bi_adj) * (mo$Level) # multiplying a matrix by a diagonal matrix of scalar (here, the affiliation K-score), scales each column by the corresponding element in diagonal.
  ioo <- t(ioo)
  iio <- cross_sum(ioo)
  iio <- iio - ((bi_adj) %*% t(bi_adj))

  # iii
  ioi   <- bi_adj * (mi$Level) # multiplying a matrix by a diagonal matrix of scalar (here, the affiliation K-score), scales each column by the corresponding element in diagonal.
  iii   <- cross_sum(ioi)

  iii   <- iii - (bi_adj %*% t(bi_adj)) #aktivitets straf
  own   <- bi_adj %*% t(bi_adj)
  own   <- own * (mi$Level)
  diag(own) <- 0
  iii   <- iii - own


  # ooo
  ioo   <- t(t(bi_adj) * (mo$Level)) # multiplying a matrix by a diagonal matrix of scalar (here, the affiliation K-score), scales each column by the corresponding element in diagonal.
  ooo   <- cross_sum(t(ioo))
  a_straf <- (t(bi_adj) %*% bi_adj)
  ooo   <- ooo - a_straf #aktivitets straf
  own   <- t(bi_adj) %*% bi_adj
  own   <- own * (mo$Level -1) #NB! skal den her faktisk trækkes fra!
  diag(own) <- 0
  ooo   <- ooo - own

  #ooi
  ioi <- bi_adj * (mi$Level) # multiplying a matrix by a diagonal matrix of scalar (here, the affiliation K-score), scales each column by the corresponding element in diagonal.
  ooi <- cross_sum(t(ioi))
  ooi <- ooi - (t(bi_adj) %*% (bi_adj))


  out <- list()
  out$iio <- iio
  out$iii <- iii
  out$ooo <- ooo
  out$ooi <- ooi
  out$iio.g <- graph_from_adjacency_matrix(iio, weighted = T, diag = F, mode = "directed")
  out$iii.g <- graph_from_adjacency_matrix(iii, weighted = T, diag = F, mode = "directed")
  out$ooo.g <- graph_from_adjacency_matrix(ooo, weighted = T, diag = F, mode = "directed")
  out$ooi.g <- graph_from_adjacency_matrix(ooi, weighted = T, diag = F, mode = "directed")
  out
}

#' Calculate the k-circle esteem scores
#'
#' @param k a list of incidence matrices as produced by \link{k.circles}
#'
#' @return a list of data.frames
#' @export circle.esteem
#'
#' @examples
circle.esteem <- function(k, esteem = esteem.graphs(k)){

  seq <- seq_along(k)
  mi <- level.membership(k, "two-mode", levels = seq) %>% filter(type == 1)
  mo <- level.membership(k, "two-mode", levels = seq) %>% filter(type == 2)

  g             <- esteem$ooo.g
  est           <- colSums(k[[1]])  %>% enframe(name = "AFFILIATION", value = "Members")
  est$Level     <- mo$Level
  est$esteem    <- diag(esteem$ooi)
  est$out.esteem       <- graph.strength(g, mode = "out")
  est$in.esteem        <- graph.strength(g, mode = "in")

  o       <- list()
  o$affil <- est


  g             <- esteem$iii.g
  est           <- rowSums(k[[1]])  %>% enframe(name = "NAME", value = "Memberships")
  est$Level     <- mi$Level
  est$esteem    <- diag(esteem$iio)
  est$out.esteem       <- graph.strength(g, mode = "out")
  est$in.esteem        <- graph.strength(g, mode = "in")

  o$ind        <- est

  o
}







