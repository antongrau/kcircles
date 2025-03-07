library(eliter)
library(tidyverse)
library(Matrix)
library(knitr)
library(RColorBrewer)
library(ggthemes)
library(scales)
library(ggrepel)
library(igraph)


x <- matrix(c(1,2,0,1,2, 1), nrow = 2, byrow = T)

# cross_sum <- function(x) {
#   M1 <- x %*% t((x != 0)) #+ t(x %*% t(x != 0))
#   M  <- (M1 + t(M1)) / 2  #- Matrix::Diagonal(diag(M1), n = length(diag(M1)))
#   M
# }
# M[b,a]

data(den)
den    <- den[den$SOURCE != "Events",]
kill   <- c("Arbejdsmiljørådet (Rådsmedlemmer)")
den    <- den %>% filter(!AFFILIATION %in% kill)
den    <- as.den(den)

incidence <- xtabs(~NAME + AFFILIATION, den, sparse = TRUE)
incidence@x[]        <- 1
circle.solution   <- k.circles(incidence, 3, check.for.nested = TRUE)
lev.solution      <- level.membership(circle.solution, mode = "two-mode")
l <- lev.solution %>% filter(type == "1")
table(l$Level)

k <- circle.solution



sum.of.k.scores <- function(k){
  seq <- seq_along(k)
  mi <- level.membership(k, "two-mode", levels = seq) %>% filter(type == 1)
  mo <- level.membership(k, "two-mode", levels = seq) %>% filter(type == 2)
  bi_adj <- k[[1]] # the bi_adjacency matrix of individuals and affiliations

  # iio
  ioo <- t(bi_adj) * (mo$Level - 1) # multiplying a matrix by a diagonal matrix of scalar (here, the affiliation K-score), scales each column by the corresponding element in diagonal.
  ioo <- t(ioo)
  iio <- cross_sum(ioo)
  #iio <- iio - ((bi_adj) %*% t(bi_adj))
  iio2 <- iio * diag(iio)

  # iii
  ioi <- bi_adj * (mi$Level - 1) # multiplying a matrix by a diagonal matrix of scalar (here, the affiliation K-score), scales each column by the corresponding element in diagonal.
  iii <- cross_sum(ioi)

  iii <- iii - ((bi_adj) %*% t(bi_adj))
  adj <- ((bi_adj %*% t(bi_adj)))
  diag(adj) <- 0
  table(adj@x)
  adj <- adj * (mi$Level - 1)
  table(adj@x)

  iii <- iii - adj
  diag(iii) <- 0
  table(iii@x)
  a <- "Thorkild Engell Jensen"
  b <- "Per Påskesen"

  adj[a, b]
  adj[b, a]

  iii[a, b]
  iii[b, a]

  # ooo
  ioo <- t(bi_adj) * (mo$Level - 1) # multiplying a matrix by a diagonal matrix of scalar (here, the affiliation K-score), scales each column by the corresponding element in diagonal.
  ioo <- t(ioo)
  ooo <- cross_sum(t(ioo))
  ooo <- ooo - (t(bi_adj) %*% (bi_adj))
  adj <- ((t(bi_adj) %*% (bi_adj))) * (mo$Level - 1)
  ooo <- ooo - adj

  #ooi
  ioi <- bi_adj * (mi$Level - 1) # multiplying a matrix by a diagonal matrix of scalar (here, the affiliation K-score), scales each column by the corresponding element in diagonal.
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

k <- circle.solution
res  <- sum.of.k.scores(k)

seq <- seq_along(k)
mi <- level.membership(k, "two-mode", levels = seq) %>% filter(type == 1)
mo <- level.membership(k, "two-mode", levels = seq) %>% filter(type == 2)



g <- res$ooo.g
est <- graph.strength(g, mode = "all") %>% enframe()
est$out <- graph.strength(g, mode = "out")
est$"in"  <- graph.strength(g, mode = "in")


est$memb_n <- colSums(k[[1]])
est$level <- mo$Level - 1
est$value_2 <- est$value - (est$level*est$memb_n)
est$dif <- est$value - est$value_2
est$esteem.o <- diag(res$ooi)
est$ml <- est$level / est$memb_n
est


g <- res$iii.g
est <- graph.strength(g, mode = "all") %>% enframe()
est$out <- graph.strength(g, mode = "out")
est$"in"  <- graph.strength(g, mode = "in")


est$memb_n <- rowSums(k[[1]])
est$level <- mi$Level - 1
est$value_2 <- est$value - (est$level*est$memb_n)
est$dif <- est$value - est$value_2
est$esteem.o <- diag(res$iio)
est$ml <- est$level / est$memb_n


# ego.plot

V(g)$level <- left_join(tibble(Name = V(g)$name), mi) |> pull(Level)

a <- "Bjarne Corydon"
a <- "Helle Thorning-Schmidt"
a <- "Ingerlise Buck"

g

eg <- g |> make_ego_graph(nodes = a)
eg <- eg[[1]]
ew <- E(eg)$weight
ew

p <- graph.plot(eg, layout = layout_with_fr(eg, weights = NA), edge.color = ew, edge.alpha = 1, edge.size = 0.2, midpoints = TRUE, vertex.size = V(eg)$level)
p + scale_color_gradient2(mid = "grey90", low = "darkred")

# Longlinks -----
load("~/Dropbox/LongLinks/R/LongLinks_Affil/saved/2_analysis_All.Rda")


circle.solution <- l.k$"2013"
k <- circle.solution
res  <- sum.of.k.scores(circle.solution)

seq <- seq_along(k)
mi <- level.membership(k, "two-mode", levels = seq) %>% filter(type == 1)
mo <- level.membership(k, "two-mode", levels = seq) %>% filter(type == 2)



g <- res$ooo.g
est <- graph.strength(g, mode = "all") %>% enframe()
est$out <- graph.strength(g, mode = "out")
est$"in"  <- graph.strength(g, mode = "in")


est$memb_n <- colSums(k[[1]])
est$level <- mo$Level - 1
est$value_2 <- est$value - (est$level*est$memb_n)
est$dif <- est$value - est$value_2
est$esteem.o <- diag(res$ooi)
est$ml <- est$level / est$memb_n


g <- res$iii.g
est <- graph.strength(g, mode = "all") %>% enframe()
est$out <- graph.strength(g, mode = "out")
est$"in"  <- graph.strength(g, mode = "in")


est$memb_n <- rowSums(k[[1]])
est$level <- mi$Level - 1
est$value_2 <- est$value - (est$level*est$memb_n)
est$dif <- est$value - est$value_2
est$esteem.o <- diag(res$iio)
est$ml <- est$level / est$memb_n
est






# JUNK _--------




seq <- seq_along(k)
mi <- level.membership(k, "two-mode", levels = seq) %>% filter(type == 1)
mo <- level.membership(k, "two-mode", levels = seq) %>% filter(type == 2)
i <- k[[1]] # the bi_adjacency matrix of individuals and affiliations

il <- i * (mi$Level - 1) # multiplying a matrix by a diagonal matrix of scalar (here, the affiliation K-score), scales each column by the corresponding element in diagonal.
il2 <- cross_sum(il) # here we compute the cross_sum A * A^T.
#il2 <- il2 - ((i) %*% t(i)) # and substracts the cross_product of the original matrix
adj <- ((i %*% t(i)) != 0) * (mi$Level - 1)
il2 <- il2 - adj
summary(il2@x)
max(adj)
il2["Peder Steen Andersen Philipp", ] |> enframe() |> View()
esteem <- il2 %>% graph_from_adjacency_matrix(., weighted = T, diag = F, mode = "directed")
est <- graph.strength(esteem, mode = "all") %>% enframe()
est$out <- graph.strength(esteem, mode = "out")
est$"in"  <- graph.strength(esteem, mode = "in")


est$memb_n <- rowSums(i)
est$level <- mi$Level - 1
est$value_2 <- est$value - (est$level*est$memb_n)
est$dif <- est$value - est$value_2
est$memberships <- membership.vector(est$name, den)

diag(il2)

plot(esteem$weight %>% sort())

s <- il["Bjarne Corydon", ] > 0

ss <- il[rowSums(il[,s]) > 0, s]

as.matrix(ss[c("Bjarne Corydon", "Sten Kristensen"),]) %>% graph_from_adjacency_matrix(., weighted = T, diag = F, mode = "undirected")  %>% as.data.frame() %>% View()
cross_sum(ss) %>% graph_from_adjacency_matrix(., weighted = T, diag = F, mode = "undirected")  %>% igraph::as_data_frame()  %>% View()

A <- matrix(1:12, nrow = 3, ncol = 4)  # Example matrix (3x4)
v <- c(2, 3, 4, 5)                     # Vector of length N (same as number of columns)

A_scaled <- A * matrix(v, nrow = nrow(A), ncol = length(v), byrow = TRUE)


sum.of.k.scores.members <- function(k, drop.last = TRUE){
  seq <- seq_along(k)
  if(identical(drop.last, TRUE)) seq <- seq[-length(seq)]

  m <- level.membership(k, "two-mode", levels = seq) %>% filter(type == 2)
  i <- k[[1]]

  stopifnot(all.equal(colnames(i), m$Name))
  il <- t(i) * (m$Level - 1)
  #il <- t(il)
  il2 <- cross_sum(il)
  g <- graph_from_adjacency_matrix(il2, weighted = T, diag = F, mode = "undirected")
  e <- g %>% igraph::as_data_frame()


  rs <- rowSums(il)
  rs
}

sum.of.k.scores.members <- function(k, drop.last = TRUE){
  seq <- seq_along(k)
  if(identical(drop.last, TRUE)) seq <- seq[-length(seq)]

  m <- level.membership(k, "two-mode", levels = seq) %>% filter(type == 2)
  i <- k[[1]]

  stopifnot(all.equal(colnames(i), m$Name))
  il <- t(i) * (m$Level - 1)
  il <- t(il)
  rs <- rowSums(il)
  rs
}


k <- l.k$"2020"

sum.of.k.scores.members <- function(k, drop.last = FALSE){
  seq <- seq_along(k)
  if(identical(drop.last, TRUE)) seq <- seq[-length(seq)]

  m <- level.membership(k, "two-mode", levels = seq) %>% filter(type == 2)
  i <- k[[1]]

  stopifnot(all.equal(colnames(i), m$Name))
  il <- t(i) * (m$Level - 1)
  il <- t(il)
  rs <- rowSums(il)
  rs
}




soks.affil  <- sum.of.k.scores(circle.solution) %>% enframe()
soks.member <- sum.of.k.scores.members(circle.solution) %>% enframe()

table(lev.solution$Level)

# Fra tags til sektor inddeling ----

sector.tags <- list()
sector.tags$"Erhvervsliv" <- c("Corporation")
sector.tags$"Erhvervsorganisationer" <- c("Business association", "Employers association")
sector.tags$"Fagbevægelse" <- c( "Unions", "Standsforening", "A-kasse", "Union controlled")
sector.tags$"Politik"      <- c("Politics", "Parliament", "Political party")
sector.tags$"Stat"         <- c("State administration", "Ministry", "State corporation", "Military",
                                "Public leaders", "Commission", "Politics", "Parliament", "State business")
sector.tags$"Videnskab og uddannelse" <- c("Science", "Education", "Universities")

# Farver til sektorer -----

sektor.pal <- c("Erhvervsliv" = "PuBu", "Fagbevægelse" = "Reds", "Videnskab og uddannelse" = "Greens",
                "Erhvervsorganisationer" = "BuPu" , "Stat" = "RdPu", "Politik" = "Greys", "Organisation" = "Oranges",
                "Kultur" = "YlOrBr")
l.colors  <- lapply(sektor.pal, brewer.pal, n = 9)

intensity <- 5

fill.scale <- c("Erhvervsliv" = l.colors$Erhvervsliv[intensity],
                "Fagbevægelse" = l.colors$Fagbevægelse[intensity + 2],
                "Videnskab og uddannelse" = l.colors$`Videnskab og uddannelse`[intensity - 1],
                "Erhvervsorganisationer" = l.colors$Erhvervsorganisationer[intensity - 2],
                "Stat" = l.colors$Stat[intensity],
                "Politik" = l.colors$Politik[intensity + 2],
                "Organisation" = l.colors$Organisation[intensity - 1],
                "Kultur" = l.colors$Kultur[intensity - 1])

disc.scale <- brewer.pal(9, "YlGnBu")
dic.scale <- c(disc.scale[1], disc.scale[length(disc.scale)])
disc.scale.3 <- disc.scale[c(1, 5, 9)]


# Data til plot -----
affil.sector       <- tags.to.sectors(den, sector.tags = sector.tags, sector.membership = TRUE, other = "Andet", mutually.exclusive = F, silent = F)
affil.sector       <- affil.sector %>% bind_rows(.id = "Sektor") %>% as_tibble() %>% select(AFFILIATION, Sektor)  %>% group_by(AFFILIATION) %>% summarise(Sektor = first(Sektor))
d                  <- affil.sector %>% select(Name = AFFILIATION, Sektor)
soks               <- bind_rows(soks.affil, soks.member) %>% rename("K-sum" = value, Name = name)

ind <- circle.solution[[9]]

g      <- ind %>% graph_from_incidence_matrix()
d      <- tibble(Name = V(g)$name) %>% left_join(d, by = "Name") %>% left_join(., soks, by = "Name")
d$type <- V(g)$type

fill <- d$Sektor
size <- d$`K-sum`
size[is.na(size)] <- min(size, na.rm = TRUE)/2

ew   <- edge.betweenness.estimate(g, cutoff = 8)
ef     <- tibble(Name = as_edgelist(g)[,2]) %>% left_join(d, by = "Name")
ec     <- ef$Sektor

# Layouts

lay    <- layout_with_fr(g, weights = 1/log(ew))
lay.labels <- tibble(Name = d$Name, X=lay[,1], Y = lay[, 2])
lay.affil <- lay.labels %>% filter(d$type == TRUE)

# Plot af two-mode netværk
p      <- graph.plot(g, lay, vertex.shape = V(g)$type, vertex.fill = fill, edge.color = ec,  vertex.size = size, edge.size = 0.2, edge.alpha = ew, edge.order = ew, norm.coords = FALSE, vertex.order = size)
p      <- p + scale_size_continuous(range = c(1, 6), guide = "none") + scale_alpha(range = c(0.3, 1), guide = "none")
p      <- p + scale_fill_manual(values = fill.scale,  aesthetics = c("color", "fill"), drop = FALSE)
p      <- p + scale_shape_manual(values = c(4, 21), guide = "none")
p      <- p + guides(fill = guide_legend(override.aes = list(shape = 21, size = 3)))
p      <- p + coord_equal()
p

p + geom_text(data = lay.labels, aes(x = X, y = Y, label = Name), check_overlap = T)
p + geom_text(data = lay.affil, aes(x = X, y = Y, label = Name), check_overlap = T)
p + geom_label_repel(data = lay.affil, aes(x = X, y = Y, label = Name))


# Centreret analyse -----

den.corp <- den[den$SOURCE == "Corporations",]
den.chair <- den.corp[den.corp$ROLE == "Chairman",]
chairs <- den.chair$NAME %>% unique()
chairs<- chairs[chairs %in% rownames(incidence)]
ind <- incidence[chairs, ]

ind <- xtabs(~NAME+AFFILIATION, den.corp, sparse = TRUE)
ind@x[]        <- 1



circle.solution   <- k.circles(ind, 2, check.for.nested = FALSE)

lev.solution      <- level.membership(circle.solution, mode = "two-mode")
table(lev.solution$Level)
View(lev.solution)
