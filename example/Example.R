library(tidyverse)
library(Matrix)
library(igraph)
library(eliter)
library(kcircles)


data(den)
den    <- den[den$SOURCE != "Events",]
kill   <- c("Arbejdsmiljørådet (Rådsmedlemmer)")
den    <- den %>% filter(!AFFILIATION %in% kill)
den    <- as.den(den)

incidence <- xtabs(~NAME + AFFILIATION, den, sparse = TRUE)
incidence@x[]        <- 1
circle.solution   <- k.circles(incidence, 2, check.for.nested = TRUE)
k                 <- circle.solution

lev.solution      <- level.membership(k, mode = "two-mode")
table(lev.solution$Level, lev.solution$type)


esteem            <- circle.esteem(k)
cluster.esteem    <- circle.esteem.for.clusters(k)

# Calculate correlates ------
library(birankr)

centrality.measures <- function(k){
incidence <- k[[1]]
esteem            <- circle.esteem(k)

g.w.i <- (incidence %*% t(incidence)) |> graph_from_adjacency_matrix(diag = FALSE, mode = "undirected", weighted = TRUE)
g.w.a <- (t(incidence) %*% incidence) |> graph_from_adjacency_matrix(diag = FALSE, mode = "undirected", weighted = TRUE)

d             <- tibble(AFFILIATION = V(g.w.a)$name,  strength = strength(g.w.a))
d$betweenness <- betweenness(g.w.a, weights = 1/E(g.w.a)$weight, cutoff = 5)
d$closeness   <- closeness(g.w.a, weights = 1/E(g.w.a)$weight, cutoff = 5)
d$eigen       <- eigen_centrality(g.w.a, weights = E(g.w.a)$weight)$vector
d$page_rank   <- page_rank(g.w.a, weights = E(g.w.a)$weight)$vector
d$coreness    <- coreness(g.w.a, mode = "all")

d$birank.HITS        <- birankr::bipartite_rank(incidence, return_mode = "columns", normalizer = "HITS")
d$birank.BiRank      <- birankr::bipartite_rank(incidence, return_mode = "columns", normalizer = "BiRank")
d$birank.CoHITS      <- birankr::bipartite_rank(incidence, return_mode = "columns", normalizer = "CoHITS")

est.a        <- esteem$affil
d            <- left_join(d, est.a)
d
}


centrality.measures.i <- function(k){

  incidence <- k[[1]]
  esteem            <- circle.esteem(k)

  g.w.i <- (incidence %*% t(incidence)) |> graph_from_adjacency_matrix(diag = FALSE, mode = "undirected", weighted = TRUE)
  g.w.a <- (t(incidence) %*% incidence) |> graph_from_adjacency_matrix(diag = FALSE, mode = "undirected", weighted = TRUE)

  d             <- tibble(NAME = V(g.w.i)$name,  strength = strength(g.w.i))
  d$betweenness <- betweenness(g.w.i, weights = 1/E(g.w.i)$weight, cutoff = 5)
  d$closeness   <- closeness(g.w.i, weights = 1/E(g.w.i)$weight, cutoff = 5)
  d$eigen       <- eigen_centrality(g.w.i, weights = E(g.w.i)$weight)$vector
  d$page_rank   <- page_rank(g.w.i, weights = E(g.w.i)$weight)$vector
  d$coreness    <- coreness(g.w.i, mode = "all")

  d$birank.HITS        <- birankr::bipartite_rank(incidence, return_mode = "rows", normalizer = "HITS")
  d$birank.BiRank      <- birankr::bipartite_rank(incidence, return_mode = "rows", normalizer = "BiRank")
  d$birank.CoHITS      <- birankr::bipartite_rank(incidence, return_mode = "rows", normalizer = "CoHITS")

  est.a        <- esteem$ind
  d            <- left_join(d, est.a)
  d
}

centrality.measures(k) |> select(-AFFILIATION) |> cor(use = "pairwise.complete", method = "spearman") |> round(2)

cm <- centrality.measures(k)
cm.den <- cm

# Longlinks ----
load("~/Dropbox/LongLinks/R/LongLinks_Affil/saved/2_analysis_All.Rda")
#k.set <- l.k[c("1920", "1940", "1960", "1980", "2000", "2013")]
k.set <- l.k[c("2013")]
k <- k.set$"2013"

last(k) |> dim()

l.cm  <- map(k.set, centrality.measures)
# l.cm$"1920" |> arrange(-esteem) |> select(-AFFILIATION) |> cor(use = "pairwise.complete", method = "spearman") |> round(2)
# l.cm$"1940" |> arrange(-esteem)
# l.cm$"1960" |> arrange(-esteem)
# l.cm$"1980" |> arrange(-esteem)
# l.cm$"2000" |> arrange(-esteem)
l.cm$"2013" |> arrange(-esteem) |> View()
l.cm$"2013" |> arrange(-esteem) |> select(-AFFILIATION) |> cor(use = "pairwise.complete", method = "spearman") |> round(2)

cm.affil <- l.cm$"2013"

f <- function(x) rank(-x) <= 200
cmx <- cm |> select(-AFFILIATION) |> transmute(across(everything(), .fns = f))
cmx <- cmx |> as.matrix()
t(cmx) %*% cmx

# Individ
l.cm  <- map(k.set, centrality.measures.i)
l.cm$"2013" |> arrange(-esteem) |> select(-NAME) |> cor(use = "pairwise.complete", method = "spearman") |> round(2)
cm <- l.cm$"2013"
cm.i <- cm

f <- function(x) rank(-x) <= 200
cmx <- cm |> select(-NAME) |> transmute(across(everything(), .fns = f))
cmx <- cmx |> as.matrix()


# Plot ----
library(corrplot)
corrplot::corrplot()

cor.den   <- cm.den |> select(esteem, Level, Members, strength, coreness, eigen, closeness, betweenness, page_rank, birank.HITS) |> cor(use = "pairwise.complete", method = "kendall")
cor.affil <- cm.affil |> select(esteem, Level, Members, strength, coreness, eigen, closeness, betweenness, page_rank, birank.HITS) |> cor(use = "pairwise.complete", method = "kendall")

corrplot(cor.den, method = "ellipse")
corrplot.mixed(cor.affil, upper = "ellipse")


cm.affil |> mutate(rank.eigen = rank(desc(eigen)), rank.esteem = rank(desc(esteem)), diff = rank(desc(eigen)) - rank(desc(esteem))) |> View()
#####



s <- c("Københavns Universitet", "Copenhagen Business School", "Aarhus Universitet", "Danmarks Tekniske Universitet", "Syddansk Universitet")
ind <- k.set$"2013"[[1]]
incidence  <- ind[, !colnames(ind) %in% s]
incidence@x[]        <- 1
circle.solution   <- k.circles(incidence, 4, check.for.nested = FALSE)
k                 <- circle.solution

lev.solution      <- level.membership(k, mode = "two-mode")
esteem            <- circle.esteem(k)
cm                <- centrality.measures(k) |> select(-AFFILIATION) |> cor(use = "pairwise.complete", method = "spearman") |> round(2)
cm




cm <- centrality.measures.i(k)
l.cm$"2013" |> arrange(-esteem) |> select(-NAME) |> cor(use = "pairwise.complete", method = "spearman") |> round(2)
cm <- l.cm$"2013"
cm.i <- cm


# d             <- tibble(NAME = V(g.w.i)$name,  strength = strength(g.w.i))
# d$betweenness <- betweenness(g.w.i, weights = 1/E(g.w.i)$weight, cutoff = 5)
# d$closeness   <- closeness(g.w.i, weights = 1/E(g.w.i)$weight, cutoff = 5)
# d$coreness    <- coreness(g.w.i, mode = "all")
# d$birank.HITS      <- birankr::bipartite_rank(incidence, return_mode = "rows", normalizer = "HITS")
# d$birank.BiRank      <- birankr::bipartite_rank(incidence, return_mode = "rows", normalizer = "BiRank")
#
# est.i        <- esteem$ind
# d            <- left_join(d, est.i)
#
# d |> select(-NAME) |> cor(use = "pairwise.complete", method = "spearman")
#


View(d)
