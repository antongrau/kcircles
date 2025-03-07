# Data and cleaning ----
library(eliter)
library(kcircles)

library(tidyverse)

# Read data ----

# Theoretical example ----
el            <- list()
el$david      <- tibble(NAME = "David", AFFILIATION = LETTERS[1:3])
el$boris      <- tibble(NAME = "Boris", AFFILIATION = LETTERS[c(1:5, 7)])
el$boris      <- tibble(NAME = "Boris", AFFILIATION = LETTERS[c(1:4)])
el$tony       <- tibble(NAME = "Tony",  AFFILIATION = LETTERS[5:6])
el$theresa    <- tibble(NAME = "Theresa", AFFILIATION = LETTERS[c(1,4, 6)])
el$gordon     <- tibble(NAME = "Gordon", AFFILIATION = LETTERS[c(2, 3, 6, 7)])

el$iso1     <- tibble(NAME = paste("iso", LETTERS[1], 1:5), AFFILIATION = LETTERS[1])
el$iso2     <- tibble(NAME = paste("iso", LETTERS[2], 1:5), AFFILIATION = LETTERS[2])
el$iso3     <- tibble(NAME = paste("iso", LETTERS[3], 1:4), AFFILIATION = LETTERS[3])
el$iso4     <- tibble(NAME = paste("iso", LETTERS[4], 1:4), AFFILIATION = LETTERS[4])
el$iso5     <- tibble(NAME = paste("iso", LETTERS[5], 1:4), AFFILIATION = LETTERS[5])
el$iso6     <- tibble(NAME = paste("iso", LETTERS[6], 1:2), AFFILIATION = LETTERS[6])

el$isoaffil <- tibble(NAME = paste("isoaffil", LETTERS[8], 1:3), AFFILIATION = LETTERS[8])

el$isocom1   <- tibble(NAME = paste("isocom", LETTERS[9], 1:3), AFFILIATION = LETTERS[9])
el$isocom2   <- tibble(NAME = paste("isocom", LETTERS[10], 1:3), AFFILIATION = LETTERS[10])
el$sam       <- tibble(NAME = "Sam", AFFILIATION = LETTERS[c(9, 10)])
el$eric      <- tibble(NAME = "Eric", AFFILIATION = LETTERS[c(9, 10)])

#el$iso7     <- tibble(NAME = paste("iso", LETTERS[3], 1:4), AFFILIATION = LETTERS[7])
el.theo     <- bind_rows(el)

# Plotting -----
incidence   <- xtabs(~NAME + AFFILIATION, el.theo, sparse = TRUE)
g           <- graph_from_biadjacency_matrix(incidence)
l.inc       <- k.circles(incidence, minimum.memberships = 2, check.for.nested = FALSE)
l.inc       <- l.inc %>% unique()

map(l.inc, dim)
lev.mem     <- level.membership(l.inc, mode = "two-mode")
esteem      <- circle.esteem(l.inc)
est         <- bind_rows(esteem$affil |> rename(name = AFFILIATION), esteem$ind |> rename(name = NAME), .id = "type")
est         <- left_join(tibble(name = V(g)$name), est)

edge.color.by.factor <- function(graph, f){
  el         <- as_edgelist(graph)
  el         <- data.frame(X = el[,1], Y = el[,2])
  d          <- data.frame(vertex.attributes(graph))
  d$f        <- f
  elx        <- left_join(el, d, by = c("X" = "name"))
  ely        <- left_join(el, d, by = c("Y" = "name"))

  out        <- cbind(elx$f, ely$f) |> apply(1, min)
  #out[out != ely$f] <- NA
  out
}

ec <- edge.color.by.factor(g, est$Level |> as_factor())
ec <- ec -1
ec <- ec |> as_factor()

p           <- graph.plot(g, text = FALSE, vertex.fill = factor(lev.mem$Level), vertex.size = V(g)$type, vertex.shape = factor(lev.mem$type), vertex.color = factor(lev.mem$Level), edge.color = ec, text.vjust = -1, edge.size = 0.2)
p           <- p + scale_shape_manual(values = c(21, 22))
p           <- p + scale_shape_manual(values = c(120, 21), guide = "none")
p           <- p + scale_size(range = c(2,6))
p           <- p + scale_color_manual(values = c("0" = "salmon", "1" = "cornflowerblue", "2" = "purple" ,  "3" = "#67000d"), na.value = "grey60") + scale_fill_manual(values = c("0" = "salmon", "1" = "cornflowerblue", "2" = "purple" ,  "3" = "#67000d"),  na.value = "grey60")
p + scale_size_manual(values = c("FALSE" = 4, "TRUE" = 4), guide = "none") + coord_fixed()




p           <- graph.plot(g, text = TRUE, vertex.fill = est$esteem, vertex.size = V(g)$type, vertex.shape = factor(lev.mem$type), vertex.color = est$esteem, text.vjust = -1, edge.size = 0.2)
p           <- p + scale_shape_manual(values = c(21, 22))
p           <- p + scale_shape_manual(values = c(120, 21), guide = "none")
p           <- p + scale_size(range = c(2,6))
p           <- p + scale_fill_viridis_c(option = "A", na.value = "grey60", direction = -1) + scale_color_viridis_c(option = "A", na.value = "grey60", direction = -1)
p + scale_size_manual(values = c("FALSE" = 4, "TRUE" = 4), guide = "none") + coord_fixed()





