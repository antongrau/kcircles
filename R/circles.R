#' K-circles decomposition
#'
#' @param incidence a sparse incidence matrix
#' @param minimum.memberships the minimum number of memberships for the individuals (rows)
#'
#' @return an object of class "k.circles"
#' @export k.circles
#'
#' @examples
#' data(den)
#' den         <- den[den$SOURCE != "Events",]
#' incidence   <- xtabs(~NAME + AFFILIATION, droplevels(den), sparse = TRUE)
#' l.inc       <- k.circles(incidence, 3, check.for.nested = TRUE)
#' level.membership(l.inc)
#' l.inc[[5]] %>% colSums() %>% sort() %>% as.matrix()
#' l.inc[[5]] %>% rowSums() %>% sort() %>% as.matrix()

k.circles <- function(incidence, minimum.memberships = 3, check.for.nested = TRUE){


  ##############
  # Tests


  # Is it a sparse matrix?
  if(!inherits(incidence, "dgCMatrix")) stop("incidence has to be a sparse matrix of class dgCMatrix. With xtabs you can set sparse to TRUE and get a valid matrix.")

  # Check for multiple memberships and odd values
  if(any(incidence@x != 1)) warning("incidence has values other than 1 and . (the sparse version of 0). Try table(incidence@x) to see them.")

  # Backup of the incidence matrix
  original.colnames           <- colnames(incidence)
  colnames(incidence)         <- 1:ncol(incidence)

  inc         <- incidence



  # k is the minimum number of members
  k  <- 0

  # l.inc is a list of incidence matrices
  l.inc       <- list()

  # j is the minimum number of memberships for any individual
  level.up <- function(inc, k, j = 3, check.for.nested = TRUE){

    test.mat <- function(inc, j, k){

      cond                       <- !is.null(dim(inc))
      if(cond){
        cond              <- any(
          c(
            any(Matrix::rowSums(inc) < j), # Is there any individuals with less than j positions
            any(Matrix::colSums(inc) < k)  # Is there any affiliations with less than k members
          ))
      }
      cond
    }

    # Levelling up
    while(test.mat(inc, j, k)){
      # Removing rows
      inc.t       <- inc[Matrix::rowSums(inc) >= j, ]
      if(is.null(dim(inc.t))) break
      inc         <- inc.t  # Keep only those members with j or more positions

      # Removing columns
      inc.t       <- inc[, Matrix::colSums(inc) >= k]  # Keep only those affiliations with more than k members
      if(is.null(dim(inc.t))) break
      inc         <- inc.t

      # Merging overlapping affiliations
      if(identical(check.for.nested, TRUE)){
        inc         <- eliter:::merge.perfect.overlap(inc, combine.labels = "&")
      }

    }
    inc
  }

  while({
    cs <-  Matrix::colSums(inc)
    if(length(cs) == 0) cs <- 0
    k <= min(cs) & ncol(inc) > minimum.memberships # While k is smaller than the lowest number of members and the number of affiliations is larger than the minimum number of memberships
  }){
    k           <- k + 1
    tmp         <- level.up(inc, k, j = minimum.memberships, check.for.nested = check.for.nested)
    inc         <- tmp
    l.inc[[k]]  <- inc

  }

  # Clean up and class
  l.inc        <- c(incidence, l.inc)
  l.inc        <- compact(l.inc)
  class(l.inc) <- append("k.circle", class(l.inc))

  # Here we reconstruct the original names and we add the merge.clusters as an attribute.

  give.the.old.names.back <- function(i){
    n            <- str_split(colnames(i), pattern = "&") %>% map(str_trim) %>% map(as.numeric)
    colnames(i)  <- map_chr(n, ~paste(original.colnames[.x], collapse = " & "))
    names(n)     <- colnames(i)
    n            <- map(n, ~ enframe(original.colnames[.x])) %>% bind_rows(.id = "name") %>% transmute(AFFILIATION = value, overlap.cluster = name)
    attr(i, "col.circle.merge.cluster") <- n
    i
  }

  l.inc <- map(l.inc, give.the.old.names.back)
  attr(l.inc, "minimum.memberships") <- minimum.memberships
  l.inc

}

merge.perfect.overlap <- function(incidence, combine.labels = "&", cutoff = 1){
  # This functions throws an error if any of the affiliations are empty

  # Goal: Merge perfectly overlapping affiliations
  # Combine their labels and remove one of the columns.
  # They merge into either the largest affiliation or to the first in the order
  # It is run when the incidence has been pruned. So affiliations with just a single member will have disappeared
  # The merged affiliation will have all its values set to 0

  adj                    <- Matrix::crossprod(incidence)
  affil.members          <- Matrix::diag(adj)
  names(affil.members)   <- rownames(adj)
  adj.s                  <- adj / affil.members
  diag(adj.s)            <- 0
  merge.ind              <- Matrix::which(adj.s >= cutoff, arr.ind = TRUE)
  merge.ind              <- tibble(row = merge.ind[, 1], col = merge.ind[, 2]) # Row and column indices

  s                      <- merge.ind %>% apply(1, sort)  %>% t()
  if(nrow(s) > 1) merge.ind <- merge.ind %>% filter(!duplicated(s)) # Check if two of equal size are there.

  merge.ind$row.name <- colnames(adj)[merge.ind$row]
  merge.ind$col.name <- colnames(adj)[merge.ind$col]



  s                      <- merge.ind$col %in% merge.ind$row  # col må ikke være i row - fordi vi må ikke slette noget der er blevet merget ind i.
  merge.ind              <- merge.ind %>% filter(!s)



  # If one of the nodes are there twice - it means that it will be merged two times. That is not a great problem as it does not influence the scores - except for merged clusters.
  # graph_from_edgelist(merge.ind[, 3:4] %>% as.matrix()) %>% plot() # Here we check if any node has an outdegree of 2 - that is not allowed!
  merge.ind             <- merge.ind[!duplicated(merge.ind$row), ]

  if(nrow(merge.ind) == 0) return(incidence)

  if(identical(combine.labels, FALSE) == FALSE){
    for(i in 1:nrow(merge.ind)){
      cr             <- merge.ind$row[i]
      cc             <- merge.ind$col[i]
      label          <- paste(colnames(incidence)[cc], combine.labels,  colnames(incidence)[cr])
      colnames(incidence)[cc] <- label
    }
  }

  add_together <- function(x,i){
    row <- merge.ind[i, ]$row
    col <-  merge.ind[i, ]$col

    change <- x[, row] > 0 & x[, col] == 0
    x[change, col] <- x[change, row]
    x
  }

  for(i in 1:nrow(merge.ind)){
    incidence <- add_together(incidence, i)

  }

  incidence[, merge.ind$row] <- 0
  drop0(incidence)
}

#' Level membership from minimal membership decomposition
#'
#' @param l.inc a list of nested incidence matrices
#'
#' @return a tibble with rownames and level membership
#' @export
#'
#' @examples
level.membership <- function(l.inc, mode = c("ind", "affil", "two-mode"), levels = seq_along(l.inc)){



  # When we merge affilations the naming of the affil and two-mode will be more complicated
  mode           <- match.arg(mode)
  l.inc          <- l.inc[levels]
  l              <- length(l.inc)

  # Membership for individuals
  membership  <- map(l.inc, rownames) %>% imap(~ tibble(Name = .x, Level = .y)) %>%
    bind_rows() %>% arrange(Name)

  mem         <- membership %>% group_by(Name) %>% summarise(Level = max(Level) -1)
  mem.ind     <- tibble(Name = rownames(l.inc[[1]])) %>% left_join(., mem, by = "Name")

  if(mode == "ind") return(mem.ind)

  # Membership for affiliations

  inc        <- l.inc[[1]]
  f          <- function(x, inc) Matrix::colSums(inc[x,]) %>% as_tibble(rownames = "Name")
  level.mem  <- map(l.inc, rownames) %>% map(., f, inc = inc) %>% set_names(0:(l-1)) %>% bind_rows(.id = "level")
  level.mem  <- level.mem %>% mutate(level = as.numeric(level)) %>% filter(value >= level)
  level.mem  <- level.mem %>% group_by(Name) %>% summarise(Level = max(level))

  mem.affil   <- tibble(Name = colnames(inc)) %>% left_join(., level.mem, by = "Name")

  if(mode == "affil") return(mem.ind)

  # Membership for two-mode
  # We assume that Igraph or other will sort first by rows and then columns

  mem.two         <- bind_rows(mem.ind, mem.affil, .id = "type")

  if(mode == "two-mode")  return(mem.two)
}


