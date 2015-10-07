#' Reorders a factor's levels to match the order given by calling unique() on 
#' the factor itself. Useful for character vectors in data frames that are
#' already ordered correctly, but when coerced to a factor, the levels are in
#' alphabetical order.
reorder_levels_by_uniq <- function(f) {
  if (!(is.factor(f))) f <- as.factor(f)
  f <- factor(f, levels=levels(f)[match(unique(f), levels(f))])
  return(f)
}

