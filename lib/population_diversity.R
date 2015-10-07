
alpha_diversity <- function(df, group.col, freq.col) {
  # Returns a variety of diversity indices, including the Gini-Simpson index,
  # the inverse Simpson index, and the Shannon-Weaver index. The SW index is 
  # calculated using the natural logarithm.
  #
  # Arguments:
  # df: a data frame where the rows are species, with a column containing the 
  #   grouping variable and another column containing the proportional 
  #   abundances of each species
  # group: the name of the column defining the grouping variable
  # freq.col: the name of the column containing the proportional abundance of
  #   each species

  call <- substitute(
    df %>% group_by(group) %>%
      summarize(gini.simpson = 1-sum(freq^2),
                inv.simpson = 1/sum(freq^2),
                shannon.weaver = -sum(freq*log(freq))),
    list(group=as.name(group.col), freq=as.name(freq.col)))
  eval(call)
}

chao2 <- function(df, group.col, count.col) {
  # Returns the Chao2 diversity estimate: S* = S_0 + f_1^2/(2*f_2), where S* is 
  # the estimated number of species in a group, S_0 is the total number of 
  # species in a group, f_1 is the number of species occuring in only "plot" or 
  # sample within the group, and f_2 is the number of species occuring in two 
  # "plots".
  #
  # Arguments:
  # df: a data frame where each row contains a species
  # group.col: the column in df that defines the groups
  # count.col: the column in df that contains incidence counts (not abundance!)
  call <- substitute(
    df %>% group_by(group) %>% 
      summarize(chao2 = length(count) + (sum(count==1)^2/(2*sum(count==2)))),
    list(group=as.name(group.col), count=as.name(count.col)))
  eval(call)
}

rarefy <- function(df, size, prob.col, spec.col, group.col) {
  # Returns a rarefied subsample of the species in df such that each group in the
  # specified group.col is rarefied down to a specified size. 
  #
  # Size should probably not be larger than the any of the within-group sizes, but
  # I'm not sure it's wrong either.
  #
  # Arguments:
  # df: a data frame where each row is a species
  # size: the desired within-group size to rarefy down to
  # prob.col: the column in df that gives the within-group proportional abundance
  #   of the species in that row
  # spec.col: the column in df that gives the name of the species in that row
  # group.col: the column in df that gives the group name
  .rarefy <- function(tmp, size, prob.col, spec.col) {
    r <- sample(tmp[[spec.col]], size=size, replace=TRUE, prob=tmp[[prob.col]])
    tmp.counts <- data.frame(table(r))
    colnames(tmp.counts) <- c(spec.col, "count")
    return(tmp.counts)
  }
  # Calls .rarefy on each group in the group.col
  rarefied.counts = ddply(.data=df, .variables=group.col, .fun=.rarefy, size, 
                          prob.col, spec.col)
  return(rarefied.counts)
}