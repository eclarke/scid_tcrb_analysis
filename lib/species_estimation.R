#' Returns a variety of estimates for the true species richness.
#' @param counts a vector of occurrences, or a matrix of counts where the rows 
#' are unique species and the columns are samples 
species_richness <- function(counts, 
                             n.samples = NULL,
                             estimators = c("chao", "jackknife")) {
  # Convert a matrix to a vector of occurrences
  if (is.matrix(counts)) {
    n.samples <- ncol(counts)
    counts <- rowSums(counts > 0) 
  } else if (is.null(n.samples)) {
    stop("Cannot determine the number of samples when given a vector of occurrences.")
  }
  
  # Commonly-used variables
  s0 <- length(counts)
  q1 <- sum(counts == 1)
  q2 <- sum(counts == 2)
  m <- n.samples
  
  if ("chao" %in% estimators) {
    if (q2 > 0) {
      chao <- s0 + (q1^2)/(2*q2)*((m-1)/m)
    } else {
      # bias-corrected Chao
      chao <- s0 + (q1*(q1-1))/(2*(q2+1))*((m-1)/m)
    }
    chao.var <- q2*(((q1/q2)/4)^4 + (q1/q2)^3 + ((q1/q2)/2)^2)
  }
}

jackknife.2 <- function(x, m) {
  n <- length(x)
  # 2nd order jackknife (Smith & Van Belle, 1984)
  jk <- function(x) {
    n <- length(x)
    q1 = sum(x==1)
    q2 = sum(x==2)
    n + (q1*(2*m-3))/m - (q2*(m-2)^2)/(m*(m-1))
  } 
  
  species <- jk(x)
  
  jk.variance <- ((n-1) / n) * sum(sapply(c(1:n), function(i) {
    (jk(x[-i]) - species)^2
  }))
  
  return(list(species=species, jk.var=jk.variance, jk.std=sqrt(jk.variance)))
}