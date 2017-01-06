

# test

require(Rcpp)
require(inline)
require(RcppArmadillo)

for (r in 1:repCount) {
    for (k in 2:maxK) {
      ## mCount is possible number of times that two sample occur in same random sample,
      ## independent of k mCount stores number of times a sample pair was sampled
      ## together.
      mCount <- fun(rep(1, length(parRes[[r]][[1]][[3]])), mCount, parRes[[r]][[1]][[3]], maxK)
      ml[[k]] <- fun(parRes[[r]][[k]], ml[[k]], parRes[[r]][[1]][[3]], maxK)

      #mCount <- connectivityMatrix(rep(1, length(parRes[[r]][[1]][[3]])), mCount, parRes[[r]][[1]][[3]])
      #ml[[k]] <- connectivityMatrix(parRes[[r]][[k]], ml[[k]], parRes[[r]][[1]][[3]])
    }
}

src <- '
  Rcpp::List cas(a);
  Rcpp::List mat(b);
  Rcpp::

  Rcpp::NumericVector clusterAssignments(a);
  Rcpp::NumericMatrix bm(b);
  arma::mat m = Rcpp::as<arma::mat>(bm);
  Rcpp::NumericVector sampleKey(c);
  Rcpp::NumericVector kmax(d);
  int kint = kmax[0];

  Rcpp::List cls(kint); // list of binary vectors of sample memebership
  int n = bm.ncol(); // number of samples

  for (int ki=1; ki < kint; ++ki) {      // for each cluster assignment possible
    Rcpp::NumericVector x(n);            //    vector for each sample in n
    for (int si=0; si < clusterAssignments.size(); ++si) {    // for each sample
      if (clusterAssignments[si] == ki) {                     // if the sample is in the cluster
        x[(sampleKey[si])-1] = 1;
      } else {
        x[(sampleKey[si])-1] = 0;
      }
    }
    cls[ki] = x;
  }

  for (int ki=2; ki <= kint; ++ki) { // for each cluster
    arma::mat xx(n,n,arma::fill::zeros);  // empty matix
    IntegerVector a(cls[(ki-1)]);
    for (int i=0; i < n; ++i) {
      for (int j=0; j < n; ++j) {
        xx(i,j) = a[i] * a[j];
      }
    }
    m = m+xx;
  }
  return(Rcpp::wrap(m));
  '
fun <- cxxfunction(signature(a = "numeric", b = "matrix", c = "numeric", d = "numeric"), src, plugin = "RcppArmadillo")

fun(c(1, 1, 1, 1, 1, 1, 2, 1), matrix(c(0), 10, 10), c(1, 2, 4, 5, 6, 7, 8, 9), 3)


connectivityMatrix <- function(clusterAssignments, m, sampleKey) {
    ## input: named vector of cluster assignments, matrix to add connectivities
    ## output: connectivity matrix

    # have a list of cluster assignments and a list of sample IDs
    names(clusterAssignments) <- sampleKey

    # for each possible cluster assignment
    #   make a vector of samples that are in that cluster
    cls <- lapply(unique(clusterAssignments), function(i) as.numeric(names(clusterAssignments[clusterAssignments %in% i])))  #list samples by clusterId

    # for each cluster
    for (i in 1:length(cls)) {

        nelts <- 1:ncol(m)
        cl <- as.numeric(nelts %in% cls[[i]])  ## produces a binary vector
        updt <- outer(cl, cl)  #product of arrays with * function; with above indicator (1/0) statement updates all cells to indicate the sample pair was observed int the same cluster;
        m <- m + updt
    }
    return(m)
}

connectivityMatrix(c(1, 1, 1, 1, 1, 1, 2, 1), matrix(c(0), 10, 10), c(1, 2, 4, 5, 6, 7, 8, 9))



library(ConsensusClusterParallelPlus)
mmm <- matrix(rnorm(1000), ncol=100)
mdist <- dist(x=mmm)
res <- ConsensusClusterParallelPlusTest(d=mdist, cores=4, returnML=T)



require(Rcpp)
require(inline)
require(RcppArmadillo)
src <- '
  Rcpp::List a(x);
  int n = a.size();
  for (int i=0; i<n; ++i) {
    NumericVector b(a[i]);
    Rcout << b << arma::endl;
  }
  '
bigfun <- cxxfunction(signature(x = "list"), src, plugin = "RcppArmadillo")

bigfun(list(c(1, 1, 1, 1, 1, 1, 2, 1), c(1,2,3), c(1, 2, 4, 5, 6, 7, 8, 9)))
