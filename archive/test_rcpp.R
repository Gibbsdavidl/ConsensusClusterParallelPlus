

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
  Rcpp::List parRes(a);            // results from the parallel runs
  Rcpp::NumericMatrix mCount(b);   // mCount matrix
  Rcpp::List ml(c);                // list of ml matrices
  Rcpp::NumericVector kmax(d);     // the kmax
  int kint = kmax[0];
  Rcpp::NumericVector rc(e);
  int repCount = rc[0];

  for (int r=0; r < repCount; ++r) {    // for each rep
      for (int k=1; r < kint; ++k) {      // for each k

        Rcpp::List extract1(parRes[r]);     // extract from parRes
        Rcpp::List extract2(extract1[0]);
        Rcpp::NumericVector clusterAssignments(extract2[k]);
        Rcpp::NumericVector sampleKey(extract2[3]);
        Rcpp::NumericVector mCountAssignments(clusterAssignments.size(), 1.0);

        Rcpp::NumericMatrix bm = ml[k];                       // extract from ml
        int n = bm.ncol(); // number of samples
        arma::mat mlmat = Rcpp::as<arma::mat>(bm);
        arma::mat mcmat = Rcpp::as<arma::mat>(mCount);

        Rcpp::List cls_mCount(kint); // list of binary vectors of sample memebership
        Rcpp::List cls_ml(kint);     // list of binary vectors of sample memebership

        // for mCount number of possible clusters is always 1
        for (int si=0; si < clusterAssignments.size(); ++si) {    // for each sample
            Rcpp::NumericVector y(n);                                //    vector for each sample in n
            if (mCountAssignments[si] == 1) {                       // if the sample is in the cluster
                y[(sampleKey[si])-1] = 1;
            } else {
                y[(sampleKey[si])-1] = 0;
            }
            cls_mCount[0] = y;
        }

        // for ml
        int thisKMax = max(clusterAssignments);
        for (int ki=1; ki < thisKMax; ++ki) {      // for each cluster assignment possible
            Rcpp::NumericVector x(n);               //    vector for each sample in n
            for (int si=0; si < clusterAssignments.size(); ++si) {    // for each sample
                if (clusterAssignments[si] == ki) {                     // if the sample is in the cluster
                    x[(sampleKey[si])-1] = 1;
                } else {
                    x[(sampleKey[si])-1] = 0;
                }
            }
            cls_ml[ki] = x;
        }

        arma::mat xx(n,n,arma::fill::zeros);  // empty matix
        IntegerVector a(cls_mCount);          // co-selected
        for (int i=0; i < n; ++i) {
            for (int j=0; j < n; ++j) {
                xx(i,j) = a[i] * a[j];
            }
        }
        mcmat = mcmat+xx;

        for (int ki=2; ki <= kint; ++ki) { // for each cluster
            arma::mat xx(n,n,arma::fill::zeros);  // empty matix
            IntegerVector a(cls_ml[(ki-1)]);
            for (int i=0; i < n; ++i) {
                for (int j=0; j < n; ++j) {
                    xx(i,j) = a[i] * a[j];
                }
            }
            mlmat = mlmat+xx;
        }
  return Rcpp::List::create(Rcpp::NumericMatrix(mcmat), Rcpp::NumericMatrix(mlmat));
  '
fun <- cxxfunction(signature(a = "list", b = "matrix", c = "list", d = "numeric", e = "numeric"), src, plugin = "RcppArmadillo")

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
      List b(a[0]);
      List c(a[1]);
      NumericVector d(b[0]);
      NumericVector e(c[0]);
      arma::vec x(d);
      arma::vec y(e);
      return List::create(d, e);
      '
    bigfun <- cxxfunction(signature(x = "list"), src, plugin = "RcppArmadillo")

    bigfun(list(list(c(1, 1, 1, 1, 1, 1, 2, 1), c(1,2,3)), list(c(1, 2, 4, 5, 6, 7, 8, 9), c(1,1,1))))
