ConsensusClusterParallelPlus <- function(d = NULL, maxK = 3, reps = 10, cores = 1,
    pItem = 0.8, pFeature = 1, clusterAlg = "hc", title = "untitled_consensus_cluster",
    innerLinkage = "average", finalLinkage = "average", distance = "pearson", ml = NULL,
    tmyPal = NULL, seed = NULL, plot = NULL, writeTable = FALSE, weightsItem = NULL,
    weightsFeature = NULL, verbose = F, corUse = "everything", returnML=F) {
    ## description: runs consensus subsamples
    require(fastcluster)

    if (is.null(seed) == TRUE) {
        seed = timeSeed = as.numeric(Sys.time())
    }
    set.seed(seed)

    # distance=ifelse( inherits(d,'dist'), attr( d, 'method' ), 'pearson' )

    if (is.null(ml) == TRUE) {

        if (!class(d) %in% c("dist", "matrix", "ExpressionSet")) {
            stop("d must be a matrix, distance object or ExpressionSet (eset object)")
        }

        if (inherits(d, "dist")) {
            ## if d is a distance matrix, fix a few things so that they don't cause problems
            ## with the analysis Note, assumption is that if d is a distance matrix, the user
            ## doesn't want to sample over the row features
            if (is.null(attr(d, "method"))) {
                attr(d, "method") <- distance <- "unknown - user-specified"
            }
            if (is.null(distance) || (distance != attr(d, "method"))) {
                distance <- attr(d, "method")
            }

            if ((!is.null(pFeature)) && (pFeature < 1)) {
                message("Cannot use the pFeatures parameter when specifying a distance matrix as the data object\n")
                pFeature <- 1
            }
            if (!is.null(weightsFeature)) {
                message("Cannot use the weightsFeature parameter when specifying a distance matrix as the data object\n")
                weightsFeature <- NULL
            }
            if (clusterAlg == "km") {
                message("Note: k-means will cluster the distance matrix you provided.  This is similar to kmdist option when suppling a data matrix")
                ## d <- as.matrix( d ) #this is now done w/in ccRun
            }
        } else {
            if (is.null(distance)) {
                ## we should never get here, but just in case
                distance <- "pearson"
            }
        }

        if ((clusterAlg == "km") && inherits(distance, "character") && (distance !=
            "euclidean")) {
            message("Note: The km (kmeans) option only supports a euclidean distance metric when supplying a data matrix.  If you want to cluster a distance matrix using k-means use the 'kmdist' option, or use a different algorithm such as 'hc' or 'pam'.  Changing distance to euclidean")
            distance <- "euclidean"
        }


        if (inherits(d, "ExpressionSet")) {
            d <- exprs(d)
        }

        ml <- ccRunPar(d = d, maxK = maxK, repCount = reps, coreCount = cores, diss = inherits(d,
            "dist"), pItem = pItem, pFeature = pFeature, innerLinkage = innerLinkage,
            clusterAlg = clusterAlg, weightsFeature = weightsFeature, weightsItem = weightsItem,
            distance = distance, verbose = verbose, corUse = corUse, returnML = returnML)
        if (returnML == T) {
          return(ml)
        }
    }
    res = list()

    ## make results directory
    if ((is.null(plot) == FALSE | writeTable) & !file.exists(paste(title, sep = ""))) {
        dir.create(paste(title, sep = ""))
    }

    ## write log file
    log <- matrix(ncol = 2, byrow = T, c("title", title, "maxK", maxK, "input matrix rows",
        ifelse(inherits(d, "matrix"), nrow(d), "dist-mat"), "input matrix columns",
        ifelse(inherits(d, "matrix"), ncol(d), ncol(as.matrix(d))), "number of bootstraps",
        reps, "item subsampling proportion", pItem, "feature subsampling proportion",
        ifelse(is.null(pFeature), 1, pFeature), "cluster algorithm", clusterAlg,
        "inner linkage type", innerLinkage, "final linkage type", finalLinkage, "correlation method",
        distance, "plot", if (is.null(plot)) NA else plot, "seed", if (is.null(seed)) NA else seed))
    colnames(log) = c("argument", "value")
    if (writeTable) {
        write.csv(file = paste(title, "/", title, ".log.csv", sep = ""), log, row.names = F)
    }
    if (is.null(plot)) {
        ## nothing
    } else if (plot == "pngBMP") {
        bitmap(paste(title, "/", "consensus%03d.png", sep = ""))
    } else if (plot == "png") {
        png(paste(title, "/", "consensus%03d.png", sep = ""))

    } else if (plot == "pdf") {
        pdf(onefile = TRUE, paste(title, "/", "consensus.pdf", sep = ""))
    } else if (plot == "ps") {
        postscript(onefile = TRUE, paste(title, "/", "consensus.ps", sep = ""))
    }

    colorList = list()
    colorM = rbind()  #matrix of colors.

    # 18 colors for marking different clusters
    thisPal <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
        "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#bd18ea",
        "#2ef4ca", "#f4cced", "#f4cc03", "#05188a", "#e5a25a", "#06f106", "#85848f",
        "#000000", "#076f25", "#93cd7f", "#4d0776", "#ffffff")

    ## plot scale
    colBreaks = NA
    if (is.null(tmyPal) == TRUE) {
        colBreaks = 10
        tmyPal = myPal(colBreaks)
    } else {
        colBreaks = length(tmyPal)
    }
    sc = cbind(seq(0, 1, by = 1/(colBreaks)))
    rownames(sc) = sc[, 1]
    sc = cbind(sc, sc)
    heatmap(sc, Colv = NA, Rowv = NA, symm = FALSE, scale = "none", col = tmyPal,
        na.rm = TRUE, labRow = rownames(sc), labCol = F, main = "consensus matrix legend")

    for (tk in 2:maxK) {
        if (verbose) {
            message(paste("consensus ", tk))
        }
        fm = ml[[tk]]
        hc = fastcluster::hclust(as.dist(1 - fm), method = finalLinkage)
        message("clustered")
        ct = cutree(hc, tk)
        names(ct) = colnames(d)
        if (class(d) == "dist") {
            names(ct) = colnames(as.matrix(d))
        }
        c = fm

        colorList = setClusterColors(res[[tk - 1]][[3]], ct, thisPal, colorList)
        pc = c
        pc = pc[hc$order, ]  #pc is matrix for plotting, same as c but is row-ordered and has names and extra row of zeros.

        if (!is.null(plot) && plot == "pngBMP") {
            pc = pc[, hc$order]  #mod for no tree
            pc = rbind(pc, 0)
            # no dendrogram if pngBMP
            oc = colorList[[1]][hc$order]  #mod for no tree
            heatmap(pc, Colv = NA, Rowv = NA, symm = FALSE, scale = "none", col = tmyPal,
                na.rm = TRUE, labRow = F, labCol = F, mar = c(5, 5), main = paste("consensus matrix k=",
                  tk, sep = ""), ColSideCol = oc)
        } else {
            pc = rbind(pc, 0)
            # former with tree:
            heatmap(pc, Colv = as.dendrogram(hc), Rowv = NA, symm = FALSE, scale = "none",
                col = tmyPal, na.rm = TRUE, labRow = F, labCol = F, mar = c(5, 5),
                main = paste("consensus matrix k=", tk, sep = ""), ColSideCol = colorList[[1]])
        }

        legend("topright", legend = unique(ct), fill = unique(colorList[[1]]), horiz = FALSE)

        res[[tk]] = list(consensusMatrix = c, consensusTree = hc, consensusClass = ct,
            ml = ml[[tk]], clrs = colorList)
        colorM = rbind(colorM, colorList[[1]])
    }
    CDF(ml)
    clusterTrackingPlot(colorM[, res[[length(res)]]$consensusTree$order])
    if (is.null(plot) == FALSE) {
        dev.off()
    }
    res[[1]] = colorM
    if (writeTable) {
        for (i in 2:length(res)) {
            write.csv(file = paste(title, "/", title, ".k=", i, ".consensusMatrix.csv",
                sep = ""), res[[i]]$consensusMatrix)
            write.table(file = paste(title, "/", title, ".k=", i, ".consensusClass.csv",
                sep = ""), res[[i]]$consensusClass, col.names = F, sep = ",")
        }
    }
    return(res)
}


ccRunPar <- function(d = d, maxK = NULL, repCount = NULL, coreCount = NULL, diss = inherits(d,
    "dist"), pItem = NULL, pFeature = NULL, innerLinkage = NULL, distance = NULL,
    clusterAlg = NULL, weightsItem = NULL, weightsFeature = NULL, verbose = NULL,
    corUse = NULL, returnML=F) {

    m = vector(mode = "list", repCount)
    ml = vector(mode = "list", maxK)
    n <- ifelse(diss, ncol(as.matrix(d)), ncol(d))

    ## mCount is possible number of times that two sample occur in same random sample,
    ## independent of k mCount stores number of times a sample pair was sampled
    ## together.
    require(Matrix)
    mCount = mConsist = Matrix(c(0), ncol = n, nrow = n, sparse = TRUE)

    ## ml
    ml[[1]] = c(0)
    for (k in 2:maxK) {
        ml[[k]] = mConsist  #initialize with list of sparse Matrix
    }

    if (is.null(distance))
        distance <- "euclidean"  ## necessary if d is a dist object and attr( d, 'method' ) == NULL

    acceptable.distance <- c("euclidean", "maximum", "manhattan", "canberra", "binary",
        "minkowski", "pearson", "spearman")

    main.dist.obj <- NULL
    if (diss) {
        main.dist.obj <- d
        ## reset the pFeature & weightsFeature params if they've been set (irrelevant if d
        ## is a dist matrix)
        if ((!is.null(pFeature)) && (pFeature < 1)) {
            message("user-supplied data is a distance matrix; ignoring user-specified pFeature parameter\n")
            pFeature <- 1  # set it to 1 to avoid problems with sampleCols
        }
        if (!is.null(weightsFeature)) {
            message("user-supplied data is a distance matrix; ignoring user-specified weightsFeature parameter\n")
            weightsFeature <- NULL  # set it to NULL to avoid problems with sampleCols
        }
    } else {
        ## d is a data matrix we're not sampling over the features
        if ((clusterAlg != "km") && (is.null(pFeature) || ((pFeature == 1) && is.null(weightsFeature)))) {
            ## only generate a main.dist.object IFF 1) d is a matrix, 2) we're not sampling
            ## the features, and 3) the algorithm isn't 'km'
            if (inherits(distance, "character")) {
                if (!distance %in% acceptable.distance & (class(try(get(distance),
                  silent = T)) != "function"))
                  stop("unsupported distance.")

                if (distance == "pearson" | distance == "spearman") {
                  main.dist.obj <- as.dist(1 - cor(d, method = distance, use = corUse))
                } else if (class(try(get(distance), silent = T)) == "function") {
                  main.dist.obj <- get(distance)(t(d))
                } else {
                  main.dist.obj <- dist(t(d), method = distance)
                }
                attr(main.dist.obj, "method") <- distance
            } else stop("unsupported distance specified.")
        } else {
            ## pFeature < 1 or a weightsFeature != NULL since d is a data matrix, the user
            ## wants to sample over the gene features, so main.dist.obj is left as NULL
        }
    }


    require(parallel)
    parsteps <- seq(from=1, to=repCount, by=coreCount)
    for (paridx in parsteps) {

        print(paridx)

        # do a set of clusterings, one for each core
        parRes <- mclapply(X = 1:coreCount, FUN = function(xi) {
            eachRep(xi, d, pItem, pFeature, weightsItem, weightsFeature, main.dist.obj,
                clusterAlg, distance, innerLinkage, maxK, verbose)
        }, mc.cores = coreCount)

        if (returnML) {
            return(parRes)
        }

        # then extract and sum up the list
        for (lp in 1:coreCount) {
            # sum up the mCount, over k, within the call to eachRep
            mCount <- parRes[[lp]][[1]] + mCount

            for (k in 2:maxK) {
                ml[[k]] <- ml[[k]] + parRes[[lp]][[k]]
            }
        }

        # then clear the memory
        parRes <- NULL
        gc()
    }

    # Convert sparse Matrices back to reg matrix
    mCount <- as.matrix(mCount)
    for (ki in 2:maxK) {
        ml[[ki]] <- as.matrix(ml[[ki]])
    }

    ## consensus fraction
    res = vector(mode = "list", maxK)
    for (k in 2:maxK) {
        ## fill in other half of matrix for tally and count.
        tmp = triangle(ml[[k]], mode = 3)
        tmpCount = triangle(mCount, mode = 3)
        res[[k]] = tmp/tmpCount
        res[[k]][which(tmpCount == 0)] = 0
    }
    message("end fraction")
    return(res)

}  # End ccRun


eachRep <- function(i, d, pItem, pFeature, weightsItem, weightsFeature, main.dist.obj,
    clusterAlg, distance, innerLinkage, maxK, verbose) {

    require(fastcluster)
    require(Matrix)

    #init
    resKs <- list()

    n <- ncol(d) ### What about distance matrices

    mCount <- Matrix(c(0), ncol = n, nrow = n, sparse=T)
    mConsist <- Matrix(c(0), ncol = n, nrow = n, sparse=T)
    ml <- vector(mode = "list", maxK)
    ml[[1]] <- mCount
    for (k in 2:maxK) {
        ml[[k]] <- mConsist  #initialize with list of sparse Matrix
    }

    ## take expression matrix sample, samples and genes
    sample_x = sampleCols(d, pItem, pFeature, weightsItem, weightsFeature)

    this_dist = NA
    if (!is.null(main.dist.obj)) {
        boot.cols <- sample_x$subcols
        this_dist <- as.matrix(main.dist.obj)[boot.cols, boot.cols]
        if (clusterAlg != "km") {
            ## if this isn't kmeans, then convert to a distance object
            this_dist <- as.dist(this_dist)
            attr(this_dist, "method") <- attr(main.dist.obj, "method")
        }
    } else {
        ## if main.dist.obj is NULL, then d is a data matrix, and either: 1) clusterAlg is
        ## 'km' 2) pFeatures < 1 or weightsFeatures have been specified, or 3) both so we
        ## can't use a main distance object and for every iteration, we will have to
        ## re-calculate either 1) the distance matrix (because we're also sampling the
        ## features as well), or 2) the submat (if using km)

        if (clusterAlg != "km") {
            if (!distance %in% acceptable.distance & (class(try(get(distance), silent = T)) !=
                "function"))
                stop("unsupported distance.")
            if ((class(try(get(distance), silent = T)) == "function")) {
                this_dist <- get(distance)(t(sample_x$submat))
            } else {
                if (distance == "pearson" | distance == "spearman") {
                  this_dist <- as.dist(1 - cor(sample_x$submat, use = corUse, method = distance))
                } else {
                  this_dist <- dist(t(sample_x$submat), method = distance)
                }
            }
            attr(this_dist, "method") <- distance
        } else {
            ## if we're not sampling the features, then grab the colslice
            if (is.null(pFeature) || ((pFeature == 1) && is.null(weightsFeature))) {
                this_dist <- d[, sample_x$subcols]
            } else {
                if (is.na(sample_x$submat)) {
                  stop("error submat is NA")
                }

                this_dist <- sample_x$submat
            }
        }
    }

    ## cluster samples for HC.
    this_cluster = NA
    if (clusterAlg == "hc") {
        this_cluster = fastcluster::hclust(this_dist, method = innerLinkage)
    }

    ## use samples for each k
    for (k in 2:maxK) {
        if (verbose) {
            message(paste("  k =", k))
        }
        if (i == 1) {
            ml[[k]] = mConsist  #initialize
        }
        this_assignment = NA
        if (clusterAlg == "hc") {
            ## prune to k for hc
            this_assignment = cutree(this_cluster, k)

        } else if (clusterAlg == "kmdist") {
            this_assignment = kmeans(this_dist, k, iter.max = 10, nstart = 1, algorithm = c("Hartigan-Wong"))$cluster

        } else if (clusterAlg == "km") {
            ## this_dist should now be a matrix corresponding to the result from sampleCols
            this_assignment <- kmeans(t(this_dist), k, iter.max = 10, nstart = 1,
                algorithm = c("Hartigan-Wong"))$cluster
        } else if (clusterAlg == "pam") {
            this_assignment <- pam(x = this_dist, k, diss = TRUE, metric = distance,
                cluster.only = TRUE)
        } else {
            ## optional cluterArg Hook.
            this_assignment <- get(clusterAlg)(this_dist, k)
        }
        # add to the tally was here
        resKs[[k]] <- this_assignment
    }  # End Ks

    samps <- sample_x[[3]]
    mc <- connectivityMatrix( rep( 1,length(samps)), mCount, samps )
    mCount <- mc

    for (ki in 2:maxK) {
        mi <- connectivityMatrix(resKs[[ki]], ml[[ki]], samps)
        ml[[ki]] <- mi
    }

    # fill up mCount matrix using resK[[1]] or sample_x
    # scale the mCount matrix by k
    #mCount[samps, samps] <- (maxK-1)
    ml[[1]] <- mCount
    return(ml)
}  # end eachRep


connectivityMatrixPar <- function( clusterAssignments, m, sampleKey){
  # m is a sparse matrix.
  ##input: named vector of cluster assignments, matrix to add connectivities
  ##output: connectivity matrix
  require(Matrix)
  names( clusterAssignments ) <- sampleKey
  cls <- lapply( unique( clusterAssignments ), function(i) as.numeric( names( clusterAssignments[ clusterAssignments %in% i ] ) ) )  #list samples by clusterId

  for ( i in 1:length( cls ) ) {
    nelts <- 1:ncol( m )
    cl <- as.numeric( nelts %in% cls[[i]] ) ## produces a binary vector
    updt <- outer( cl, cl ) #product of arrays with * function; with above indicator (1/0) statement updates all cells to indicate the sample pair was observed int the same cluster;
    updt <- Matrix(updt, sparse=T)
    m <- m + updt
  }
  return(m)
}
