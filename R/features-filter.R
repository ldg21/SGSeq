##' Filter previously predicted features using more stringent criteria.
##' 
##' Initial predictions with \code{predictTxFeatures} must have been performed
##' with \code{include_counts = TRUE} and \code{retain_coverage = TRUE},
##' so that predicted features contain elementMetadata columns \dQuote{N},
##' \dQuote{N_splicesite} and \dQuote{coverage}.
##'
##' @title Filter predicted features
##' @inheritParams predictTxFeaturesPerSample
##' @param features \code{TxFeatures} object with predicted features,
##'   including elementMetadata columns \dQuote{N}, \dQuote{N_splicesite} and
##'   \dQuote{coverage}.
##' @return \code{TxFeatures} object with filtered features
##' @keywords internal
##' @author Leonard Goldstein

filterFeatures <- function(features, paired_end, read_length, frag_length,
    lib_size, min_junction_count = NULL, alpha, psi, beta, gamma)
{

    if (is.null(min_junction_count)) {
    
        min_junction_count <- convertFpkmToCount(alpha, paired_end,
            read_length, frag_length, lib_size)

    }

    J <- filterJunctions(features, min_junction_count, psi)
    I <- filterExonsInternal(features, J, beta)
    X <- filterExonsTerminal(features, J, "F", gamma)
    Y <- filterExonsTerminal(features, J, "L", gamma)
    
    filtered <- c(J, I, X, Y)
    filtered <- sort(filtered)
    
    return(filtered)
    
}

filterJunctions <- function(features, min_junction_count, psi)
{

    J <- features[type(features) == "J"]
    J <- J[mcols(J)$N >= min_junction_count]
    J <- J[mcols(J)$N >= psi * max(mcols(J)$N_splicesite)]

    return(J)

}

filterExonsInternal <- function(features, junctions, beta)
{
    
    D <- flank(junctions, -1, TRUE)
    A <- flank(junctions, -1, FALSE)    
    I <- features[type(features) == "I"]
    I <- I[flank(I, -1, TRUE) %in% A & flank(I, -1, FALSE) %in% D]

    if (length(I) > 0) {
    
        I <- I[min(mcols(I)$coverage) >= beta * min(mcols(I)$N_splicesite)]

    }

    return(I)

}

filterExonsTerminal <- function(features, junctions, type = c("F", "L"),
    gamma)
{

    type <- match.arg(type)
    
    S <- flank(junctions, -1, switch(type, "F" = TRUE, "L" = FALSE))
    Z <- features[type(features) == type]
    Z <- Z[flank(Z, -1, switch(type, "F" = FALSE, "L" = TRUE)) %in% S]

    if (length(Z) == 0) { return(Z) }
    
    ranges <- mcols(Z)$coverage >= gamma * mcols(Z)$N_splicesite
    rl <- runLength(ranges)
    
    if (type == "F") {

        w <- plast(rl)
        Z <- flank(Z, -w, FALSE)
        ir <- IRanges(end = elementLengths(ranges), width = w)
        
    }
    if (type == "L") {

        w <- pfirst(rl)
        Z <- flank(Z, -w, TRUE)
        ir <- IRanges(start = 1, width = w)

    }

    i <- setNames(split(ir, seq_along(ir)), NULL)
    mcols(Z)$coverage <- mcols(Z)$coverage[i]

    return(Z)

}

##' @title Remove exons with no flanking splice junctions
##' @inheritParams filterTerminalExons
##' @return \code{TxFeatures} object with filtered features,
##'   or indices of retained features if \code{return_index = TRUE}
##' @keywords internal
##' @author Leonard Goldstein

removeExonsIsolated <- function(features, return_index = FALSE)
{

    i_J <- which(type(features) == "J")
    i_I <- which(type(features) == "I")
    i_X <- which(type(features) == "F")
    i_Y <- which(type(features) == "L")
    
    J <- features[i_J]
    I <- features[i_I]
    X <- features[i_X]
    Y <- features[i_Y]
    
    D <- flank(J, -1, TRUE)
    A <- flank(J, -1, FALSE)

    remove <- vector()

    if (length(i_X) > 0) {
    
        hits <- findOverlaps(flank(X, -1, FALSE), D)
        i <- which(!seq_along(i_X) %in% queryHits(hits))
        remove <- c(remove, i_X[i])
        
    }
    if (length(i_Y) > 0) {
    
        hits <- findOverlaps(flank(Y, -1, TRUE), A)
        i <- which(!seq_along(i_Y) %in% queryHits(hits))
        remove <- c(remove, i_Y[i])
        
    }
    if (length(i_I) > 0) {
        
        hits_A <- findOverlaps(flank(I, -1, TRUE), A)
        hits_D <- findOverlaps(flank(I, -1, FALSE), D)
        i <- which(!seq_along(i_I) %in%
            intersect(queryHits(hits_A), queryHits(hits_D)))
        remove <- c(remove, i_I[i])

    }

    index <- which(!seq_along(features) %in% remove)

    if (return_index) {

        return(index)

    } else {
                
        return(features[index])

    }
    
}

##' Terminal exons that share their splice site with an internal exon are
##' filtered based on the overhang with respect to the internal exon. 
##'
##' \code{predictTxFeatures} predicts flanking terminal exons for each
##' identified splice junction. Thus splice junctions are guaranteed to
##' have flanking exons, even after filtering exons during merging with
##' \code{mergeTxFeatures}. However, many predicted terminal exons share
##' a splice site with predicted internal exons and are often contained
##' within them. Many of these predictions are unlikely to be real terminal
##' exons and are excluded before further analysis.
##' 
##' @title Filter terminal exons that share splice sites with internal exons
##' @param features \code{TxFeatures} object
##' @param min_overhang For terminal exons sharing a splice site with an
##'   internal exon, minimum overhang required for terminal exon to be
##'   included. Use \code{NA} to remove all terminal exons sharing a splice
##'   site with an internal exon.
##' @param return_index Logical indicating whether indices of retained
##'   features should be returned instead of filtered features
##' @return \code{TxFeatures} object with filtered features,
##'   or indices of retained features if \code{return_index = TRUE}
##' @examples
##' txf_filtered <- filterTerminalExons(txf)
##' @author Leonard Goldstein

filterTerminalExons <- function(features, min_overhang = NA,
    return_index = FALSE)
{

    f_I <- features[type(features) == "I"]

    remove <- vector()
    
    for (type in c("F", "L")) {
    
        i_Z <- which(type(features) == type)
        f_Z <- features[i_Z]
        
        if (length(f_I) > 0 && length(f_Z) > 0) {

            start <- switch(type, "F" = FALSE, "L" = TRUE)
            s_Z <- flank(f_Z, -1, start = start)
            s_I <- flank(f_I, -1, start = start)
            hits <- findOverlaps(s_Z, s_I)
            
            if (!is.na(min_overhang) && length(hits) > 0) {
                
                w <- width(f_Z[queryHits(hits)]) -
                    width(f_I[subjectHits(hits)])
                q_min_w <- tapply(w, queryHits(hits), min)
                q <- as.integer(names(which(q_min_w < min_overhang)))
                
            } else {

                q <- unique(queryHits(hits))

            }

            remove <- c(remove, i_Z[q])
                        
        }

    }

    index <- which(!seq_along(features) %in% remove)

    if (return_index) {

        return(index)

    } else {
        
        return(features[index])

    }

}
