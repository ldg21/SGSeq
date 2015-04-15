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
##' @inheritParams processTerminalExons
##' @return \code{TxFeatures} object with filtered features
##' @keywords internal
##' @author Leonard Goldstein

removeExonsIsolated <- function(features)
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

    filtered <- features[!seq_along(features) %in% remove]

    return(filtered)
    
}

##' Predicted terminal exons are processed as described under Details.
##' 
##' Processing of terminal exon predictions is done in two steps:
##' (1) terminal exons that share a splice site with an internal exon are
##' filtered, and (2) remaining terminal exons that overlap other exons
##' are trimmed.
##' 
##' \code{predictTxFeatures} predicts flanking terminal exons for each
##' identified splice junction. This ensures that each splice junction
##' has a flanking exon after merging with \code{mergeTxFeatures}.
##' This approach results in many predicted terminal exons that
##' share a splice site with predicted internal exons (often contained
##' within them or with a short overhang due to incorrect alignments).
##' Most of these are not real terminal exons and are filtered before
##' further analysis. Filtering based on the overhang is controlled with
##' argument \code{min_overhang}.
##'
##' Some of the remaining predicted terminal exons overlap other exons
##' such that their unspliced boundary shows a short overhang with
##' respect to a spliced boundary of the overlapping exon. Often these
##' exon extensions into an intron are due to incorrect alignments.
##' Terminal exons with overhang smaller than \code{min_overhang} are
##' trimmed such that their trimmmed unspliced boundary coincides with
##' the spliced boundary of the overlapping exon. 
##' 
##' @title Process predicted terminal exons
##' @param features \code{TxFeatures} object
##' @param min_overhang Minimum overhang required to suppress filtering or
##'   trimming of predicted terminal exons (see Details). Use \code{NA} to
##'   exclude all terminal exons sharing a splice with an internal exon
##'   and trim all remaining terminal exons overlapping other exons.
##' @return \code{TxFeatures} object with processed features
##' @examples
##' txf_processed <- processTerminalExons(txf)
##' @author Leonard Goldstein

processTerminalExons <- function(features, min_overhang = NA)
{

    ## (1) filter terminal exons that share a splice site with an internal exon
  
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

    filtered <- features[!seq_along(features) %in% remove]

    ## (2) trim terminal exons overlapping other exons
    
    processed <- filtered
  
    for (type in c("F", "L")) {
    
        i_Z <- which(type(filtered) == type)
        f_Z <- filtered[i_Z]

        type_2 <- c("I", switch(type, "F" = "L", "L" = "F"))
        f_X <- filtered[type(filtered) %in% type_2]
        
        if (length(f_Z) > 0 && length(f_X) > 0) {

            hits <- findOverlaps(f_Z, f_X)

            start <- switch(type, "F" = TRUE, "L" = FALSE)
            u_Z <- flank(f_Z, -1, start)
            s_X <- flank(f_X, -1, start)

            hits_2 <- findOverlaps(u_Z[queryHits(hits)],
                f_X[subjectHits(hits)])
            hits_2 <- hits_2[queryHits(hits_2) == subjectHits(hits_2)]
            hits <- hits[!seq_along(hits) %in% queryHits(hits_2)]

            if (length(hits) > 0) {
            
                w <- abs(start(u_Z)[queryHits(hits)] -
                    start(s_X)[subjectHits(hits)])
                
                q_min_w <- tapply(w, queryHits(hits), min)

                if (!is.na(min_overhang)) {
                
                    q <- as.integer(names(which(q_min_w < min_overhang)))
                
                } else {
                
                    q <- unique(queryHits(hits))

                }

                processed[i_Z][q] <- flank(filtered[i_Z][q],
                    -1 * width(filtered[i_Z][q]) + q_min_w[as.character(q)],
                    !start)

            }
                
        }

    }

    return(processed)

}
