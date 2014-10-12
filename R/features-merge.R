##' Merge features, typically after feature prediction in multiple samples.
##'
##' Merged features are the union of splice junctions and internal exons.
##' For terminal exons with shared spliced boundary, the longest exon is
##' retained.
##' 
##' @title Merge redundant features
##' @param ... one or more \code{TxFeatures} objects, or a single
##'   list of \code{TxFeatures} objects
##' @param min_n_sample Minimum number of samples a feature must be
##'   observed in to be included
##' @return \code{TxFeatures} object with merged features
##' @examples
##' txf_merged <- mergeTxFeatures(txf, txf)
##' @author Leonard Goldstein

mergeTxFeatures <- function(..., min_n_sample = 1)
{

    dots <- list(...)

    if (length(dots) == 1 && is(dots[[1]], "list")) {

        x <- dots[[1]]

    } else {

        x <- dots

    }
    
    if (all(sapply(x, is, "TxFeatures") | sapply(x, is.null)) &&
        !all(sapply(x, is.null))) {

        x <- x[which(elementLengths(x) > 0)]
        x <- lapply(x, dropMcols)
        features <- do.call(c, x)
        
    } else {

        stop("... must be one or more TxFeatures or a single
            list of TxFeatures")

    }

    if (is.null(features)) { return(TxFeatures()) }
    
    ## obtain unique splice junctions and internal exons
    J <- selectFeatures(features, "J", min_n_sample)
    I <- selectFeatures(features, "I", min_n_sample)
    
    ## merge terminal exons
    Z <- mergeExonsTerminal(features, min_n_sample)

    ## combine features
    features <- c(J, I, Z)

    if (is.null(features)) { return(TxFeatures()) }

    ## filter isolated exons (without flanking splice junctions)
    if (min_n_sample > 1) { features <- removeExonsIsolated(features) }
            
    ## sort features
    features <- sort(features)

    return(features)

}

selectFeatures <- function(features, type, min_n_sample = 1)
{

    index <- which(type(features) == type)

    if (length(index) > 0) {

        features <- features[index]
        co <- gr2co(features)
        co_n <- table(co)
        i <- which(co %in% names(which(co_n >= min_n_sample)))
        features <- unique(features[i])
        
    } else {

        si <- seqinfo(features)
        features <- TxFeatures()
        seqinfo(features) <- si
        
    }
    
    return(features)

}

mergeExonsTerminal <- function(features, min_n_sample = 1)
{

    index <- which(type(features) %in% c("F", "L"))

    if (length(index) > 0) {
    
        features <- features[index]
        splicesite <- feature2name(features, collapse_terminal = TRUE)
        splicesite_n <- table(splicesite)
        i <- which(splicesite %in% names(which(splicesite_n >= min_n_sample)))
        features <- features[i]
        splicesite <- splicesite[i]
        splicesite <- factor(splicesite)
        splicesite_i <- split(seq_along(features), splicesite)
        splicesite_w <- split(width(features), splicesite)
        splicesite_i <- mapply(function(i, w) { i[which.max(w)] },
            i = splicesite_i, w = splicesite_w, SIMPLIFY = TRUE)
        exons <- features[splicesite_i]

        for (ann in c("txName", "geneName")) {

            exons_ann <- collapseCharacterList(slot(features, ann), splicesite)
            slot(exons, ann) <- setNames(exons_ann, NULL)
            
        }
            
    } else {

        exons <- TxFeatures()
        seqinfo(exons) <- seqinfo(features)
        
    }
    
    return(exons)

}

mergeSGFeatures <- function(...)
{

    dots <- list(...)

    if (length(dots) == 1 && is(dots[[1]], "list")) {

        x <- dots[[1]]

    } else {

        x <- dots

    }
    
    if (all(sapply(x, is, "SGFeatures") | sapply(x, is.null)) &&
        !all(sapply(x, is.null))) {

        x <- x[which(elementLengths(x) > 0)]
        x <- lapply(x, dropMcols)
        x <- do.call(c, x)
        
    } else {

        stop("... must be one or more SGFeatures objects or a single
            list of SGFeatures objects")

    }

    features <- granges(x)
    mcols(features)$type <- as.character(type(x))

    junctions <- features[mcols(features)$type == "J"]
    exons <- features[mcols(features)$type == "E"]

    features <- processFeatures(exons, junctions)
    features <- addFeatureID(features)
    features <- addGeneID(features)
    features <- SGFeatures(features)
    
    return(features)
    
}
