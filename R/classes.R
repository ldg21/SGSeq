## Validity checks

validTxFeatures <- function(object)
{

    slot_values <- list(
        strand = c("+", "-"),
        type = typeLevels("TxFeatures"))
    validSlotValues(object, slot_values)

}

validSGFeatures <- function(object)
{

    slot_values <- list(
        strand = c("+", "-"),
        type = typeLevels("SGFeatures"))
    validSlotValues(object, slot_values)

}

validTxSegments <- function(object)
{

    mcol_type <- c(
        from = "character",
        to = "character",
        type = "character",
        splice5p = "logical",
        splice3p = "logical",
        featureID = "character",
        segmentID = "integer",
        geneID = "integer",
        txName = "CharacterList",
        geneName = "CharacterList")
    validMcols(object, mcol_type)
    
}

validTxVariants <- function(object)
{

    mcol_type <- c(
        from = "character",
        to = "character",
        type = "character",
        featureID = "character",
        segmentID = "character",
        geneID = "integer",
        eventID = "integer",
        variantID = "integer",
        closed3p = "logical",
        closed5p = "logical",
        featureID5p = "IntegerList",
        featureID3p = "IntegerList",
        variantType = "CharacterList",
        variantName = "character",
        txName = "CharacterList",
        geneName = "CharacterList")
    validMcols(object, mcol_type)

}

validSGFeatureCounts <- function(object)
{

    assay_type <- c(
        counts = "integer",
        FPKM = "numeric")
    validAssays(object, assay_type)
    
}

validTxVariantCounts <- function(object)
{

    assay_type <- c(
        countsVariant5p = "integer",
        countsVariant3p = "integer",
        countsTotal5p = "integer",
        countsTotal3p = "integer",
        variantFreq = "numeric")
    validAssays(object, assay_type)
    
}

## Validity checks - helper functions

validSlotValues <- function(object, slot_values)
{

    specified <- names(slot_values)
    slots <- lapply(specified, slot, object = object)
    i <- which(!mapply(validValues, slots, slot_values))

    if (length(i) > 0) {

        return(paste("invalid values for slot(s)",
            paste(dQuote(specified[i]), collapse = ",")))

    }
    
    return(TRUE)

}

validValues <- function(x, values)
{

    return(all(x %in% values))

}

typeLevels <- function(class)
{

    switch(class,
        TxFeatures = c("J", "I", "F", "L", "U"),
        SGFeatures = c("J", "E", "D", "A"))
    
}

validMcols <- function(object, mcol_type)
{

    i <- which(!names(mcol_type) %in% names(mcols(object)))
    
    if (length(i) > 0) {

        return(paste("missing elementMetadata column(s)",
            paste(dQuote(names(mcol_type)[i]), collapse = ",")))
        
    }
    
    i <- which(!mapply(is, mcols(object)[names(mcol_type)], mcol_type))
    
    if (length(i) > 0) {

        return(paste("invalid type for elementMetadata column(s)",
            paste(dQuote(names(mcol_type)[i]), collapse = ",")))
        
    }

    return(TRUE)

}

validAssays <- function(object, assay_type)
{

    i <- which(!names(assay_type) %in% names(assays(object)))
    
    if (length(i) > 0) {

        return(paste("missing assay(s)",
            paste(dQuote(names(assay_type)[i]), collapse = ",")))
        
    }
    
    i <- which(!mapply(function(x, t) { is(as.vector(x), t) },
        assays(object)[names(assay_type)], assay_type))
    
    if (length(i) > 0) {

        return(paste("invalid type for assay(s)",
            paste(dQuote(names(assay_type)[i]), collapse = ",")))
        
    }

    return(TRUE)

}

## Class definitions

setClass(
    Class = "TxFeatures",
    slots = c(
        type = "factor",
        txName = "CharacterList",
        geneName = "CharacterList"),
    contains = "GRanges",
    validity = validTxFeatures
)

setClass(
    Class = "SGFeatures",
    slots = c(
        type = "factor",
        splice5p = "logical",
        splice3p = "logical",
        featureID = "integer",
        geneID = "integer",
        txName = "CharacterList",
        geneName = "CharacterList"),
    contains = "GRanges",
    validity = validSGFeatures
)

setClass(
     Class = "TxSegments",
     slots = c(unlistData = "SGFeatures"),
     contains = "GRangesList",
     validity = validTxSegments
)

setClass(
    Class = "TxVariants",
    slots = c(unlistData = "SGFeatures"),
    contains = "GRangesList",
    validity = validTxVariants
)

setClass(
    Class = "SGFeatureCounts",
    slots = c(rowData = "SGFeatures"),
    contains = "SummarizedExperiment",
    validity = validSGFeatureCounts
)

setClass(
    Class = "TxVariantCounts",
    slots = c(rowData = "TxVariants"),
    contains = "SummarizedExperiment",
    validity = validTxVariantCounts
)

setClassUnion("Features", c("TxFeatures", "SGFeatures"))

setClassUnion("Paths", c("TxSegments", "TxVariants"))

setClassUnion("Counts", c("SGFeatureCounts", "TxVariantCounts"))

## Methods for extraColumnSlotNames

setMethod(GenomicRanges:::extraColumnSlotNames, "TxFeatures",
    function(x) { c("type", "txName", "geneName") })

setMethod(GenomicRanges:::extraColumnSlotNames, "SGFeatures",
    function(x) { c("type", "splice5p", "splice3p", "featureID", "geneID",
        "txName", "geneName") })

## Constructor functions

##' Creates an instance of S4 class \code{TxFeatures} for storing
##' transcript features.
##' 
##' \code{TxFeatures} extend \code{GRanges} with column slot \code{type}
##' specifying feature type. \code{type} is a factor with levels
##' \code{J} (splice junction), \code{I} (internal exon),
##' \code{F} (5' terminal exon), \code{L} (3' terminal exon),
##' \code{U} (unspliced transcript).
##'
##' \code{txName} and \code{geneName} are CharacterLists storing
##' transcript and gene annotation, respectively.
##' 
##' @title Constructor function for S4 class \code{TxFeatures}
##' @param x \code{GRanges} with known strand (\dQuote{+}, \dQuote{-})
##' @param type Character vector or factor taking values in
##'   \code{J}, \code{I}, \code{F}, \code{L}, \code{U}
##' @param txName \code{CharacterList} of transcript names or \code{NULL}
##' @param geneName \code{CharacterList} of gene names or \code{NULL}
##' @return A \code{TxFeatures} object
##' @examples
##' gr <- GRanges(1, IRanges(101, 200), "+")
##' txf <- TxFeatures(gr, type = "J")
##' @author Leonard Goldstein

TxFeatures <- function(x, type = mcols(x)$type, txName = mcols(x)$txName,
    geneName = mcols(x)$geneName)
{

    if (missing(x)) {

        x <- GRanges()
        type <- factor(levels = typeLevels("TxFeatures"))
        txName <- CharacterList()
        geneName <- CharacterList()

    } else {
    
        if (is.character(type))
            type <- factor(type, levels = typeLevels("TxFeatures"))
        
        if (is.null(txName))
            txName <- CharacterList(vector("list", length(x)))
        
        if (is.null(geneName))
            geneName <- CharacterList(vector("list", length(x)))

    }
        
    new("TxFeatures", granges(x), type = type, txName = txName,
        geneName = geneName)

}

##' Creates an instance of S4 class \code{SGFeatures} for storing
##' splice graph features.
##' 
##' \code{SGFeatures} extend \code{GRanges} with column slot \code{type}
##' specifying feature type. \code{type} is a factor with levels
##' \code{J} (splice junction), \code{E} (exon bin),
##' \code{D} (splice donor), \code{A} (splice acceptor).
##'
##' \code{splice5p} and \code{splice3p} are logical vectors indicating
##' whether reads extending across the 5' and 3' boundaries of an
##' exon bin must be spliced at the boundary to be considered compatible
##' with the exon bin.
##'
##' \code{featureID} and \code{geneID} are integer vectors representing
##' unique identifiers for features and genes (connected components in
##' the splice graph).
##' 
##' \code{txName} and \code{geneName} are CharacterLists storing
##' transcript and gene annotation, respectively.
##' 
##' @title Constructor function for S4 class \code{SGFeatures}
##' @param x \code{GRanges} with known strand (\dQuote{+}, \dQuote{-})
##' @param type Character vector or factor taking values in \code{J},
##'   \code{E}, \code{D}, \code{A}
##' @param splice5p Logical vector indicating whether reads extending across
##'   the 5' boundary of an exon bin must be spliced at the boundary
##' @param splice3p Logical vector indicating whether reads extending across
##'   the 3' boundary of an exon bin must be spliced at the boundary
##' @param featureID Integer vector of feature IDs
##' @param geneID Integer vector of gene IDs
##' @param txName \code{CharacterList} of transcript names or \code{NULL}
##' @param geneName \code{CharacterList} of gene names or \code{NULL}
##' @return An \code{SGFeatures} object
##' @examples
##' sgf <- SGFeatures()
##' @author Leonard Goldstein

SGFeatures <- function(x, type = mcols(x)$type,
    splice5p = mcols(x)$splice5p, splice3p = mcols(x)$splice3p,
    featureID = mcols(x)$featureID, geneID = mcols(x)$geneID, 
    txName = mcols(x)$txName, geneName = mcols(x)$geneName)
{

    if (missing(x)) {

        x <- GRanges()
        type <- factor(levels = typeLevels("SGFeatures"))
        splice5p <- logical()
        splice3p <- logical()
        featureID <- integer()
        geneID <- integer()
        txName <- CharacterList()
        geneName <- CharacterList()
        
    } else {    
    
        if (is.character(type))
            type <- factor(type, levels = typeLevels("SGFeatures"))
        
        if (is.null(txName))
            txName <- CharacterList(vector("list", length(x)))
        
        if (is.null(geneName))
            geneName <- CharacterList(vector("list", length(x)))

    }
    
    new("SGFeatures", granges(x), type = type, splice5p = splice5p,
        splice3p = splice3p, featureID = featureID, geneID = geneID,
        txName = txName, geneName = geneName)
        
}

##' Creates an instance of S4 class \code{TxSegments} for storing
##' transcript segments. 
##' 
##' @title Constructor function for S4 class \code{TxSegments}
##' @param x \code{GRangesList} of \code{SGFeatures} with appropriate
##'   outer elementMetadata columns
##' @return A \code{TxSegments} object
##' @keywords internal
##' @author Leonard Goldstein

TxSegments <- function(x)
{

    if (missing(x)) {

        x <- GRangesList()
        x@unlistData <- SGFeatures()
        x@elementMetadata <- DataFrame(
            from = character(),
            to = character(),
            type = character(),
            splice5p = logical(),
            splice3p = logical(),
            featureID = character(),
            segmentID = integer(),
            geneID = integer(),
            txName = CharacterList(),
            geneName = CharacterList())

    } else {
        
        n <- length(x)
        
        if (is.null(mcols(x)$txName))
            mcols(x)$txName <- CharacterList(vector("list", n))
    
        if (is.null(mcols(x)$geneName))
            mcols(x)$geneName <- CharacterList(vector("list", n))

    }
    
    new("TxSegments", x)

}

##' Creates an instance of S4 class \code{TxVariants} for storing
##' transcript variants. 
##' 
##' @title Constructor function for S4 class \code{TxVariants}
##' @param x \code{GRangesList} of \code{SGFeatures} with appropriate
##'   outer elementMetadata columns
##' @return A \code{TxVariants} object
##' @examples
##' txv <- TxVariants()
##' @author Leonard Goldstein

TxVariants <- function(x)
{

    if (missing(x)) {

        x <- GRangesList()
        x@unlistData <- SGFeatures()
        x@elementMetadata <- DataFrame(
            from = character(),
            to = character(),
            type = character(),
            featureID = character(),
            segmentID = character(),
            geneID = integer(),
            eventID = integer(),
            variantID = integer(),
            closed3p = logical(),
            closed5p = logical(),
            featureID5p = IntegerList(),
            featureID3p = IntegerList(),
            variantType = CharacterList(),
            variantName = character(),
            txName = CharacterList(),
            geneName = CharacterList())

    } else {
    
        n <- length(x)
        
        if (is.null(mcols(x)$txName))
            mcols(x)$txName <- CharacterList(vector("list", n))
        
        if (is.null(mcols(x)$geneName))
            mcols(x)$geneName <- CharacterList(vector("list", n))
        
        if (is.null(mcols(x)$variantType))
            mcols(x)$variantType <- CharacterList(vector("list", n))
        
        if (is.null(mcols(x)$variantName))
            mcols(x)$variantName <- vector("character", n)

    }

    new("TxVariants", x)

}

##' Creates an instance of S4 class \code{SGFeatureCounts} for storing
##' compatible splice graph feature counts.
##' 
##' @title Constructor function for S4 class \code{SGFeatureCounts}
##' @param x \code{SummarizedExperiment} with \code{SGFeatures}
##'   as \code{rowData} and assays \dQuote{counts}, \dQuote{FPKM}
##' @return An \code{SGFeatureCounts} object
##' @examples
##' sgfc <- SGFeatureCounts()
##' @author Leonard Goldstein

SGFeatureCounts <- function(x)
{

    if (missing(x)) {

        assays <- list(
            counts = matrix(integer(), 0, 0),
            FPKM = matrix(numeric(), 0, 0))

        x <- SummarizedExperiment(assays, rowData = SGFeatures())
        
    }

    new("SGFeatureCounts", x)

}

##' Creates an instance of S4 class \code{TxVariantCounts} for storing
##' representative transcript variant counts.
##' 
##' @title Constructor function for S4 class \code{SGFeatureCounts}
##' @param x \code{SummarizedExperiment} with \code{TxVariants}
##'   as \code{rowData} and appropriate assays
##' @return A \code{TxVariantCounts} object
##' @examples
##' txvc <- TxVariantCounts()
##' @author Leonard Goldstein

TxVariantCounts <- function(x)
{

    if (missing(x)) {

        assays <- list(
            countsVariant5p = matrix(integer(), 0, 0),
            countsVariant3p = matrix(integer(), 0, 0),
            countsTotal5p = matrix(integer(), 0, 0),
            countsTotal3p = matrix(integer(), 0, 0),
            variantFreq = matrix(numeric(), 0, 0))

        x <- SummarizedExperiment(assays, rowData = TxVariants())

    }
    
    new("TxVariantCounts", x)

}
