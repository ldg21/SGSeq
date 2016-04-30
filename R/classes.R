## Validity checks

validTxFeatures <- function(object)
{

    slot_values <- list(
        strand = c("+", "-"),
        type = typeLevels("TxFeatures"))
    c(validExtraColumnSlotLengths(object),
        validExtraColumnSlotValues(object, slot_values))

}

validSGFeatures <- function(object)
{

    slot_values <- list(
        strand = c("+", "-"),
        type = typeLevels("SGFeatures"))
    c(validExtraColumnSlotLengths(object),
        validExtraColumnSlotValues(object, slot_values))

}

validSGSegments <- function(object)
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

validSGVariants <- function(object)
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
        closed5p = "logical",
        closed3p = "logical",
        closed5pEvent = "logical",
        closed3pEvent = "logical",
        featureID5p = "IntegerList",
        featureID3p = "IntegerList",
        featureID5pEvent = "IntegerList",
        featureID3pEvent = "IntegerList",
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

validSGVariantCounts <- function(object)
{

    assay_type <- c(
        countsVariant5p = "integer",
        countsVariant3p = "integer",
        countsEvent5p = "integer",
        countsEvent3p = "integer",
        variantFreq = "numeric")
    validAssays(object, assay_type)

}

## Validity checks - helper functions

validExtraColumnSlotLengths <- function(object)
{

    l <- length(object)
    slots <- GenomicRanges:::extraColumnSlotNames(object)
    i <- which(elementNROWS(lapply(slots, slot, object = object)) != l)

    if (length(i) > 0) {

        return(paste("invalid length for slot(s)",
            paste(dQuote(slots[i]), collapse = ", ")))

    }

}

validExtraColumnSlotValues <- function(object, slot_values)
{

    specified <- names(slot_values)
    slots <- lapply(specified, slot, object = object)
    i <- which(!mapply(validValues, slots, slot_values))

    if (length(i) > 0) {

        return(paste("invalid values for slot(s)",
            paste(dQuote(specified[i]), collapse = ", ")))

    }

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

        return(paste("missing metadata column(s)",
            paste(dQuote(names(mcol_type)[i]), collapse = ", ")))

    }

    i <- which(!mapply(is, mcols(object)[names(mcol_type)], mcol_type))

    if (length(i) > 0) {

        return(paste("invalid type for metadata column(s)",
            paste(dQuote(names(mcol_type)[i]), collapse = ", ")))

    }

}

validAssays <- function(object, assay_type)
{

    i <- which(!names(assay_type) %in% names(assays(object)))

    if (length(i) > 0) {

        return(paste("missing assay(s)",
            paste(dQuote(names(assay_type)[i]), collapse = ", ")))

    }

    i <- which(!mapply(function(x, t) { is(as.vector(x), t) },
        assays(object)[names(assay_type)], assay_type))

    if (length(i) > 0) {

        return(paste("invalid type for assay(s)",
            paste(dQuote(names(assay_type)[i]), collapse = ", ")))

    }

}

## Class definitions

setClass(
    Class = "TxFeatures",
    slots = c(
        type = "factor",
        txName = "CharacterList",
        geneName = "CharacterList"),
    contains = "GRanges"
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
    contains = "GRanges"
)

setClass(
    Class = "SGSegments",
    slots = c(unlistData = "SGFeatures"),
    contains = "GRangesList"
)

setClass(
    Class = "SGVariants",
    slots = c(unlistData = "SGFeatures"),
    contains = "GRangesList"
)

setClass(
    Class = "SGFeatureCounts",
    slots = c(rowRanges = "SGFeatures"),
    contains = "RangedSummarizedExperiment"
)

setClass(
    Class = "SGVariantCounts",
    slots = c(rowRanges = "SGVariants"),
    contains = "RangedSummarizedExperiment"
)

setClassUnion("Features", c("TxFeatures", "SGFeatures"))
setClassUnion("Paths", c("SGSegments", "SGVariants"))
setClassUnion("Counts", c("SGFeatureCounts", "SGVariantCounts"))

## validity checks

setValidity2("TxFeatures", validTxFeatures)
setValidity2("SGFeatures", validSGFeatures)
setValidity2("SGSegments", validSGSegments)
setValidity2("SGVariants", validSGVariants)
setValidity2("SGFeatureCounts", validSGFeatureCounts)
setValidity2("SGVariantCounts", validSGVariantCounts)

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
##' \code{TxFeatures} extends \code{GRanges} with column slot \code{type}
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
##' @param type Character vector or factor, taking values in
##'   \code{J}, \code{I}, \code{F}, \code{L}, \code{U}
##' @param txName \code{CharacterList} of transcript names or \code{NULL}
##' @param geneName \code{CharacterList} of gene names or \code{NULL}
##' @return \code{TxFeatures} object
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

    slots <- c("type", "txName", "geneName")
    mcols(x) <- mcols(x)[!names(mcols(x)) %in% slots]

    new("TxFeatures", x, type = type, txName = txName, geneName = geneName)

}

##' Creates an instance of S4 class \code{SGFeatures} for storing
##' splice graph features.
##'
##' \code{SGFeatures} extends \code{GRanges} with column slot \code{type}
##' specifying feature type. \code{type} is a factor with levels
##' \code{J} (splice junction), \code{E} (exon bin),
##' \code{D} (splice donor), \code{A} (splice acceptor).
##'
##' \code{splice5p} and \code{splice3p} are logical vectors indicating
##' mandatory splices at the 5' and 3' end of an exon bin, respectively.
##' These are used to determine whether reads extending across the
##' 5' and 3' boundaries of an exon bin must be spliced at the boundary
##' to be considered compatible with the exon bin.
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
##' @param splice5p Logical vector indicating a mandatory splice at the
##'   5' end of an exon bin (determining whether reads extending across
##'   the 5' boundary must be spliced to be considered compatible)
##' @param splice3p Logical vector indicating a mandatory splice at the
##'   3' end of an exon bin (determining whether reads extending across
##'   the 3' boundary must be spliced to be considered compatible)
##' @param featureID Integer vector of feature IDs
##' @param geneID Integer vector of gene IDs
##' @param txName \code{CharacterList} of transcript names or \code{NULL}
##' @param geneName \code{CharacterList} of gene names or \code{NULL}
##' @return \code{SGFeatures} object
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

##' Creates an instance of S4 class \code{SGSegments} for storing
##' splice graph segments.
##'
##' @title Constructor function for S4 class \code{SGSegments}
##' @param x \code{GRangesList} of \code{SGFeatures} with appropriate
##'   outer metadata columns
##' @return \code{SGSegments} object
##' @keywords internal
##' @author Leonard Goldstein

SGSegments <- function(x)
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

    new("SGSegments", x)

}

##' Creates an instance of S4 class \code{SGVariants} for storing
##' splice variants.
##'
##' @title Constructor function for S4 class \code{SGVariants}
##' @param x \code{GRangesList} of \code{SGFeatures} with appropriate
##'   outer metadata columns
##' @return \code{SGVariants} object
##' @examples
##' sgv <- SGVariants()
##' @author Leonard Goldstein

SGVariants <- function(x)
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
            closed5p = logical(),
            closed3p = logical(),
            closed5pEvent = logical(),
            closed3pEvent = logical(),
            featureID5p = IntegerList(),
            featureID3p = IntegerList(),
            featureID5pEvent = IntegerList(),
            featureID3pEvent = IntegerList(),
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

    new("SGVariants", x)

}

##' Creates an instance of S4 class \code{SGFeatureCounts} for storing
##' compatible splice graph feature counts.
##'
##' @title Constructor function for S4 class \code{SGFeatureCounts}
##' @param x \code{RangedSummarizedExperiment} with \code{SGFeatures}
##'   as \code{rowRanges} and assays \dQuote{counts}, \dQuote{FPKM}
##' @return \code{SGFeatureCounts} object
##' @examples
##' sgfc <- SGFeatureCounts()
##' @author Leonard Goldstein

SGFeatureCounts <- function(x)
{

    if (missing(x)) {

        assays <- list(
            counts = matrix(integer(), 0, 0),
            FPKM = matrix(numeric(), 0, 0))

        x <- SummarizedExperiment(assays, rowRanges = SGFeatures())

    }

    new("SGFeatureCounts", x)

}

##' Creates an instance of S4 class \code{SGVariantCounts} for storing
##' splice variant counts.
##'
##' @title Constructor function for S4 class \code{SGFeatureCounts}
##' @param x \code{RangedSummarizedExperiment} with \code{SGVariants}
##'   as \code{rowRanges} and appropriate assays
##' @return \code{SGVariantCounts} object
##' @examples
##' sgvc <- SGVariantCounts()
##' @author Leonard Goldstein

SGVariantCounts <- function(x)
{

    if (missing(x)) {

        assays <- list(
            countsVariant5p = matrix(integer(), 0, 0),
            countsVariant3p = matrix(integer(), 0, 0),
            countsEvent5p = matrix(integer(), 0, 0),
            countsEvent3p = matrix(integer(), 0, 0),
            variantFreq = matrix(numeric(), 0, 0))

        x <- SummarizedExperiment(assays, rowRanges = SGVariants())

    }

    new("SGVariantCounts", x)

}

##' Update object created with previous version of SGSeq.
##'
##' @title Update object
##' @param object Object to be updated
##' @param verbose Should a warning message be generated
##' @return Updated object
##' @author Leonard Goldstein
##' @name updateObject
NULL

##' @rdname updateObject
setMethod("updateObject", "SGVariants", function(object, verbose) {

    current <- names(mcols(object))

    required <- c(
        "from",
        "to",
        "type",
        "featureID",
        "segmentID",
        "geneID",
        "eventID",
        "variantID",
        "closed5p",
        "closed3p",
        "closed5pEvent",
        "closed3pEvent",
        "featureID5p",
        "featureID3p",
        "featureID5pEvent",
        "featureID3pEvent",
        "variantType",
        "variantName",
        "txName",
        "geneName")

    if (all(required %in% current)) return(object)

    if (verbose) {

        msg <- c(
        "SGVariants object was created with previous version of SGSeq.",
        "Please use updateObject() to update. Note that for each event",
        "all variants must be included for updated object to be valid.")
        msg <- paste(msg, collapse = "\n")
        warning(msg, call. = FALSE)

    }

    object <- as(object, "GRangesList")
    object_eventID <- as.factor(mcols(object)$eventID)

    ## add closed5pEvent
    event_closed5p <- tapply(mcols(object)$closed5p, object_eventID, all)
    mcols(object)$closed5pEvent <- as.logical(event_closed5p[match(
        object_eventID, names(event_closed5p))])

    ## add closed3pEvent
    event_closed3p <- tapply(mcols(object)$closed3p, object_eventID, all)
    mcols(object)$closed3pEvent <- as.logical(event_closed3p[match(
        object_eventID, names(event_closed3p))])

    ## add featureID5pEvent
    event_featureID5p <- split(unlist(mcols(object)$featureID5p),
       object_eventID[togroup0(mcols(object)$featureID5p)])
    event_featureID5p <- unique(IntegerList(event_featureID5p))
    mcols(object)$featureID5pEvent <- event_featureID5p[match(
        object_eventID, names(event_featureID5p))]

    ## add featureID3pEvent
    event_featureID3p <- split(unlist(mcols(object)$featureID3p),
       object_eventID[togroup0(mcols(object)$featureID3p)])
    event_featureID3p <- unique(IntegerList(event_featureID3p))
    mcols(object)$featureID3pEvent <- event_featureID3p[match(
        object_eventID, names(event_featureID3p))]

    ## reorder columns
    j <- which(!names(mcols(object)) %in% required)
    if (length(j) > 0) {

        mcols(object) <- cbind(mcols(object)[required], mcols(object)[j])

    } else {

        mcols(object) <- mcols(object)[required]

    }

    ## create new SGVariants object
    object <- SGVariants(object)

    return(object)

})

##' @rdname updateObject
setMethod("updateObject", "SGVariantCounts", function(object, verbose) {

    current <- names(assays(object))

    required <- c(
        "countsVariant5p",
        "countsVariant3p",
        "countsEvent5p",
        "countsEvent3p",
        "variantFreq")

    if (all(required %in% current)) return(object)

    if (verbose) {

        msg <- c(
        "SGVariantCounts object was created with previous version of SGSeq.",
        "Please use updateObject() to update. Note that for each event all",
        "variants must be included for updated object to be valid.")
        msg <- paste(msg, collapse = "\n")
        warning(msg, call. = FALSE)

    }

    object <- as(object, "RangedSummarizedExperiment")

    ## update SGVariants
    rowRanges(object) <- updateObject(rowRanges(object))

    ## remove assay countsTotal
    assays(object)$countsTotal <- NULL

    ## update assay names
    old2new <- c(
        countsTotal5p = "countsEvent5p",
        countsTotal3p = "countsEvent3p",
        countsVariant = "countsVariant5pOr3p")
    j <- which(names(assays(object)) %in% names(old2new))
    names(assays(object))[j] <- old2new[names(assays(object))[j]]

    ## create new SGVariantCounts object
    object <- SGVariantCounts(object)

    return(object)

})
