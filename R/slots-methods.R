## all classes

##' @rdname slots
setMethod("type", "Features",
    function(x) { x@type })

##' @rdname slots
setMethod("type", "Paths",
    function(x) { mcols(x)$type })

##' @rdname slots
setMethod("type", "Counts",
    function(x) { type(rowRanges(x)) })

##' @rdname slots
setReplaceMethod("type", "Features",
    function(x, value) { x@type <- value; x })

##' @rdname slots
setReplaceMethod("type", "Paths",
    function(x, value) { mcols(x)$type <- value; x })

##' @rdname slots
setReplaceMethod("type", "Counts",
    function(x, value) { type(rowRanges(x)) <- value; x })

##' @rdname slots
setMethod("txName", "Features",
    function(x) { x@txName })

##' @rdname slots
setMethod("txName", "Paths",
    function(x) { mcols(x)$txName })

##' @rdname slots
setMethod("txName", "Counts",
    function(x) { txName(rowRanges(x)) })

##' @rdname slots
setReplaceMethod("txName", "Features",
    function(x, value) { x@txName <- value; x })

##' @rdname slots
setReplaceMethod("txName", "Paths",
    function(x, value) { mcols(x)$txName <- value; x })

##' @rdname slots
setReplaceMethod("txName", "Counts",
    function(x, value) { txName(rowRanges(x)) <- value; x })

##' @rdname slots
setMethod("geneName", "Features",
    function(x) { x@geneName })

##' @rdname slots
setMethod("geneName", "Paths",
    function(x) { mcols(x)$geneName })

##' @rdname slots
setMethod("geneName", "Counts",
    function(x) { geneName(rowRanges(x)) })

##' @rdname slots
setReplaceMethod("geneName", "Features",
    function(x, value) { x@geneName <- value; x })

##' @rdname slots
setReplaceMethod("geneName", "Paths",
    function(x, value) { mcols(x)$geneName <- value; x })

##' @rdname slots
setReplaceMethod("geneName", "Counts",
    function(x, value) { geneName(rowRanges(x)) <- value; x })

## SGFeatures, SGSegments, SGVariants

##' @rdname slots
setMethod("featureID", "SGFeatures",
    function(x) { x@featureID })

##' @rdname slots
setMethod("featureID", "Paths",
    function(x) { mcols(x)$featureID })

##' @rdname slots
setMethod("featureID", "Counts",
    function(x) { featureID(rowRanges(x)) })

##' @rdname slots
setReplaceMethod("featureID", "SGFeatures",
    function(x, value) { x@featureID <- value; x })

##' @rdname slots
setReplaceMethod("featureID", "Paths",
    function(x, value) { mcols(x)$featureID <- value; x })

##' @rdname slots
setReplaceMethod("featureID", "Counts",
    function(x, value) { featureID(rowRanges(x)) <- value; x })

##' @rdname slots
setMethod("geneID", "SGFeatures",
    function(x) { x@geneID })

##' @rdname slots
setMethod("geneID", "Paths",
    function(x) { mcols(x)$geneID })

##' @rdname slots
setMethod("geneID", "Counts",
    function(x) { geneID(rowRanges(x)) })

##' @rdname slots
setReplaceMethod("geneID", "SGFeatures",
    function(x, value) { x@geneID <- value; x })

##' @rdname slots
setReplaceMethod("geneID", "Paths",
    function(x, value) { mcols(x)$geneID <- value; x })

##' @rdname slots
setReplaceMethod("geneID", "Counts",
    function(x, value) { geneID(rowRanges(x)) <- value; x })

## SGFeatures, SGSegments

##' @rdname slots
setMethod("splice5p", "SGFeatures",
    function(x) { x@splice5p })

##' @rdname slots
setMethod("splice5p", "SGSegments",
    function(x) { mcols(x)$splice5p })

##' @rdname slots
setMethod("splice5p", "SGFeatureCounts",
    function(x) { splice5p(rowRanges(x)) })

##' @rdname slots
setReplaceMethod("splice5p", "SGFeatures",
    function(x, value) { x@splice5p <- value; x })

##' @rdname slots
setReplaceMethod("splice5p", "SGSegments",
    function(x, value) { mcols(x)$splice5p <- value; x })

##' @rdname slots
setReplaceMethod("splice5p", "SGFeatureCounts",
    function(x, value) { splice5p(rowRanges(x)) <- value; x })

##' @rdname slots
setMethod("splice3p", "SGFeatures",
    function(x) { x@splice3p })

##' @rdname slots
setMethod("splice3p", "SGSegments",
    function(x) { mcols(x)$splice3p })

##' @rdname slots
setMethod("splice3p", "SGFeatureCounts",
    function(x) { splice3p(rowRanges(x)) })

##' @rdname slots
setReplaceMethod("splice3p", "SGFeatures",
    function(x, value) { x@splice3p <- value; x })

##' @rdname slots
setReplaceMethod("splice3p", "SGSegments",
    function(x, value) { mcols(x)$splice3p <- value; x })

##' @rdname slots
setReplaceMethod("splice3p", "SGFeatureCounts",
    function(x, value) { splice3p(rowRanges(x)) <- value; x })

## SGSegments, SGVariants

##' @rdname slots
setMethod("segmentID", "Paths",
    function(x) { mcols(x)$segmentID })

##' @rdname slots
setMethod("segmentID", "SGVariantCounts",
    function(x) { segmentID(rowRanges(x)) })

##' @rdname slots
setReplaceMethod("segmentID", "Paths",
    function(x, value) { mcols(x)$segmentID <- value; x })

##' @rdname slots
setReplaceMethod("segmentID", "SGVariantCounts",
    function(x, value) { segmentID(rowRanges(x)) <- value; x })

##' @rdname slots
setMethod("from", "Paths",
    function(x) { mcols(x)$from })

##' @rdname slots
setMethod("from", "SGVariantCounts",
    function(x) { from(rowRanges(x)) })

##' @rdname slots
setReplaceMethod("from", "Paths",
    function(x, value) { mcols(x)$from <- value; x })

##' @rdname slots
setReplaceMethod("from", "SGVariantCounts",
    function(x, value) { from(rowRanges(x)) <- value; x })

##' @rdname slots
setMethod("to", "Paths",
    function(x) { mcols(x)$to })

##' @rdname slots
setMethod("to", "SGVariantCounts",
    function(x) { to(rowRanges(x)) })

##' @rdname slots
setReplaceMethod("to", "Paths",
    function(x, value) { mcols(x)$to <- value; x })

##' @rdname slots
setReplaceMethod("to", "SGVariantCounts",
    function(x, value) { to(rowRanges(x)) <- value; x })

## SGVariants

##' @rdname slots
setMethod("eventID", "SGVariants",
    function(x) { mcols(x)$eventID })

##' @rdname slots
setMethod("eventID", "SGVariantCounts",
    function(x) { eventID(rowRanges(x)) })

##' @rdname slots
setReplaceMethod("eventID", "SGVariants",
    function(x, value) { mcols(x)$eventID <- value; x })

##' @rdname slots
setReplaceMethod("eventID", "SGVariantCounts",
    function(x, value) { eventID(rowRanges(x)) <- value; x })

##' @rdname slots
setMethod("variantID", "SGVariants",
    function(x) { mcols(x)$variantID })

##' @rdname slots
setMethod("variantID", "SGVariantCounts",
    function(x) { variantID(rowRanges(x)) })

##' @rdname slots
setReplaceMethod("variantID", "SGVariants",
    function(x, value) { mcols(x)$variantID <- value; x })

##' @rdname slots
setReplaceMethod("variantID", "SGVariantCounts",
    function(x, value) { variantID(rowRanges(x)) <- value; x })

##' @rdname slots
setMethod("closed5p", "SGVariants",
    function(x) { mcols(x)$closed5p })

##' @rdname slots
setMethod("closed5p", "SGVariantCounts",
    function(x) { closed5p(rowRanges(x)) })

##' @rdname slots
setReplaceMethod("closed5p", "SGVariants",
    function(x, value) { mcols(x)$closed5p <- value; x })

##' @rdname slots
setReplaceMethod("closed5p", "SGVariantCounts",
    function(x, value) { closed5p(rowRanges(x)) <- value; x })

##' @rdname slots
setMethod("closed3p", "SGVariants",
    function(x) { mcols(x)$closed3p })

##' @rdname slots
setMethod("closed3p", "SGVariantCounts",
    function(x) { closed3p(rowRanges(x)) })

##' @rdname slots
setReplaceMethod("closed3p", "SGVariants",
    function(x, value) { mcols(x)$closed3p <- value; x })

##' @rdname slots
setReplaceMethod("closed3p", "SGVariantCounts",
    function(x, value) { closed3p(rowRanges(x)) <- value; x })

##' @rdname slots
setMethod("closed5pEvent", "SGVariants",
    function(x) { mcols(x)$closed5pEvent })

##' @rdname slots
setMethod("closed5pEvent", "SGVariantCounts",
    function(x) { closed5pEvent(rowRanges(x)) })

##' @rdname slots
setReplaceMethod("closed5pEvent", "SGVariants",
    function(x, value) { mcols(x)$closed5pEvent <- value; x })

##' @rdname slots
setReplaceMethod("closed5pEvent", "SGVariantCounts",
    function(x, value) { closed5pEvent(rowRanges(x)) <- value; x })

##' @rdname slots
setMethod("closed3pEvent", "SGVariants",
    function(x) { mcols(x)$closed3pEvent })

##' @rdname slots
setMethod("closed3pEvent", "SGVariantCounts",
    function(x) { closed3pEvent(rowRanges(x)) })

##' @rdname slots
setReplaceMethod("closed3pEvent", "SGVariants",
    function(x, value) { mcols(x)$closed3pEvent <- value; x })

##' @rdname slots
setReplaceMethod("closed3pEvent", "SGVariantCounts",
    function(x, value) { closed3pEvent(rowRanges(x)) <- value; x })

##' @rdname slots
setMethod("variantName", "SGVariants",
    function(x) { mcols(x)$variantName })

##' @rdname slots
setMethod("variantName", "SGVariantCounts",
    function(x) { variantName(rowRanges(x)) })

##' @rdname slots
setReplaceMethod("variantName", "SGVariants",
    function(x, value) { mcols(x)$variantName <- value; x })

##' @rdname slots
setReplaceMethod("variantName", "SGVariantCounts",
    function(x, value) { variantName(rowRanges(x)) <- value; x })

##' @rdname slots
setMethod("variantType", "SGVariants",
    function(x) { mcols(x)$variantType })

##' @rdname slots
setMethod("variantType", "SGVariantCounts",
    function(x) { variantType(rowRanges(x)) })

##' @rdname slots
setReplaceMethod("variantType", "SGVariants",
    function(x, value) { mcols(x)$variantType <- value; x })

##' @rdname slots
setReplaceMethod("variantType", "SGVariantCounts",
    function(x, value) { variantType(rowRanges(x)) <- value; x })

##' @rdname slots
setMethod("featureID5p", "SGVariants",
    function(x) { mcols(x)$featureID5p })

##' @rdname slots
setMethod("featureID5p", "SGVariantCounts",
    function(x) { featureID5p(rowRanges(x)) })

##' @rdname slots
setReplaceMethod("featureID5p", "SGVariants",
    function(x, value) { mcols(x)$featureID5p <- value; x })

##' @rdname slots
setReplaceMethod("featureID5p", "SGVariantCounts",
    function(x, value) { featureID5p(rowRanges(x)) <- value; x })

##' @rdname slots
setMethod("featureID3p", "SGVariants",
    function(x) { mcols(x)$featureID3p })

##' @rdname slots
setMethod("featureID3p", "SGVariantCounts",
    function(x) { featureID3p(rowRanges(x)) })

##' @rdname slots
setReplaceMethod("featureID3p", "SGVariants",
    function(x, value) { mcols(x)$featureID3p <- value; x })

##' @rdname slots
setReplaceMethod("featureID3p", "SGVariantCounts",
    function(x, value) { featureID3p(rowRanges(x)) <- value; x })

##' @rdname slots
setMethod("featureID5pEvent", "SGVariants",
    function(x) { mcols(x)$featureID5pEvent })

##' @rdname slots
setMethod("featureID5pEvent", "SGVariantCounts",
    function(x) { featureID5pEvent(rowRanges(x)) })

##' @rdname slots
setReplaceMethod("featureID5pEvent", "SGVariants",
    function(x, value) { mcols(x)$featureID5pEvent <- value; x })

##' @rdname slots
setReplaceMethod("featureID5pEvent", "SGVariantCounts",
    function(x, value) { featureID5pEvent(rowRanges(x)) <- value; x })

##' @rdname slots
setMethod("featureID3pEvent", "SGVariants",
    function(x) { mcols(x)$featureID3pEvent })

##' @rdname slots
setMethod("featureID3pEvent", "SGVariantCounts",
    function(x) { featureID3pEvent(rowRanges(x)) })

##' @rdname slots
setReplaceMethod("featureID3pEvent", "SGVariants",
    function(x, value) { mcols(x)$featureID3pEvent <- value; x })

##' @rdname slots
setReplaceMethod("featureID3pEvent", "SGVariantCounts",
    function(x, value) { featureID3pEvent(rowRanges(x)) <- value; x })

## SGFeatureCounts

##' @rdname assays
setMethod("counts", "SGFeatureCounts",
    function(object) { assay(object, "counts") })

##' @rdname assays
setReplaceMethod("counts", "SGFeatureCounts",
    function(object, value) { assays(object)$counts <- value; object })

##' @rdname assays
setMethod("FPKM", "SGFeatureCounts",
    function(object) { assay(object, "FPKM") })

##' @rdname assays
setReplaceMethod("FPKM", "SGFeatureCounts",
    function(object, value) { assays(object)$FPKM <- value; object })

## SGVariantCounts

##' @rdname assays
setMethod("counts", "SGVariantCounts",
    function(object, option = "variant5pOr3p") {

        valid <- c("variant5p", "variant3p", "event5p", "event3p",
            "variant5pOr3p")

        if (missing(option) || !option %in% valid) {

            msg <- paste("argument option must be one of\n",
                paste(paste0("\"", valid, "\""), collapse = ", "))
            stop(msg)

        }

        if (option == "variant5pOr3p" &&
            !"countsVariant5pOr3p" %in% names(assays(object))) {

            msg <- paste("object is missing assay \"countsVariant5pOr3p\",",
                "run getSGVariantCounts() first")
            stop(msg)

        }

        assay_name <- switch(option,
            variant5p = "countsVariant5p",
            variant3p = "countsVariant3p",
            event5p = "countsEvent5p",
            event3p = "countsEvent3p",
            variant5pOr3p = "countsVariant5pOr3p")


        assay(object, assay_name)

})

##' @rdname assays
setReplaceMethod("counts", "SGVariantCounts",
    function(object, option, value) {

        valid <- c("variant5p", "variant3p", "event5p", "event3p",
            "variant5pOr3p")

        if (missing(option) || !option %in% valid) {

            msg <- paste("argument option must be one of\n",
                paste(paste0("\"", valid, "\""), collapse = ", "))
            stop(msg)

        }

        if (option == "variant5pOr3p" &&
            !"countsVariant5pOr3p" %in% names(assays(object))) {

            msg <- paste("object is missing assay \"countsVariant5pOr3p\",",
                "run getSGVariantCounts() first")
            stop(msg)

        }

        assay_name <- switch(option,
            variant5p = "countsVariant5p",
            variant3p = "countsVariant3p",
            event5p = "countsEvent5p",
            event3p = "countsEvent3p",
            variant5pOr3p = "countsVariant5pOr3p")

        assays(object)[[assay_name]] <- value

        object

})

##' @rdname assays
setMethod("FPKM", "SGVariantCounts",
    function(object, option, min_anchor = 1) {

        valid <- c("variant5p", "variant3p", "event5p", "event3p")

        if (missing(option) || !option %in% valid) {

            msg <- paste("argument option must be one of\n",
                paste(paste0("\"", valid, "\""), collapse = ", "))
            stop(msg)

        }

        getScaledCounts(object = object, min_anchor = min_anchor,
            counts = counts(object, option))

})

##' @rdname assays
setMethod("variantFreq", "SGVariantCounts",
    function(object) { assay(object, "variantFreq") })

##' @rdname assays
setReplaceMethod("variantFreq", "SGVariantCounts",
    function(object, value) { assays(object)$variantFreq <- value; object })
