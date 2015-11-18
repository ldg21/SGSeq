## all classes

##' @rdname slots
setMethod("type", "Features",
    function(object) { object@type })

##' @rdname slots
setMethod("type", "Paths",
    function(object) { mcols(object)$type })

##' @rdname slots
setMethod("type", "Counts",
    function(object) { type(rowRanges(object)) })

##' @rdname slots
setReplaceMethod("type", "Features",
    function(object, value) { object@type <- value; object })

##' @rdname slots
setReplaceMethod("type", "Paths",
    function(object, value) { mcols(object)$type <- value; object })

##' @rdname slots
setReplaceMethod("type", "Counts",
    function(object, value) { type(rowRanges(object)) <- value; object })

##' @rdname slots
setMethod("txName", "Features",
    function(object) { object@txName })

##' @rdname slots
setMethod("txName", "Paths",
    function(object) { mcols(object)$txName })

##' @rdname slots
setMethod("txName", "Counts",
    function(object) { txName(rowRanges(object)) })

##' @rdname slots
setReplaceMethod("txName", "Features",
    function(object, value) { object@txName <- value; object })

##' @rdname slots
setReplaceMethod("txName", "Paths",
    function(object, value) { mcols(object)$txName <- value; object })

##' @rdname slots
setReplaceMethod("txName", "Counts",
    function(object, value) { txName(rowRanges(object)) <- value; object })

##' @rdname slots
setMethod("geneName", "Features",
    function(object) { object@geneName })

##' @rdname slots
setMethod("geneName", "Paths",
    function(object) { mcols(object)$geneName })

##' @rdname slots
setMethod("geneName", "Counts",
    function(object) { geneName(rowRanges(object)) })

##' @rdname slots
setReplaceMethod("geneName", "Features",
    function(object, value) { object@geneName <- value; object })

##' @rdname slots
setReplaceMethod("geneName", "Paths",
    function(object, value) { mcols(object)$geneName <- value; object })

##' @rdname slots
setReplaceMethod("geneName", "Counts",
    function(object, value) { geneName(rowRanges(object)) <- value; object })

## SGFeatures, SGSegments, SGVariants

##' @rdname slots
setMethod("featureID", "SGFeatures",
    function(object) { object@featureID })

##' @rdname slots
setMethod("featureID", "Paths",
    function(object) { mcols(object)$featureID })

##' @rdname slots
setMethod("featureID", "Counts",
    function(object) { featureID(rowRanges(object)) })

##' @rdname slots
setReplaceMethod("featureID", "SGFeatures",
    function(object, value) { object@featureID <- value; object })

##' @rdname slots
setReplaceMethod("featureID", "Paths",
    function(object, value) { mcols(object)$featureID <- value; object })

##' @rdname slots
setReplaceMethod("featureID", "Counts",
    function(object, value) { featureID(rowRanges(object)) <- value; object })

##' @rdname slots
setMethod("geneID", "SGFeatures",
    function(object) { object@geneID })

##' @rdname slots
setMethod("geneID", "Paths",
    function(object) { mcols(object)$geneID })

##' @rdname slots
setMethod("geneID", "Counts",
    function(object) { geneID(rowRanges(object)) })

##' @rdname slots
setReplaceMethod("geneID", "SGFeatures",
    function(object, value) { object@geneID <- value; object })

##' @rdname slots
setReplaceMethod("geneID", "Paths",
    function(object, value) { mcols(object)$geneID <- value; object })

##' @rdname slots
setReplaceMethod("geneID", "Counts",
    function(object, value) { geneID(rowRanges(object)) <- value; object })

## SGFeatures, SGSegments

##' @rdname slots
setMethod("splice5p", "SGFeatures",
    function(object) { object@splice5p })

##' @rdname slots
setMethod("splice5p", "SGSegments",
    function(object) { mcols(object)$splice5p })

##' @rdname slots
setMethod("splice5p", "SGFeatureCounts",
    function(object) { splice5p(rowRanges(object)) })

##' @rdname slots
setReplaceMethod("splice5p", "SGFeatures",
    function(object, value) { object@splice5p <- value; object })

##' @rdname slots
setReplaceMethod("splice5p", "SGSegments",
    function(object, value) { mcols(object)$splice5p <- value; object })

##' @rdname slots
setReplaceMethod("splice5p", "SGFeatureCounts",
    function(object, value) { splice5p(rowRanges(object)) <- value; object })

##' @rdname slots
setMethod("splice3p", "SGFeatures",
    function(object) { object@splice3p })

##' @rdname slots
setMethod("splice3p", "SGSegments",
    function(object) { mcols(object)$splice3p })

##' @rdname slots
setMethod("splice3p", "SGFeatureCounts",
    function(object) { splice3p(rowRanges(object)) })

##' @rdname slots
setReplaceMethod("splice3p", "SGFeatures",
    function(object, value) { object@splice3p <- value; object })

##' @rdname slots
setReplaceMethod("splice3p", "SGSegments",
    function(object, value) { mcols(object)$splice3p <- value; object })

##' @rdname slots
setReplaceMethod("splice3p", "SGFeatureCounts",
    function(object, value) { splice3p(rowRanges(object)) <- value; object })

## SGSegments, SGVariants

##' @rdname slots
setMethod("segmentID", "Paths",
    function(object) { mcols(object)$segmentID })

##' @rdname slots
setMethod("segmentID", "SGVariantCounts",
    function(object) { segmentID(rowRanges(object)) })

##' @rdname slots
setReplaceMethod("segmentID", "Paths",
    function(object, value) { mcols(object)$segmentID <- value; object })

##' @rdname slots
setReplaceMethod("segmentID", "SGVariantCounts",
    function(object, value) { segmentID(rowRanges(object)) <- value; object })

##' @rdname slots
setMethod("from", "Paths",
    function(object) { mcols(object)$from })

##' @rdname slots
setMethod("from", "SGVariantCounts",
    function(object) { from(rowRanges(object)) })

##' @rdname slots
setReplaceMethod("from", "Paths",
    function(object, value) { mcols(object)$from <- value; object })

##' @rdname slots
setReplaceMethod("from", "SGVariantCounts",
    function(object, value) { from(rowRanges(object)) <- value; object })

##' @rdname slots
setMethod("to", "Paths",
    function(object) { mcols(object)$to })

##' @rdname slots
setMethod("to", "SGVariantCounts",
    function(object) { to(rowRanges(object)) })

##' @rdname slots
setReplaceMethod("to", "Paths",
    function(object, value) { mcols(object)$to <- value; object })

##' @rdname slots
setReplaceMethod("to", "SGVariantCounts",
    function(object, value) { to(rowRanges(object)) <- value; object })

## SGVariants

##' @rdname slots
setMethod("eventID", "SGVariants",
    function(object) { mcols(object)$eventID })

##' @rdname slots
setMethod("eventID", "SGVariantCounts",
    function(object) { eventID(rowRanges(object)) })

##' @rdname slots
setReplaceMethod("eventID", "SGVariants",
    function(object, value) { mcols(object)$eventID <- value; object })

##' @rdname slots
setReplaceMethod("eventID", "SGVariantCounts",
    function(object, value) { eventID(rowRanges(object)) <- value; object })

##' @rdname slots
setMethod("variantID", "SGVariants",
    function(object) { mcols(object)$variantID })

##' @rdname slots
setMethod("variantID", "SGVariantCounts",
    function(object) { variantID(rowRanges(object)) })

##' @rdname slots
setReplaceMethod("variantID", "SGVariants",
    function(object, value) { mcols(object)$variantID <- value; object })

##' @rdname slots
setReplaceMethod("variantID", "SGVariantCounts",
    function(object, value) { variantID(rowRanges(object)) <- value; object })

##' @rdname slots
setMethod("closed5p", "SGVariants",
    function(object) { mcols(object)$closed5p })

##' @rdname slots
setMethod("closed5p", "SGVariantCounts",
    function(object) { closed5p(rowRanges(object)) })

##' @rdname slots
setReplaceMethod("closed5p", "SGVariants",
    function(object, value) { mcols(object)$closed5p <- value; object })

##' @rdname slots
setReplaceMethod("closed5p", "SGVariantCounts",
    function(object, value) { closed5p(rowRanges(object)) <- value; object })

##' @rdname slots
setMethod("closed3p", "SGVariants",
    function(object) { mcols(object)$closed3p })

##' @rdname slots
setMethod("closed3p", "SGVariantCounts",
    function(object) { closed3p(rowRanges(object)) })

##' @rdname slots
setReplaceMethod("closed3p", "SGVariants",
    function(object, value) { mcols(object)$closed3p <- value; object })

##' @rdname slots
setReplaceMethod("closed3p", "SGVariantCounts",
    function(object, value) { closed3p(rowRanges(object)) <- value; object })

##' @rdname slots
setMethod("closed5pEvent", "SGVariants",
    function(object) { mcols(object)$closed5pEvent })

##' @rdname slots
setMethod("closed5pEvent", "SGVariantCounts",
    function(object) { closed5pEvent(rowRanges(object)) })

##' @rdname slots
setReplaceMethod("closed5pEvent", "SGVariants",
    function(object, value) { mcols(object)$closed5pEvent <- value; object })

##' @rdname slots
setReplaceMethod("closed5pEvent", "SGVariantCounts",
    function(object, value) {
        closed5pEvent(rowRanges(object)) <- value
        object
    })

##' @rdname slots
setMethod("closed3pEvent", "SGVariants",
    function(object) { mcols(object)$closed3pEvent })

##' @rdname slots
setMethod("closed3pEvent", "SGVariantCounts",
    function(object) { closed3pEvent(rowRanges(object)) })

##' @rdname slots
setReplaceMethod("closed3pEvent", "SGVariants",
    function(object, value) { mcols(object)$closed3pEvent <- value; object })

##' @rdname slots
setReplaceMethod("closed3pEvent", "SGVariantCounts",
    function(object, value) {
        closed3pEvent(rowRanges(object)) <- value
        object
    })

##' @rdname slots
setMethod("variantName", "SGVariants",
    function(object) { mcols(object)$variantName })

##' @rdname slots
setMethod("variantName", "SGVariantCounts",
    function(object) { variantName(rowRanges(object)) })

##' @rdname slots
setReplaceMethod("variantName", "SGVariants",
    function(object, value) { mcols(object)$variantName <- value; object })

##' @rdname slots
setReplaceMethod("variantName", "SGVariantCounts",
    function(object, value) {
        variantName(rowRanges(object)) <- value
        object
    })

##' @rdname slots
setMethod("variantType", "SGVariants",
    function(object) { mcols(object)$variantType })

##' @rdname slots
setMethod("variantType", "SGVariantCounts",
    function(object) { variantType(rowRanges(object)) })

##' @rdname slots
setReplaceMethod("variantType", "SGVariants",
    function(object, value) { mcols(object)$variantType <- value; object })

##' @rdname slots
setReplaceMethod("variantType", "SGVariantCounts",
    function(object, value) {
        variantType(rowRanges(object)) <- value
        object
    })

##' @rdname slots
setMethod("featureID5p", "SGVariants",
    function(object) { mcols(object)$featureID5p })

##' @rdname slots
setMethod("featureID5p", "SGVariantCounts",
    function(object) { featureID5p(rowRanges(object)) })

##' @rdname slots
setReplaceMethod("featureID5p", "SGVariants",
    function(object, value) { mcols(object)$featureID5p <- value; object })

##' @rdname slots
setReplaceMethod("featureID5p", "SGVariantCounts",
    function(object, value) {
        featureID5p(rowRanges(object)) <- value
        object
    })

##' @rdname slots
setMethod("featureID3p", "SGVariants",
    function(object) { mcols(object)$featureID3p })

##' @rdname slots
setMethod("featureID3p", "SGVariantCounts",
    function(object) { featureID3p(rowRanges(object)) })

##' @rdname slots
setReplaceMethod("featureID3p", "SGVariants",
    function(object, value) { mcols(object)$featureID3p <- value; object })

##' @rdname slots
setReplaceMethod("featureID3p", "SGVariantCounts",
    function(object, value) {
        featureID3p(rowRanges(object)) <- value
        object
    })

##' @rdname slots
setMethod("featureID5pEvent", "SGVariants",
    function(object) { mcols(object)$featureID5pEvent })

##' @rdname slots
setMethod("featureID5pEvent", "SGVariantCounts",
    function(object) { featureID5pEvent(rowRanges(object)) })

##' @rdname slots
setReplaceMethod("featureID5pEvent", "SGVariants",
    function(object, value) {
        mcols(object)$featureID5pEvent <- value
        object
    })

##' @rdname slots
setReplaceMethod("featureID5pEvent", "SGVariantCounts",
    function(object, value) {
        featureID5pEvent(rowRanges(object)) <- value
        object
    })

##' @rdname slots
setMethod("featureID3pEvent", "SGVariants",
    function(object) { mcols(object)$featureID3pEvent })

##' @rdname slots
setMethod("featureID3pEvent", "SGVariantCounts",
    function(object) { featureID3pEvent(rowRanges(object)) })

##' @rdname slots
setReplaceMethod("featureID3pEvent", "SGVariants",
    function(object, value) {
        mcols(object)$featureID3pEvent <- value
        object
    })

##' @rdname slots
setReplaceMethod("featureID3pEvent", "SGVariantCounts",
    function(object, value) {
        featureID3pEvent(rowRanges(object)) <- value
        object
    })

## SGFeatureCounts

##' @rdname assays
setMethod("counts", "SGFeatureCounts",
    function(object) { assay(object, "counts") })

##' @rdname assays
setReplaceMethod("counts", "SGFeatureCounts",
    function(object, value) {
        assays(object)$counts <- value
        object
    })

##' @rdname assays
setMethod("FPKM", "SGFeatureCounts",
    function(object) { assay(object, "FPKM") })

##' @rdname assays
setReplaceMethod("FPKM", "SGFeatureCounts",
    function(object, value) {
        assays(object)$FPKM <- value
        object
    })

## SGVariantCounts

##' @rdname assays
setMethod("counts5p", "SGVariantCounts",
    function(object) { assay(object, "counts5p") })

##' @rdname assays
setReplaceMethod("counts5p", "SGVariantCounts",
    function(object, value) {
        assays(object)$counts5p <- value
        object
    })

##' @rdname assays
setMethod("counts3p", "SGVariantCounts",
    function(object) { assay(object, "counts3p") })

##' @rdname assays
setReplaceMethod("counts3p", "SGVariantCounts",
    function(object, value) {
        assays(object)$counts3p <- value
        object
    })

##' @rdname assays
setMethod("counts5pEvent", "SGVariantCounts",
    function(object) { assay(object, "counts5pEvent") })

##' @rdname assays
setReplaceMethod("counts5pEvent", "SGVariantCounts",
    function(object, value) {
        assays(object)$counts5pEvent <- value
        object
    })

##' @rdname assays
setMethod("counts3pEvent", "SGVariantCounts",
    function(object) { assay(object, "counts3pEvent") })

##' @rdname assays
setReplaceMethod("counts3pEvent", "SGVariantCounts",
    function(object, value) {
        assays(object)$counts3pEvent <- value
        object
    })

##' @rdname assays
setMethod("variantFreq", "SGVariantCounts",
    function(object) { assay(object, "variantFreq") })

##' @rdname assays
setReplaceMethod("variantFreq", "SGVariantCounts",
    function(object, value) {
        assays(object)$variantFreq <- value
        object
    })

##' @rdname assays
setMethod("counts", "SGVariantCounts",
    function(object) { assay(object, "counts") })

##' @rdname assays
setReplaceMethod("counts", "SGVariantCounts",
    function(object, value) {
        assays(object)$counts <- value
        object
    })
