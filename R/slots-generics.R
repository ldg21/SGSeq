##' Accessor and replacement functions for metadata columns.
##'
##' S4 classes defined in the \code{SGSeq} package contain metadata columns
##' that store information for each element in the object. For example, class
##' \code{TxFeatures} contains a column \code{type} that indicates feature
##' type. The specific columns contained in an object depend on its class.
##'
##' @title Accessing and replacing metadata columns
##' @param x Object containing metadata column
##' @param value Replacement value
##' @return Content of metadata column for accessor functions or updated
##' object for replacement functions.
##' @examples
##' head(type(txf_ann))
##' head(type(sgf_ann))
##' @author Leonard Goldstein
##' @name slots
NULL

##' Accessor and replacement functions for assay data.
##'
##' Functions \code{counts} and \code{FPKM} are used to extract counts and
##' FPKM values from \code{SGFeatureCounts} and \code{SGVariantCounts}
##' objects. Function \code{variantFreq} is used to access relative usage
##' estimates from \code{SGVariantCounts} objects.
##'
##' @title Accessing and replacing assay data
##' @param object Object containing assay data
##' @param ... Arguments passed to method for \code{SGVariantCounts} objects
##' @param option For \code{SGVariantCounts} objects, \code{option} specifies
##' whether output should be based on counts of reads compatible with the
##' variant at the start (\dQuote{variant5p}), end (\dQuote{variant3p}) or
##' either (\dQuote{variant}), or whether output should be based on counts of
##' reads compatible with any variant belonging to the same event
##' (\dQuote{event5p} or \dQuote{event3p})
##' @param min_anchor For \code{SGVariantCounts} objects, \code{min_anchor}
##' specifies the minimum anchor length when computing FPKM values
##' @param value Replacement value
##' @return Assay data for accessor functions or updated object for
##' replacement functions.
##' @examples
##' x <- counts(sgfc_pred)
##' y <- FPKM(sgfc_pred)
##' u <- counts(sgvc_pred, "variant5p")
##' v <- FPKM(sgvc_pred, "variant5p")
##' @author Leonard Goldstein
##' @name assays
NULL

## all classes

## NOTE generic type() imported from Biostrings

##' @rdname slots
setGeneric("type<-",
    function(x, value) standardGeneric("type<-"))

##' @rdname slots
setGeneric("txName",
    function(x) standardGeneric("txName"))

##' @rdname slots
setGeneric("txName<-",
    function(x, value) standardGeneric("txName<-"))

##' @rdname slots
setGeneric("geneName",
    function(x) standardGeneric("geneName"))

##' @rdname slots
setGeneric("geneName<-",
    function(x, value) standardGeneric("geneName<-"))

## SGFeatures, SGSegments, SGVariants

##' @rdname slots
setGeneric("featureID",
    function(x) standardGeneric("featureID"))

##' @rdname slots
setGeneric("featureID<-",
    function(x, value) standardGeneric("featureID<-"))

##' @rdname slots
setGeneric("geneID",
    function(x) standardGeneric("geneID"))

##' @rdname slots
setGeneric("geneID<-",
    function(x, value) standardGeneric("geneID<-"))

## SGFeatures, SGSegments

##' @rdname slots
setGeneric("splice5p",
    function(x) standardGeneric("splice5p"))

##' @rdname slots
setGeneric("splice5p<-",
    function(x, value) standardGeneric("splice5p<-"))

##' @rdname slots
setGeneric("splice3p",
    function(x) standardGeneric("splice3p"))

##' @rdname slots
setGeneric("splice3p<-",
    function(x, value) standardGeneric("splice3p<-"))

## SGSegments, SGVariants

## NOTE generic from() imported from S4Vectors

##' @rdname slots
setGeneric("from<-",
    function(x, value) standardGeneric("from<-"))

## NOTE generic to() imported from S4Vectors

##' @rdname slots
setGeneric("to<-",
    function(x, value) standardGeneric("to<-"))

##' @rdname slots
setGeneric("segmentID",
    function(x) standardGeneric("segmentID"))

##' @rdname slots
setGeneric("segmentID<-",
    function(x, value) standardGeneric("segmentID<-"))

## SGVariants

##' @rdname slots
setGeneric("variantID",
    function(x) standardGeneric("variantID"))

##' @rdname slots
setGeneric("variantID<-",
    function(x, value) standardGeneric("variantID<-"))

##' @rdname slots
setGeneric("eventID",
    function(x) standardGeneric("eventID"))

##' @rdname slots
setGeneric("eventID<-",
    function(x, value) standardGeneric("eventID<-"))

##' @rdname slots
setGeneric("closed5p",
    function(x) standardGeneric("closed5p"))

##' @rdname slots
setGeneric("closed5p<-",
    function(x, value) standardGeneric("closed5p<-"))

##' @rdname slots
setGeneric("closed3p",
    function(x) standardGeneric("closed3p"))

##' @rdname slots
setGeneric("closed3p<-",
    function(x, value) standardGeneric("closed3p<-"))

##' @rdname slots
setGeneric("closed5pEvent",
    function(x) standardGeneric("closed5pEvent"))

##' @rdname slots
setGeneric("closed5pEvent<-",
    function(x, value) standardGeneric("closed5pEvent<-"))

##' @rdname slots
setGeneric("closed3pEvent",
    function(x) standardGeneric("closed3pEvent"))

##' @rdname slots
setGeneric("closed3pEvent<-",
    function(x, value) standardGeneric("closed3pEvent<-"))

##' @rdname slots
setGeneric("variantType",
    function(x) standardGeneric("variantType"))

##' @rdname slots
setGeneric("variantType<-",
    function(x, value) standardGeneric("variantType<-"))

##' @rdname slots
setGeneric("variantName",
    function(x) standardGeneric("variantName"))

##' @rdname slots
setGeneric("variantName<-",
    function(x, value) standardGeneric("variantName<-"))

##' @rdname slots
setGeneric("featureID5p",
    function(x) standardGeneric("featureID5p"))

##' @rdname slots
setGeneric("featureID5p<-",
    function(x, value) standardGeneric("featureID5p<-"))

##' @rdname slots
setGeneric("featureID3p",
    function(x) standardGeneric("featureID3p"))

##' @rdname slots
setGeneric("featureID3p<-",
    function(x, value) standardGeneric("featureID3p<-"))

##' @rdname slots
setGeneric("featureID5pEvent",
    function(x) standardGeneric("featureID5pEvent"))

##' @rdname slots
setGeneric("featureID5pEvent<-",
    function(x, value) standardGeneric("featureID5pEvent<-"))

##' @rdname slots
setGeneric("featureID3pEvent",
    function(x) standardGeneric("featureID3pEvent"))

##' @rdname slots
setGeneric("featureID3pEvent<-",
    function(x, value) standardGeneric("featureID3pEvent<-"))

## SGFeatureCounts

## NOTE generic counts() imported from BiocGenerics

## NOTE generic 'counts<-'() imported from BiocGenerics

##' @rdname assays
setGeneric("FPKM",
    function(object, ...) standardGeneric("FPKM"))

##' @rdname assays
setGeneric("FPKM<-",
    function(object, ..., value) standardGeneric("FPKM<-"))

## SGVariantCounts

##' @rdname assays
setGeneric("variantFreq",
    function(object) standardGeneric("variantFreq"))

##' @rdname assays
setGeneric("variantFreq<-",
    function(object, value) standardGeneric("variantFreq<-"))
