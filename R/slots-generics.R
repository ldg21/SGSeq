##' Accessor and replacement functions for metadata columns.
##'
##' S4 classes defined in the \code{SGSeq} package contain metadata columns
##' that store information for each element in the object. For example, class
##' \code{TxFeatures} contains a column \code{type} that indicates feature
##' type. The specific columns contained in an object depend on its class.
##'
##' To facilitate accessing and modifying metadata columns, for each column
##' there exists a function with name identical to the column name
##' that can be used to access and modify it (see examples).
##'
##' @title Accessing and replacing metadata columns
##' @param object Object containing metadata column
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
##' \code{Counts} objects defined in the \code{SGSeq} package contain
##' different types of assay data. For example, class \code{SGFeatureCounts}
##' contains assays \code{counts} and \code{FPKM}.
##'
##' To facilitate accessing and modifying assays, for each assay
##' there exists a function with name identical to the assay name
##' that can be used to access and modify it (see examples).
##'
##' @title Accessing and replacing assay data
##' @param object Object containing assay data
##' @param value Replacement value
##' @return Assay data for accessor functions, updated object for replacement
##' functions.
##' @examples
##' x <- counts(sgfc_pred)
##' y <- FPKM(sgfc_pred)
##' @author Leonard Goldstein
##' @name assays
NULL

## all classes

##' @rdname slots
setGeneric("type",
    function(object) standardGeneric("type"))

##' @rdname slots
setGeneric("type<-",
    function(object, value) standardGeneric("type<-"))

##' @rdname slots
setGeneric("txName",
    function(object) standardGeneric("txName"))

##' @rdname slots
setGeneric("txName<-",
    function(object, value) standardGeneric("txName<-"))

##' @rdname slots
setGeneric("geneName",
    function(object) standardGeneric("geneName"))

##' @rdname slots
setGeneric("geneName<-",
    function(object, value) standardGeneric("geneName<-"))

## SGFeatures, SGSegments, SGVariants

##' @rdname slots
setGeneric("featureID",
    function(object) standardGeneric("featureID"))

##' @rdname slots
setGeneric("featureID<-",
    function(object, value) standardGeneric("featureID<-"))

##' @rdname slots
setGeneric("geneID",
    function(object) standardGeneric("geneID"))

##' @rdname slots
setGeneric("geneID<-",
    function(object, value) standardGeneric("geneID<-"))

## SGFeatures, SGSegments

##' @rdname slots
setGeneric("splice5p",
    function(object) standardGeneric("splice5p"))

##' @rdname slots
setGeneric("splice5p<-",
    function(object, value) standardGeneric("splice5p<-"))

##' @rdname slots
setGeneric("splice3p",
    function(object) standardGeneric("splice3p"))

##' @rdname slots
setGeneric("splice3p<-",
    function(object, value) standardGeneric("splice3p<-"))

## SGSegments, SGVariants

##' @rdname slots
setGeneric("from",
    function(object) standardGeneric("from"))

##' @rdname slots
setGeneric("from<-",
    function(object, value) standardGeneric("from<-"))

##' @rdname slots
setGeneric("to",
    function(object) standardGeneric("to"))

##' @rdname slots
setGeneric("to<-",
    function(object, value) standardGeneric("to<-"))

##' @rdname slots
setGeneric("segmentID",
    function(object) standardGeneric("segmentID"))

##' @rdname slots
setGeneric("segmentID<-",
    function(object, value) standardGeneric("segmentID<-"))

## SGVariants

##' @rdname slots
setGeneric("variantID",
    function(object) standardGeneric("variantID"))

##' @rdname slots
setGeneric("variantID<-",
    function(object, value) standardGeneric("variantID<-"))

##' @rdname slots
setGeneric("eventID",
    function(object) standardGeneric("eventID"))

##' @rdname slots
setGeneric("eventID<-",
    function(object, value) standardGeneric("eventID<-"))

##' @rdname slots
setGeneric("closed5p",
    function(object) standardGeneric("closed5p"))

##' @rdname slots
setGeneric("closed5p<-",
    function(object, value) standardGeneric("closed5p<-"))

##' @rdname slots
setGeneric("closed3p",
    function(object) standardGeneric("closed3p"))

##' @rdname slots
setGeneric("closed3p<-",
    function(object, value) standardGeneric("closed3p<-"))

##' @rdname slots
setGeneric("closed5pEvent",
    function(object) standardGeneric("closed5pEvent"))

##' @rdname slots
setGeneric("closed5pEvent<-",
    function(object, value) standardGeneric("closed5pEvent<-"))

##' @rdname slots
setGeneric("closed3pEvent",
    function(object) standardGeneric("closed3pEvent"))

##' @rdname slots
setGeneric("closed3pEvent<-",
    function(object, value) standardGeneric("closed3pEvent<-"))

##' @rdname slots
setGeneric("variantType",
    function(object) standardGeneric("variantType"))

##' @rdname slots
setGeneric("variantType<-",
    function(object, value) standardGeneric("variantType<-"))

##' @rdname slots
setGeneric("variantName",
    function(object) standardGeneric("variantName"))

##' @rdname slots
setGeneric("variantName<-",
    function(object, value) standardGeneric("variantName<-"))

##' @rdname slots
setGeneric("featureID5p",
    function(object) standardGeneric("featureID5p"))

##' @rdname slots
setGeneric("featureID5p<-",
    function(object, value) standardGeneric("featureID5p<-"))

##' @rdname slots
setGeneric("featureID3p",
    function(object) standardGeneric("featureID3p"))

##' @rdname slots
setGeneric("featureID3p<-",
    function(object, value) standardGeneric("featureID3p<-"))

##' @rdname slots
setGeneric("featureID5pEvent",
    function(object) standardGeneric("featureID5pEvent"))

##' @rdname slots
setGeneric("featureID5pEvent<-",
    function(object, value) standardGeneric("featureID5pEvent<-"))

##' @rdname slots
setGeneric("featureID3pEvent",
    function(object) standardGeneric("featureID3pEvent"))

##' @rdname slots
setGeneric("featureID3pEvent<-",
    function(object, value) standardGeneric("featureID3pEvent<-"))

## SGFeatureCounts

##' @rdname assays
setGeneric("FPKM",
    function(object) standardGeneric("FPKM"))

##' @rdname assays
setGeneric("FPKM<-",
    function(object, value) standardGeneric("FPKM<-"))

## SGVariantCounts

##' @rdname assays
setGeneric("counts5p",
    function(object) standardGeneric("counts5p"))

##' @rdname assays
setGeneric("counts5p<-",
    function(object, value) standardGeneric("counts5p<-"))

##' @rdname assays
setGeneric("counts3p",
    function(object) standardGeneric("counts3p"))

##' @rdname assays
setGeneric("counts3p<-",
    function(object, value) standardGeneric("counts3p<-"))

##' @rdname assays
setGeneric("counts5pEvent",
    function(object) standardGeneric("counts5pEvent"))

##' @rdname assays
setGeneric("counts5pEvent<-",
    function(object, value) standardGeneric("counts5pEvent<-"))

##' @rdname assays
setGeneric("counts3pEvent",
    function(object) standardGeneric("counts3pEvent"))

##' @rdname assays
setGeneric("counts3pEvent<-",
    function(object, value) standardGeneric("counts3pEvent<-"))

##' @rdname assays
setGeneric("variantFreq",
    function(object) standardGeneric("variantFreq"))

##' @rdname assays
setGeneric("variantFreq<-",
    function(object, value) standardGeneric("variantFreq<-"))
