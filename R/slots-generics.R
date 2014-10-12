##' Accessor and replacement functions for column slots.
##'
##' S4 classes defined in the \code{SGSeq} package contain columns
##' with information for each element in the object. For example, class
##' \code{TxFeatures} contains a column \code{type} that indicates feature
##' type. The specific columns contained in an object depend on its class.
##'
##' To facilitate accessing and modifying columns, for each column
##' there exists a function, with name identical to the column name,
##' that can be used to access and modify it (see examples).
##' 
##' @title Accessing and replacing column slots
##' @param object Object containing column slot
##' @param value Replacement value
##' @return Column value for accessor functions, updated object for
##' replacement functions.
##' @examples
##' head(type(txf))
##' head(type(sgf))
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
##' there exists a function, with name identical to the assay name,
##' that can be used to access and modify it (see examples).
##'
##' @title Accessing and replacing assay data
##' @param object Object containing assay data
##' @param value Replacement value
##' @return Assay data for accessor functions, updated object for replacement
##' functions.
##' @examples
##' x <- counts(sgfc)
##' y <- FPKM(sgfc)
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

## SGFeatures, TxSegments, TxVariants

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

## SGFeatures, TxSegments

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

## TxSegments, TxVariants

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

## TxVariants

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

## SGFeatureCounts

##' @rdname assays
setGeneric("FPKM",
    function(object) standardGeneric("FPKM"))

##' @rdname assays
setGeneric("FPKM<-",
    function(object, value) standardGeneric("FPKM<-"))

## TxVariantCounts

##' @rdname assays
setGeneric("countsVariant5p",
    function(object) standardGeneric("countsVariant5p"))

##' @rdname assays
setGeneric("countsVariant5p<-",
    function(object, value) standardGeneric("countsVariant5p<-"))

##' @rdname assays
setGeneric("countsVariant3p",
    function(object) standardGeneric("countsVariant3p"))

##' @rdname assays
setGeneric("countsVariant3p<-",
    function(object, value) standardGeneric("countsVariant3p<-"))

##' @rdname assays
setGeneric("countsTotal5p",
    function(object) standardGeneric("countsTotal5p"))

##' @rdname assays
setGeneric("countsTotal5p<-",
    function(object, value) standardGeneric("countsTotal5p<-"))

##' @rdname assays
setGeneric("countsTotal3p",
    function(object) standardGeneric("countsTotal3p"))

##' @rdname assays
setGeneric("countsTotal3p<-",
    function(object, value) standardGeneric("countsTotal3p<-"))

##' @rdname assays
setGeneric("variantFreq",
    function(object) standardGeneric("variantFreq"))

##' @rdname assays
setGeneric("variantFreq<-",
    function(object, value) standardGeneric("variantFreq<-"))
