##' Example sample information
##' 
##' Sample information for example BAM files included in the
##' \code{SGSeq} package.
##'
##' @format A \code{data.frame} with columns \dQuote{sample_name},
##' \dQuote{file_bam}, \dQuote{paired_end}, \dQuote{read_length},
##' \dQuote{frag_length} and \dQuote{lib_size}.
##' @keywords internal
##' @author Leonard Goldstein
##' @name si
NULL

##' Example genomic region of interest
##' 
##' FBXO31 gene locus, based on UCSC knownGene annotation.
##'
##' @format A \code{GRanges} object.
##' @keywords internal
##' @author Leonard Goldstein
##' @name gr
NULL

##' Example transcript features (annotation-based)
##' 
##' Transcript features for FBXO31, based on UCSC knownGene annotation.
##'
##' @format A \code{TxFeatures} object.
##' @keywords internal
##' @author Leonard Goldstein
##' @name txf_ann
NULL

##' Example transcript features (predicted)
##' 
##' Transcript features for FBXO31, predicted from example BAM files.
##'
##' @format A \code{TxFeatures} object.
##' @keywords internal
##' @author Leonard Goldstein
##' @name txf_pred
NULL

##' Example splice graph features (annotation-based)
##'
##' Splice graph features for FBXO31, based on UCSC knownGene annotation.
##'
##' @format An \code{SGFeatures} object.
##' @keywords internal
##' @author Leonard Goldstein
##' @name sgf_ann
NULL

##' Example splice graph features (predicted)
##'
##' Splice graph features for FBXO31, predicted from example BAM files.
##'
##' @format An \code{SGFeatures} object.
##' @keywords internal
##' @author Leonard Goldstein
##' @name sgf_pred
NULL

##' Example splice graph feature counts (annotation-based)
##'
##' Compatible counts and FPKMs for FBXO31 splice graph features,
##' based on UCSC knownGene annotation.
##'
##' @format An \code{SGFeatureCounts} object.
##' @keywords internal
##' @author Leonard Goldstein
##' @name sgfc_ann
NULL

##' Example splice graph feature counts (predicted)
##'
##' Compatible counts and FPKMs for FBXO31 splice graph features,
##' predicted from example BAM files.
##'
##' @format An \code{SGFeatureCounts} object.
##' @keywords internal
##' @author Leonard Goldstein
##' @name sgfc_pred
NULL

##' Example splice variants (annotation-based)
##'
##' Splice variants for FBXO31, based on UCSC knownGene annotation.
##'
##' @format A \code{SGVariants} object.
##' @keywords internal
##' @author Leonard Goldstein
##' @name sgv_ann
NULL

##' Example splice variants (predicted)
##'
##' Splice variants for FBXO31, predicted from example BAM files.
##'
##' @format A \code{SGVariants} object.
##' @keywords internal
##' @author Leonard Goldstein
##' @name sgv_pred
NULL

##' Example splice variant counts (annotated)
##'
##' Splice variant counts and frequencies for FBXO31.
##' Splice variants are based on UCSC knownGene annotation.
##'
##' @format A \code{SGVariantCounts} object.
##' @keywords internal
##' @author Leonard Goldstein
##' @name sgvc_ann
NULL

##' Example splice variant counts (predicted)
##'
##' Splice variant counts and frequencies for FBXO31.
##' Splice variants were predicted from example BAM files.
##'
##' @format A \code{SGVariantCounts} object.
##' @keywords internal
##' @author Leonard Goldstein
##' @name sgvc_pred
NULL

##' Example splice variant counts (annotated) from BAM files
##'
##' Splice variant counts and frequencies for FBXO31.
##' Splice variants are based on UCSC knownGene annotation.
##' Counts were obtained from BAM files.
##'
##' @format A \code{SGVariantCounts} object.
##' @keywords internal
##' @author Leonard Goldstein
##' @name sgvc_ann_from_bam
NULL

##' Example splice variant counts (predicted) from BAM files
##'
##' Splice variant counts and frequencies for FBXO31.
##' Splice variants were predicted from example BAM files.
##' Counts were obtained from BAM files.
##'
##' @format A \code{SGVariantCounts} object.
##' @keywords internal
##' @author Leonard Goldstein
##' @name sgvc_pred_from_bam
NULL
