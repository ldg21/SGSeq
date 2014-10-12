##' Example sample information
##' 
##' Sample information for example BAM files included in the
##' \code{SGSeq} package.
##'
##' @format A \code{data.frame} with columns \dQuote{sample_name},
##' \dQuote{file_bam}, \dQuote{paired_end}, \dQuote{read_length},
##' \dQuote{frag_length} and \dQuote{lib_size}.
##' @name si
NULL

##' Example region of interest
##' 
##' FBXO31 gene locus, based on UCSC knownGene annotation.
##'
##' @format A \code{GRanges} object.
##' @name gr
NULL

##' Example transcript features
##' 
##' Transcript features for FBXO31, based on UCSC knownGene annotation.
##'
##' @format A \code{TxFeatures} object.
##' @name txf
NULL

##' Example splice graph features
##'
##' Splice graph features for FBXO31, predicted from example BAM files.
##'
##' @format An \code{SGFeatures} object.
##' @name sgf
NULL

##' Example splice graph feature counts
##'
##' Compatible counts and FPKMs for predicted FBXO31 splice graph features.
##'
##' @format An \code{SGFeatureCounts} object.
##' @name sgfc
NULL

##' Example transcript variants
##'
##' Transcript variants for FBXO31, based on splice graph features predicted
##' from example BAM files.
##'
##' @format A \code{TxVariants} object.
##' @name txv
NULL

##' Example transcript variant counts
##'
##' Transcript variant counts and frequencies for FBXO31. Splice graph
##' features were predicted and counts obtained from example BAM files.
##'
##' @format A \code{TxVariantCounts} object.
##' @name txvc
NULL
