##' High-level function for the prediction and quantification of
##' splice junctions, exon bins and splice sites from BAM files.
##'
##' If alignment information is not included in \code{sample_info},
##' it is obtained directly from BAM files with \code{\link{getBamInfo}}.
##'
##' Splice junctions and exons are predicted from BAM files with
##' \code{\link{predictTxFeatures}}.
##'
##' Known features can be provided as \code{TxFeatures} or
##' \code{SGFeatures} via argument \code{features}. 
##'
##' If \code{features} is not \code{NULL} and \code{predict} is
##' \code{TRUE}, known features are augmented with predictions.
##'
##' Known and/or predicted transcript features are converted to splice
##' graph features. For details, see \code{\link{convertToSGFeatures}}.
##' 
##' Optionally, splice graph features can be annotated with respect to
##' a \code{TxFeatures} object provided via argument \code{annotation}.
##' For details, see the help page for function \code{\link{annotate}}.
##'
##' Finally, compatible fragment counts for splice graph features are
##' obtained from BAM files with \code{\link{getSGFeatureCounts}}.
##' 
##' @title Analysis of splice graph features from BAM files
##' @inheritParams getBamInfo
##' @inheritParams predictTxFeatures
##' @inheritParams predictTxFeaturesPerSample
##' @inheritParams mergeTxFeatures
##' @inheritParams processTerminalExons
##' @param features \code{TxFeatures} or \code{SGFeatures} object
##' @param predict Logical indicating whether transcript
##'   features should be predicted from BAM files
##' @param alpha Minimum FPKM required for a splice junction to be included
##' @param annotation \code{TxFeatures} object used for annotation
##' @return \code{SGFeatureCounts} object
##' @examples
##' dir <- system.file("extdata", package = "SGSeq")
##' si$file_bam <- file.path(dir, "bams", si$file_bam)
##' sgfc <- analyzeFeatures(si, gr)
##' @author Leonard Goldstein

analyzeFeatures <- function(sample_info, which = NULL,
    features = NULL, predict = is.null(features), 
    alpha = 2, psi = 0.1, beta = 0.2, gamma = 0.2,
    min_n_sample = 1, min_overhang = NA, annotation = NULL, 
    cores_per_sample = 1, BPPARAM = MulticoreParam(1))
{
    
    if (!validSampleInfo(sample_info))
        stop("sample_info must be a data.frame including
            character columns sample_name, file_bam")
    
    if (is.null(features) && !predict)
        stop("cannot have features NULL and predict FALSE")

    if (!is.null(features) && !is(features, "Features"))
        stop("features must be a TxFeatures or SGFeatures object")

    if (!is.null(annotation) && !is(annotation, "TxFeatures"))
        stop("annotation must be a TxFeatures object")
    
    if (!validBamInfo(sample_info)) {

        message("Obtain BAM info...")        
        sample_info <- getBamInfo(
            sample_info = sample_info,
            BPPARAM = BPPARAM)

    }

    if (predict) {

        message("Predict features...")
        predicted <- predictTxFeatures(
            sample_info = sample_info,
            which = which,
            alpha = alpha,
            psi = psi,
            beta = beta,
            gamma = gamma,
            min_n_sample = min_n_sample,
            min_overhang = min_overhang,
            cores_per_sample = cores_per_sample,
            BPPARAM = BPPARAM)

        if (!is.null(features)) {

            if (is(features, "TxFeatures")) {

                message("Merge features...")
                features <- mergeTxFeatures(predicted, features)
                
            } else {
            
                message("Process features...")
                predicted <- convertToSGFeatures(predicted)
                
                message("Merge features...")
                features <- mergeSGFeatures(predicted, features)

            }
            
        } else {

            features <- predicted

        }

    }

    if (is(features, "TxFeatures")) {

        message("Process features...")
        features <- convertToSGFeatures(features)
        
    }
        
    if (!is.null(annotation)) {

        message("Annotate features...")
        features <- annotate(features, annotation)

    }
    
    message("Obtain counts...")
    counts <- getSGFeatureCounts(
        sample_info = sample_info,
        features = features,
        cores_per_sample = cores_per_sample,
        BPPARAM = BPPARAM)
    
    return(counts)
    
}

##' Obtain paired-end status, median aligned read length, 
##' median aligned insert size and library size from BAM file.
##'
##' Alignment information can be inferred from a subset of BAM records
##' by setting the number of records via argument \code{yieldSize}.
##' Note that library size can only be obtained if \code{yieldSize} is {NULL}.
##' 
##' @title Obtain alignment information from BAM files
##' @param sample_info \code{data.frame} with sample information including
##'   mandatory character columns \dQuote{sample_name} and \dQuote{file_bam}.
##' @param yieldSize Number of records used for obtaining alignment
##'   information, or \code{NULL} for all records
##' @param BPPARAM \code{BiocParallelParam} for processing samples in
##'   parallel, defaults to \code{MulticoreParam(1)}
##' @return \code{sample_info} with additional columns \dQuote{paired_end},
##'   \dQuote{read_length}, \dQuote{frag_length}, and \dQuote{lib_size}
##'   if \code{yieldSize} is \code{NULL}
##' @examples 
##' dir <- system.file("extdata", package = "SGSeq")
##' si$file_bam <- file.path(dir, "bams", si$file_bam)
##' si <- si[, c("sample_name", "file_bam")]
##' si_complete <- getBamInfo(si)
##' @author Leonard Goldstein

getBamInfo <- function(sample_info, yieldSize = NULL,
    BPPARAM = MulticoreParam(1))
{

    if (!validSampleInfo(sample_info))
        stop("sample_info must be a data.frame including
            character columns sample_name, file_bam")

    list_bamInfo <- bplapply(
        sample_info$file_bam,
        getBamInfoPerSample,
        yieldSize = yieldSize,
        BPPARAM = BPPARAM
    )
    
    bamInfo <- do.call(rbind, list_bamInfo)

    for (col in names(bamInfo)) {
        
        sample_info[[col]] <- bamInfo[[col]]
        
    }
    
    return(sample_info)
    
}

##' Transcript features are predicted for each sample and merged across
##' samples. Subsequently, terminal exons are filtered and trimmed
##' (if applicable). For details, see the help pages for
##' \code{\link{predictTxFeaturesPerSample}}, \code{\link{mergeTxFeatures}},
##' and \code{\link{processTerminalExons}}.
##' 
##' @title Transcript feature prediction from BAM files
##' @inheritParams predictTxFeaturesPerSample
##' @inheritParams mergeTxFeatures
##' @param sample_info \code{data.frame} with sample information.
##'   Required columns are \dQuote{sample_name}, \dQuote{file_bam},
##'   \dQuote{paired_end}, \dQuote{read_length}, \dQuote{frag_length}
##'   and \dQuote{lib_size}. Alignment information can be obtained with
##'   function \code{getBamInfo}.
##' @param min_overhang After merging, terminal exons are processed.
##'   For terminal exons sharing a splice site with an internal exon,
##'   minimum overhang required for terminal exons to be included.
##'   For remaining terminal exons overlapping other exons, minimum 
##'   overhang required to suppress trimming. Use \code{NA} to remove all
##'   terminal exons sharing a splice with an internal exon and trim all
##'   remaining terminal exons overlapping other exons. Use \code{NULL}
##'   to disable processing (disabling processing is useful if results are
##'   subsequently merged with other predictions and processing is
##'   postponed until after the merging step).
##' @param cores_per_sample Number of cores per sample
##' @param BPPARAM \code{BiocParallelParam} for processing samples in parallel,
##'   defaults to \code{MulticoreParam(1)}
##' @return A \code{TxFeatures} object
##' @examples
##' dir <- system.file("extdata", package = "SGSeq")
##' si$file_bam <- file.path(dir, "bams", si$file_bam)
##' txf <- predictTxFeatures(si, gr)
##' @author Leonard Goldstein

predictTxFeatures <- function(sample_info, which = NULL,
    alpha = 2, psi = 0, beta = 0.2, gamma = 0.2,
    min_junction_count = NULL, min_n_sample = 1, min_overhang = NA,
    cores_per_sample = 1, BPPARAM = MulticoreParam(1))
{

    if (!validSampleInfo(sample_info))
        stop("sample_info must be a data.frame including
            character columns sample_name, file_bam")

    if (!validBamInfo(sample_info))
        stop("Incomplete sample_info")
    
    list_features <- bpmapply(
        predictTxFeaturesPerSample,
        file_bam = sample_info$file_bam,
        paired_end = sample_info$paired_end,
        read_length = sample_info$read_length,
        frag_length = sample_info$frag_length,
        lib_size = sample_info$lib_size,
        MoreArgs = list(
            which = which,
            alpha = alpha,
            psi = psi,
            beta = beta,
            gamma = gamma,
            min_junction_count = min_junction_count,
            include_counts = FALSE,
            retain_coverage = FALSE,
            cores = cores_per_sample),
        USE.NAMES = FALSE,
        BPPARAM = BPPARAM
    )
    
    features <- mergeTxFeatures(list_features, min_n_sample = min_n_sample)

    if (!is.null(min_overhang)) {

        features <- processTerminalExons(features, min_overhang)

    }

    return(features)
    
}

##' Compatible counts are obtained for each sample and combined into
##' an \code{SGFeatureCounts} object. 
##'
##' @title Compatible counts for splice graph features from BAM files
##' @inheritParams getSGFeatureCountsPerSample
##' @inheritParams predictTxFeatures
##' @return An \code{SGFeatureCounts} object
##' @examples
##' dir <- system.file("extdata", package = "SGSeq")
##' si$file_bam <- file.path(dir, "bams", si$file_bam)
##' sgfc <- getSGFeatureCounts(si, sgf)
##' @author Leonard Goldstein

getSGFeatureCounts <- function(sample_info, features, cores_per_sample = 1,
    BPPARAM = MulticoreParam(1))
{

    if (!validSampleInfo(sample_info))
        stop("sample_info must be a data.frame including
            character columns sample_name, file_bam")

    if (!validBamInfo(sample_info))
        stop("Incomplete sample_info")
    
    if (!is(features, "SGFeatures"))        
        stop("features must be an SGFeatures object")
    
    list_N <- bpmapply(
        getSGFeatureCountsPerSample,
        file_bam = sample_info$file_bam,
        paired_end = sample_info$paired_end,
        MoreArgs = list(
            features = features,
            cores = cores_per_sample),
        SIMPLIFY = FALSE,
        USE.NAMES = FALSE,
        BPPARAM = BPPARAM
    )

    N <- do.call(cbind, list_N)    

    counts <- makeSGFeatureCounts(
        rowData = features,
        colData = sample_info,
        counts = N)
    
    return(counts)
        
}

##' High-level function for the analysis of transcript variants from 
##' splice graph features. Transcript variants are identified with
##' \code{\link{findTxVariants}}. Representative counts and estimated
##' variant frequencies are obtained with \code{\link{getTxVariantCounts}}.
##' 
##' @title Analysis of transcript variants
##' @inheritParams findTxVariants
##' @param object \code{SGFeatureCounts} object
##' @return A \code{TxVariantCounts} object
##' @examples
##' txvc <- analyzeVariants(sgfc)
##' @author Leonard Goldstein

analyzeVariants <- function(object, maxnvariant = 20, cores = 1)
{

    if (!is(object, "SGFeatureCounts")) 
        stop("object must be an SGFeatureCounts object")

    variants <- findTxVariants(
        features = rowData(object),
        maxnvariant = maxnvariant,
        cores = cores)

    counts <- getTxVariantCounts(object, variants)

    return(counts)

}

validSampleInfo <- function(object)
{

    if (!is.data.frame(object)) {
      
        return(FALSE)

    }
      
    col_type <- c("sample_name" = "character", "file_bam" = "character")

    if (!all(names(col_type) %in% names(object))) {

        return(FALSE)

    }

    if (!all(mapply(is, object[names(col_type)], col_type))) {

        return(FALSE)

    }

    return(TRUE)

}
    
validBamInfo <- function(object)
{

    required_cols <- c("paired_end", "read_length", "frag_length", "lib_size")
    
    if (!all(required_cols %in% names(object))) {

        return(FALSE)

    }

    return(TRUE)

}
