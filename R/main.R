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
##' path <- system.file("extdata", package = "SGSeq")
##' si$file_bam <- file.path(path, "bams", si$file_bam)
##' sgfc <- analyzeFeatures(si, gr)
##' @author Leonard Goldstein

analyzeFeatures <- function(sample_info, which = NULL,
    features = NULL, predict = is.null(features), 
    alpha = 2, psi = 0.1, beta = 0.2, gamma = 0.2,
    min_n_sample = 1, min_overhang = NA, annotation = NULL,
    max_complexity = 20, verbose = FALSE,
    cores_per_sample = 1, BPPARAM = MulticoreParam(1))
{

    checkSampleInfo(sample_info)
    
    if (is.null(features) && !predict)
        stop("cannot have features NULL and predict FALSE")

    if (!is.null(features) && !is(features, "Features"))
        stop("features must be a TxFeatures or SGFeatures object")

    if (!is.null(annotation) && !is(annotation, "TxFeatures"))
        stop("annotation must be a TxFeatures object")
    
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
            max_complexity = max_complexity,
            verbose = verbose,
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
        verbose = verbose,
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
##' @inheritParams predictTxFeaturesPerSample
##' @param sample_info Data frame with sample information including
##'   mandatory character columns \dQuote{sample_name} and \dQuote{file_bam}.
##' @param yieldSize Number of records used for obtaining alignment
##'   information, or \code{NULL} for all records
##' @param BPPARAM \code{BiocParallelParam} for processing samples in
##'   parallel, defaults to \code{MulticoreParam(1)}
##' @return \code{sample_info} with additional columns \dQuote{paired_end},
##'   \dQuote{read_length}, \dQuote{frag_length}, and \dQuote{lib_size}
##'   if \code{yieldSize} is \code{NULL}
##' @examples 
##' path <- system.file("extdata", package = "SGSeq")
##' si$file_bam <- file.path(path, "bams", si$file_bam)
##' si <- si[, c("sample_name", "file_bam")]
##' si_complete <- getBamInfo(si)
##' @author Leonard Goldstein

getBamInfo <- function(sample_info, yieldSize = NULL, verbose = FALSE,
    BPPARAM = MulticoreParam(1))
{

    checkSampleInfo(sample_info, FALSE)

    list_bamInfo <- bpmapply(
        getBamInfoPerSample,
        file_bam = sample_info$file_bam,
        sample_name = sample_info$sample_name,
        MoreArgs = list(
            yieldSize = yieldSize,
            verbose = verbose),
        SIMPLIFY = FALSE,
        BPPARAM = BPPARAM
    )

    checkApplyResultsForErrors(
        list_bamInfo,
        "getBamInfoPerSample",
        sample_info$sample_name)

    bamInfo <- do.call(rbind, list_bamInfo)

    checkBamInfo(bamInfo)

    cols <- c("paired_end", "read_length", "frag_length", "lib_size")
    cols <- cols[cols %in% names(bamInfo)]
    
    for (col in cols) {
        
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
##' @param sample_info Data frame with sample information.
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
##' path <- system.file("extdata", package = "SGSeq")
##' si$file_bam <- file.path(path, "bams", si$file_bam)
##' txf <- predictTxFeatures(si, gr)
##' @author Leonard Goldstein

predictTxFeatures <- function(sample_info, which = NULL,
    alpha = 2, psi = 0, beta = 0.2, gamma = 0.2,
    min_junction_count = NULL, max_complexity = 20,
    min_n_sample = 1, min_overhang = NA, verbose = FALSE,
    cores_per_sample = 1, BPPARAM = MulticoreParam(1))
{

    checkSampleInfo(sample_info)
    
    list_features <- bpmapply(
        predictTxFeaturesPerSample,
        file_bam = sample_info$file_bam,
        paired_end = sample_info$paired_end,
        read_length = sample_info$read_length,
        frag_length = sample_info$frag_length,
        lib_size = sample_info$lib_size,
        sample_name = sample_info$sample_name,
        MoreArgs = list(
            which = which,
            alpha = alpha,
            psi = psi,
            beta = beta,
            gamma = gamma,
            min_junction_count = min_junction_count,
            include_counts = FALSE,
            retain_coverage = FALSE,
            junctions_only = FALSE,
            max_complexity = max_complexity,
            verbose = verbose,
            cores = cores_per_sample),
        SIMPLIFY = FALSE,
        USE.NAMES = FALSE,
        BPPARAM = BPPARAM
    )

    checkApplyResultsForErrors(
        list_features,
        "predictTxFeaturesPerSample",
        sample_info$sample_name)

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
##' @inheritParams predictTxFeaturesPerSample
##' @param counts_only Logical indicating only counts should be returned
##' @return An \code{SGFeatureCounts} object or integer matrix of counts
##'   if \code{counts_only = TRUE}
##' @examples
##' path <- system.file("extdata", package = "SGSeq")
##' si$file_bam <- file.path(path, "bams", si$file_bam)
##' sgfc <- getSGFeatureCounts(si, sgf)
##' @author Leonard Goldstein

getSGFeatureCounts <- function(sample_info, features, counts_only = FALSE,
    cores_per_sample = 1, verbose = FALSE, BPPARAM = MulticoreParam(1))
{

    checkSampleInfo(sample_info)

    if (!is(features, "SGFeatures"))        
        stop("features must be an SGFeatures object")

    list_counts <- bpmapply(
        getSGFeatureCountsPerSample,
        file_bam = sample_info$file_bam,
        paired_end = sample_info$paired_end,
        sample_name = sample_info$sample_name,
        MoreArgs = list(
            features = features,
            retain_coverage = FALSE,
            verbose = verbose,
            cores = cores_per_sample),
        SIMPLIFY = FALSE,
        USE.NAMES = FALSE,
        BPPARAM = BPPARAM
    )

    checkApplyResultsForErrors(
        list_counts,
        "getSGFeatureCountsPerSample",
        sample_info$sample_name)

    counts <- do.call(cbind, list_counts)

    if (counts_only) return(counts)
    
    sgfc <- makeSGFeatureCounts(
        rowRanges = features,
        colData = sample_info,
        counts = counts)
    
    return(sgfc)
        
}

##' High-level function for the analysis of transcript variants from 
##' splice graph features. Transcript variants are identified with
##' \code{\link{findSGVariants}}. Representative counts and estimated
##' variant frequencies are obtained with \code{\link{getSGVariantCounts}}.
##' 
##' @title Analysis of transcript variants
##' @inheritParams findSGVariants
##' @param object \code{SGFeatureCounts} object
##' @return A \code{SGVariantCounts} object
##' @examples
##' sgvc <- analyzeVariants(sgfc)
##' @author Leonard Goldstein

analyzeVariants <- function(object, maxnvariant = 20, cores = 1)
{

    if (!is(object, "SGFeatureCounts")) 
        stop("object must be an SGFeatureCounts object")

    variants <- findSGVariants(
        features = rowRanges(object),
        maxnvariant = maxnvariant,
        cores = cores)

    counts <- getSGVariantCounts(object, variants, cores)

    return(counts)

}

checkSampleInfo <- function(sample_info, complete = TRUE)
{

    col_type <- c(
         sample_name = "character",
         file_bam = "character",
         paired_end = "logical",
         read_length = "numeric",
         frag_length = "numeric",
         lib_size = "numeric")

    if (!complete) col_type <- col_type[1:2]

    if (!is(sample_info, "data.frame") && !is(sample_info, "DataFrame")) {

        err <- paste("sample_info must be a data frame with columns",
            paste(names(col_type), collapse = ", "))
        stop(err, call. = FALSE)
      
    }

    missing <- !names(col_type) %in% names(sample_info)

    if (any(missing)) {
      
        err <- paste("sample_info is missing column(s)",
            paste(names(col_type)[missing], collapse = ", "))
        stop(err, call. = FALSE)

    }

    invalid <- !mapply(is, sample_info[names(col_type)], col_type)

    if (any(invalid)) {
    
        err <- paste("sample_info column(s)",
            paste(names(col_type)[invalid], collapse = ", "), "\n",
            "must be of type", paste(col_type[invalid], collapse = ", "))
        stop(err, call. = FALSE)
        
    }

}    

checkBamInfo <- function(bam_info)
{

    if (!all(bam_info$XS)) {
      
        files <- bam_info$file_bam[!bam_info$XS]
        files <- paste0("'", files, "'")

        msg <- paste("Custom tag 'XS' not found in BAM file:\n  ", files)
        msg <- paste(msg, collapse = "\n")
        msg <- paste0(msg, "\n\n",
            "Spliced alignments must include a custom tag 'XS'.\n",
            "Compatible BAM files can be obtained with an alignment\n",
            "program for RNA-seq data (e.g. GSNAP, STAR, TopHat).")
        
        warning(msg, call. = FALSE)

    }
  
}
