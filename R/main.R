##' High-level function for the prediction and quantification of
##' splice junctions, exon bins and splice sites from BAM files.
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
##' @inheritParams predictTxFeatures
##' @inheritParams predictTxFeaturesPerSample
##' @inheritParams mergeTxFeatures
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
    max_complexity = 20, verbose = FALSE, cores = 1)
{

    checkSampleInfo(sample_info)
    
    if (is.null(features) && !predict) {
      
        stop("cannot have features NULL and predict FALSE")

    }

    if (!is.null(features) && !is(features, "Features")) {
      
        stop("features must be a TxFeatures or SGFeatures object")

    }

    if (!is.null(annotation) && !is(annotation, "TxFeatures")) {
      
        stop("annotation must be a TxFeatures object")

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
            max_complexity = max_complexity,
            verbose = verbose,
            cores = cores)

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
        cores = cores)
    
    return(counts)
    
}

##' Obtain paired-end status, median aligned read length, 
##' median aligned insert size and library size from BAM files.
##'
##' Library information can be inferred from a subset of BAM records
##' by setting the number of records via argument \code{yieldSize}.
##' Note that library size is only obtained if \code{yieldSize} is {NULL}.
##' 
##' @title Obtain library information from BAM files
##' @param sample_info Data frame with sample information including
##'   mandatory character columns \dQuote{sample_name} and \dQuote{file_bam}.
##' @param yieldSize Number of records used for obtaining library
##'   information, or \code{NULL} for all records
##' @param cores Number of cores available for parallel processing
##' @return \code{sample_info} with additional columns \dQuote{paired_end},
##'   \dQuote{read_length}, \dQuote{frag_length}, and \dQuote{lib_size}
##'   if \code{yieldSize} is \code{NULL}
##' @examples 
##' path <- system.file("extdata", package = "SGSeq")
##' si$file_bam <- file.path(path, "bams", si$file_bam)
##' si <- si[, c("sample_name", "file_bam")]
##' si_complete <- getBamInfo(si)
##' @author Leonard Goldstein

getBamInfo <- function(sample_info, yieldSize = NULL, cores = 1)
{

    checkSampleInfo(sample_info, FALSE)

    list_bamInfo <- mcmapply(
        getBamInfoPerSample,
        file_bam = sample_info$file_bam,
        sample_name = sample_info$sample_name,
        MoreArgs = list(
            yieldSize = yieldSize),
        SIMPLIFY = FALSE,
        mc.preschedule = FALSE,
        mc.cores = cores
    )

    checkApplyResultsForErrors(
        list_bamInfo,
        "getBamInfoPerSample",
        sample_info$sample_name,
        "character")

    bamInfo <- do.call(rbind, list_bamInfo)

    checkBamInfo(bamInfo)

    cols <- c("paired_end", "read_length", "frag_length", "lib_size")
    cols <- cols[cols %in% names(bamInfo)]
    
    for (col in cols) {
        
        sample_info[[col]] <- bamInfo[[col]]
        
    }
    
    return(sample_info)
    
}

##' Splice junctions and exons are predicted for each sample and merged
##' across samples. Terminal exons are filtered and trimmed, if applicable.
##' For details, see the help pages for
##' \code{\link{predictTxFeaturesPerSample}}, \code{\link{mergeTxFeatures}},
##' and \code{\link{processTerminalExons}}.
##' 
##' @title Splice junction and exon prediction from BAM files
##' @inheritParams predictTxFeaturesPerSample
##' @inheritParams mergeTxFeatures
##' @param sample_info Data frame with sample information.
##'   Required columns are \dQuote{sample_name}, \dQuote{file_bam},
##'   \dQuote{paired_end}, \dQuote{read_length}, \dQuote{frag_length}
##'   and \dQuote{lib_size}. Library information can be obtained with
##'   function \code{getBamInfo}.
##' @param min_overhang Minimum overhang required to suppress filtering or
##'   trimming of predicted terminal exons (see the manual page for
##'   \code{processTerminalExons}). Use \code{NULL} to disable processing
##'   (disabling processing is useful if results are subsequently merged
##'   with other predictions and processing is postponed until after the
##'   merging step).
##' @param cores Number of cores available for parallel processing
##' @return A \code{TxFeatures} object
##' @examples
##' path <- system.file("extdata", package = "SGSeq")
##' si$file_bam <- file.path(path, "bams", si$file_bam)
##' txf <- predictTxFeatures(si, gr)
##' @author Leonard Goldstein

predictTxFeatures <- function(sample_info, which = NULL,
    alpha = 2, psi = 0.1, beta = 0.2, gamma = 0.2,
    min_junction_count = NULL, max_complexity = 20,
    min_n_sample = 1, min_overhang = NA, verbose = FALSE,
    cores = 1)
{

    checkSampleInfo(sample_info)

    cores <- setCores(cores, nrow(sample_info))

    list_features <- mcmapply(
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
            cores = cores$per_sample),
        SIMPLIFY = FALSE,
        USE.NAMES = FALSE,
        mc.preschedule = FALSE,
        mc.cores = cores$n_sample
    )

    checkApplyResultsForErrors(
        list_features,
        "predictTxFeaturesPerSample",
        sample_info$sample_name,
        "character")

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
##' sgfc <- getSGFeatureCounts(si, sgf_pred)
##' @author Leonard Goldstein

getSGFeatureCounts <- function(sample_info, features, counts_only = FALSE,
    verbose = FALSE, cores = 1)
{

    checkSampleInfo(sample_info)

    if (!is(features, "SGFeatures")) {
      
        stop("features must be an SGFeatures object")

    }

    cores <- setCores(cores, nrow(sample_info))

    list_counts <- mcmapply(
        getSGFeatureCountsPerSample,
        file_bam = sample_info$file_bam,
        paired_end = sample_info$paired_end,
        sample_name = sample_info$sample_name,
        MoreArgs = list(
            features = features,
            retain_coverage = FALSE,
            verbose = verbose,
            cores = cores$per_sample),
        SIMPLIFY = FALSE,
        USE.NAMES = FALSE,
        mc.preschedule = FALSE,
        mc.cores = cores$n_sample
    )

    checkApplyResultsForErrors(
        list_counts,
        "getSGFeatureCountsPerSample",
        sample_info$sample_name,
        "character")

    counts <- do.call(cbind, list_counts)

    if (counts_only) return(counts)
    
    sgfc <- makeSGFeatureCounts(
        rowRanges = features,
        colData = sample_info,
        counts = counts)
    
    return(sgfc)
        
}

##' High-level function for the analysis of splice variants from 
##' splice graph features. Splice variants are identified with
##' \code{\link{findSGVariants}}. Representative counts are obtained
##' and variant frequencies estimated with \code{\link{getSGVariantCounts}}.
##' 
##' @title Analysis of splice variants
##' @inheritParams findSGVariants
##' @inheritParams getSGVariantCounts
##' @param object \code{SGFeatureCounts} object
##' @return An \code{SGVariantCounts} object
##' @examples
##' sgvc <- analyzeVariants(sgfc_pred)
##' @author Leonard Goldstein

analyzeVariants <- function(object, maxnvariant = 20, min_denominator = NA,
    cores = 1)
{

    if (!is(object, "SGFeatureCounts")) {
      
        stop("object must be an SGFeatureCounts object")

    }

    variants <- findSGVariants(
        features = rowRanges(object),
        maxnvariant = maxnvariant,
        cores = cores)

    counts <- getSGVariantCounts(variants, object,
        min_denominator = min_denominator, cores = cores)

    return(counts)

}

##' For splice variants obtain counts of compatible fragments
##' extending across the start or end of each variant.
##' Counts can be obtained from an \code{SGFeatureCounts} object
##' or from BAM files. Only one of the two arguments \code{object}
##' and \code{sample_info} must be specified. Splice variant
##' frequencies are estimated based on representive counts.
##' 
##' @title Representative counts and frequency estimates for
##'   splice variants
##' @inheritParams predictTxFeatures
##' @inheritParams predictTxFeaturesPerSample
##' @param variants \code{SGVariants} object
##' @param features \code{SGFeatures} object that must include all features
##'   included in featureID5p(variants) and featureID3p(variants)
##' @param object \code{SGFeatureCounts} object
##' @param min_denominator Integer specifying minimum denominator when
##'   calculating variant frequencies. If the denominator is smaller than
##'   \code{min_denominator}, variant frequencies are set to \code{NA}.
##'   If \code{NA}, all variant frequencies are returned.
##' @param cores Number of cores available for parallel processing
##' @return An \code{SGVariantCounts} object
##' @examples
##' sgvc_from_sgfc <- getSGVariantCounts(sgv_pred, sgfc_pred)
##' path <- system.file("extdata", package = "SGSeq")
##' si$file_bam <- file.path(path, "bams", si$file_bam)
##' sgvc_from_bam <- getSGVariantCounts(sgv_pred,
##'   features = sgf_pred, sample_info = si)
##' @author Leonard Goldstein

getSGVariantCounts <- function(variants, object = NULL, features = NULL,
    sample_info = NULL, min_denominator = NA, verbose = FALSE, cores = 1)
{

    if (!is(variants, "SGVariants")) {
      
        stop("'variants' must be an SGVariants object")

    }

    if ((is.null(object) && is.null(sample_info)) ||
        (!is.null(object) && !is.null(sample_info))) {

        stop("Either 'object' or 'sample_info' must not be NULL")

    }

    if (!is.null(sample_info) && is.null(features)) {

        stop("For use with 'sample_info', 'features' must be provided.")

    }    
    
    if (any(table(eventID(variants)) == 1)) {

        msg <- paste(c(
            "Detected events with a single variant.",
            "'variants' must include all variants for a given event to obtain",
            "accurate countsTotal5p, countsTotal3p and variantFreq."),
            collapse = "\n")              
        warning(msg, call. = FALSE)

    }

    if (!is.null(object)) {
    
        sgvc <- getSGVariantCountsFromSGFeatureCounts(
            variants = variants,
            object = object,
            min_denominator = min_denominator, 
            cores = cores)

    } else {

        sgvc <- getSGVariantCountsFromBamFiles(
            variants = variants,
            features = features,
            sample_info = sample_info,
            counts_only = FALSE,
            min_denominator = min_denominator,
            verbose = verbose,
            cores = cores)

    }

    return(sgvc)
    
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

setCores <- function(cores, n_sample)
{
  
    n <- as.integer(max(floor(cores/n_sample), 1))
    s <- as.integer(floor(cores/n))
    list(n_sample = s, per_sample = n, total = s * n)
    
}
