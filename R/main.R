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
    alpha = 2, psi = 0, beta = 0.2, gamma = 0.2, min_junction_count = NULL,
    min_anchor = 1, min_n_sample = 1, min_overhang = NA, annotation = NULL,
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
            min_junction_count = min_junction_count,
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
        min_anchor = min_anchor,
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
##'   mandatory columns \dQuote{sample_name} and \dQuote{file_bam}.
##'   Column \dQuote{sample_name} must be a character vector. Column
##'   \dQuote{file_bam} can be a character vector or \code{BamFileList}.
##' @param yieldSize Number of records used for obtaining library
##'   information, or \code{NULL} for all records
##' @param cores Number of cores available for parallel processing
##' @return \code{sample_info} with additional columns \dQuote{paired_end},
##'   \dQuote{read_length}, \dQuote{frag_length}, and \dQuote{lib_size}
##'   if \code{yieldSize} is \code{NULL}
##' @examples
##' path <- system.file("extdata", package = "SGSeq")
##' si$file_bam <- file.path(path, "bams", si$file_bam)
##'
##' ## data.frame as sample_info and character vector as file_bam
##' si <- si[, c("sample_name", "file_bam")]
##' si_complete <- getBamInfo(si)
##'
##' ## DataFrame as sample_info and BamFileList as file_bam
##' DF <- DataFrame(si)
##' DF$file_bam <- BamFileList(DF$file_bam)
##' DF_complete <- getBamInfo(DF)
##' @author Leonard Goldstein

getBamInfo <- function(sample_info, yieldSize = NULL, cores = 1)
{

    checkSampleInfo(sample_info, FALSE)

    list_bamInfo <- mcmapply(
        getBamInfoPerSample,
        file_bam = sample_info$file_bam,
        sample_name = sample_info$sample_name,
        MoreArgs = list(yieldSize = yieldSize),
        SIMPLIFY = FALSE,
        mc.preschedule = setPreschedule(cores),
        mc.cores = cores
    )

    checkApplyResultsForErrors(
        list_bamInfo,
        "getBamInfoPerSample",
        sample_info$sample_name,
        "character")

    bamInfo <- do.call(rbindDfsWithoutRowNames, list_bamInfo)

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
##' @return \code{TxFeatures} object
##' @examples
##' path <- system.file("extdata", package = "SGSeq")
##' si$file_bam <- file.path(path, "bams", si$file_bam)
##' txf <- predictTxFeatures(si, gr)
##' @author Leonard Goldstein

predictTxFeatures <- function(sample_info, which = NULL,
    alpha = 2, psi = 0, beta = 0.2, gamma = 0.2,
    min_junction_count = NULL, min_anchor = 1, max_complexity = 20,
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
            min_anchor = min_anchor,
            include_counts = FALSE,
            retain_coverage = FALSE,
            junctions_only = FALSE,
            max_complexity = max_complexity,
            verbose = verbose,
            cores = cores$per_sample),
        SIMPLIFY = FALSE,
        USE.NAMES = FALSE,
        mc.preschedule = setPreschedule(cores$n_sample),
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
##' @return code{SGFeatureCounts} object, or integer matrix of counts
##'   if \code{counts_only = TRUE}
##' @examples
##' path <- system.file("extdata", package = "SGSeq")
##' si$file_bam <- file.path(path, "bams", si$file_bam)
##' sgfc <- getSGFeatureCounts(si, sgf_pred)
##' @author Leonard Goldstein

getSGFeatureCounts <- function(sample_info, features, min_anchor = 1,
    counts_only = FALSE, verbose = FALSE, cores = 1)
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
            min_anchor = min_anchor,
            retain_coverage = FALSE,
            verbose = verbose,
            cores = cores$per_sample),
        SIMPLIFY = FALSE,
        USE.NAMES = FALSE,
        mc.preschedule = setPreschedule(cores$n_sample),
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
        counts = counts,
        min_anchor = min_anchor)

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
##' @inheritParams predictTxFeaturesPerSample
##' @param object \code{SGFeatureCounts} object
##' @return \code{SGVariantCounts} object
##' @examples
##' sgvc <- analyzeVariants(sgfc_pred)
##' @author Leonard Goldstein

analyzeVariants <- function(object, maxnvariant = 20, include = "default",
    min_denominator = NA, min_anchor = 1, cores = 1)
{

    if (!is(object, "SGFeatureCounts")) {

        stop("object must be an SGFeatureCounts object")

    }

    variants <- findSGVariants(
        features = rowRanges(object),
        maxnvariant = maxnvariant,
        include = include,
        cores = cores)

    counts <- getSGVariantCounts(variants, object,
        min_denominator = min_denominator,
        min_anchor = min_anchor, cores = cores)

    return(counts)

}

##' For splice variants, obtain counts of compatible fragments spanning
##' the start and/or end of each variant.
##' Counts can be obtained from an \code{SGFeatureCounts} object
##' or from BAM files. Only one of the two arguments \code{feature_counts}
##' or \code{sample_info} must be specified. Local estimates of relative
##' usage are calculated at the start and/or end of each splice variant.
##' For splice variants with relative usage estimates at both start and end,
##' these are combined by taking a weighted mean, where weights are
##' proportional to the total number of reads spanning the respective
##' boundary.
##'
##' @title Representative counts and frequency estimates for
##'   splice variants
##' @inheritParams predictTxFeatures
##' @inheritParams predictTxFeaturesPerSample
##' @param variants \code{SGVariants} object
##' @param feature_counts \code{SGFeatureCounts} object
##' @param min_denominator Integer specifying minimum denominator when
##'   calculating variant frequencies. The total number of boundary-spanning
##'   reads must be equal to or greater than \code{min_denominator} for at
##'   least one event boundary. Otherwise estimates are set to \code{NA}.
##'   If \code{NA}, all estimates are returned.
##' @param cores Number of cores available for parallel processing
##' @return \code{SGVariantCounts} object
##' @examples
##' sgvc_from_sgfc <- getSGVariantCounts(sgv_pred, sgfc_pred)
##' path <- system.file("extdata", package = "SGSeq")
##' si$file_bam <- file.path(path, "bams", si$file_bam)
##' sgvc_from_bam <- getSGVariantCounts(sgv_pred, sample_info = si)
##' @author Leonard Goldstein

getSGVariantCounts <- function(variants, feature_counts = NULL,
    sample_info = NULL, min_denominator = NA, min_anchor = 1,
    verbose = FALSE, cores = 1)
{

    if (!is(variants, "SGVariants")) {

        stop("'variants' must be an SGVariants object")

    }

    variants <- updateObject(variants, verbose = TRUE)

    if ((is.null(feature_counts) && is.null(sample_info)) ||
        (!is.null(feature_counts) && !is.null(sample_info))) {

        stop("Either 'feature_counts' or 'sample_info' must not be NULL")

    }

    if (!is.null(feature_counts)) {

        f5p <- featureID5pEvent(variants)
        f3p <- featureID3pEvent(variants)
        fid <- unique(c(unlist(f5p), unlist(f3p)))

        if (!all(fid %in% featureID(feature_counts))) {

            stop("feature_counts is missing required features")

        }

        sgvc <- getSGVariantCountsFromSGFeatureCounts(
            variants = variants,
            sgfc = feature_counts,
            min_denominator = min_denominator,
            cores = cores)

    } else {

        sgvc <- getSGVariantCountsFromBamFiles(
            variants = variants,
            sample_info = sample_info,
            counts_only = FALSE,
            min_denominator = min_denominator,
            min_anchor = min_anchor,
            verbose = verbose,
            cores = cores)

    }

    return(sgvc)

}

checkSampleInfo <- function(sample_info, complete = TRUE)
{

    col_type <- list(
         sample_name = "character",
         file_bam = c("character", "BamFileList"),
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

    invalid <- !mapply(isOr, sample_info[names(col_type)], col_type)

    if (any(invalid)) {

        err <- paste("sample_info column(s)",
            paste(names(col_type)[invalid], collapse = ", "), "\n",
            "must be of type",
            paste(unstrsplit(col_type[invalid], "/"), collapse = ", "))
        stop(err, call. = FALSE)

    }

}

checkBamInfo <- function(bam_info)
{

    if (!all(bam_info$XS)) {

        files <- bam_info$file_bam[!bam_info$XS]
        if (is(files, "BamFileList")) files <- path(files)
        files <- paste0("'", files, "'")
        msg <- paste("Custom tag 'XS' not found in BAM file:\n  ", files)
        msg <- paste(msg, collapse = "\n")
        msg <- paste0(msg, "\n\n",
            "Spliced alignments must include a custom tag 'XS'.\n",
            "Compatible BAM files can be obtained with an alignment\n",
            "program for RNA-seq data (e.g. GSNAP, HISAT, STAR, TopHat).")
        warning(msg, call. = FALSE)

    }

}

setCores <- function(cores, n_sample)
{

    n <- as.integer(max(floor(cores/n_sample), 1))
    s <- as.integer(floor(cores/n))
    list(n_sample = s, per_sample = n, total = s * n)

}

setPreschedule <- function(cores)
{

    as.integer(cores) == 1L

}
