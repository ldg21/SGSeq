##' Obtain counts of compatible fragments for splice graph features.
##'
##' @title Compatible fragment counts for splice graph features
##' @inheritParams predictTxFeaturesPerSample
##' @param features \code{SGFeatures} object
##' @return Numeric vector of compatible fragment counts
##' @keywords internal
##' @author Leonard Goldstein

getSGFeatureCountsPerSample <- function(features, file_bam, paired_end,
    sample_name, min_anchor, retain_coverage, verbose, cores)
{

    list_range <- range(features)
    hits <- findOverlaps(features, list_range)
    list_index <- split(queryHits(hits), subjectHits(hits))
    list_features <- split(features[queryHits(hits)], subjectHits(hits))

    list_counts <- mclapply(
        list_features,
        getSGFeatureCountsPerStrand,
        file_bam = file_bam,
        paired_end = paired_end,
        sample_name = sample_name,
        min_anchor = min_anchor,
        retain_coverage = retain_coverage,
        verbose = verbose,
        mc.preschedule = setPreschedule(cores),
        mc.cores = cores)

    checkApplyResultsForErrors(
        list_counts,
        "getSGFeatureCountsPerStrand",
        gr2co(list_range),
        "try-error")

    if (retain_coverage) {

        counts <- do.call(rbind, list_counts)
        counts <- counts[order(unlist(list_index)), ]

    } else {

        counts <- unlist(list_counts, use.names = FALSE)
        counts <- counts[order(unlist(list_index))]

    }

    generateCompleteMessage(sample_name)

    return(counts)

}

getSGFeatureCountsPerStrand <- function(features, file_bam, paired_end,
    sample_name, min_anchor, retain_coverage, verbose)
{

    which <- range(features)

    if (length(which) > 1) {

        stop("features must be on the same seqlevel and strand")

    }

    ir <- extractRangesFromFeatures(features)

    ## obtain exons and introns
    gap <- readGap(file_bam, paired_end, which, sample_name, verbose)
    frag_exonic <- gap$frag_exonic
    frag_intron <- gap$frag_intron

    ## extract feature type and spliced boundaries
    type <- mcols(ir)$type
    spliceL <- mcols(ir)$spliceL
    spliceR <- mcols(ir)$spliceR

    i_J <- which(type == "J")
    i_E <- which(type == "E")
    i_S <- which(type %in% c("spliceL", "spliceR"))

    N <- rep(NA_integer_, length(ir))

    if (length(i_J) > 0) {

        N[i_J] <- junctionCompatible(ir[i_J], frag_exonic, frag_intron,
            min_anchor)

    }

    if (length(i_E) > 0) {

        E_index <- exonCompatible(ir[i_E], spliceL[i_E], spliceR[i_E],
            frag_exonic, frag_intron, FALSE)
        N[i_E] <- elementNROWS(E_index)

    }

    if (length(i_S) > 0) {

        N[i_S] <- splicesiteOverlap(ir[i_S],
            sub("splice", "", type[i_S], fixed = TRUE),
            frag_exonic, frag_intron, min_anchor, "unspliced")

    }

    if (retain_coverage) {

        counts <- DataFrame(N = N)
        counts$N_splicesite <- IntegerList(vector("list", nrow(counts)))
        counts$coverage <- RleList(IntegerList(vector("list", nrow(counts))))

        if (length(i_J) > 0) {

            counts$N_splicesite[i_J] <- splicesiteCounts(ir[i_J],
                frag_exonic, frag_intron, min_anchor, "junction", "all")

        }

        if (length(i_E) > 0) {

            counts$N_splicesite[i_E] <- splicesiteCounts(ir[i_E],
                frag_exonic, frag_intron, min_anchor, "exon", "spliced")
            counts$coverage[i_E] <- exonCoverage(ir[i_E], E_index,
                frag_exonic)

        }

        if (as.character(strand(which)) == "-") {

            counts$N_splicesite <- endoapply(counts$N_splicesite, rev)
            counts$coverage <- endoapply(counts$coverage, rev)

        }

    } else {

        counts <- N

    }

    if (verbose) generateCompleteMessage(paste(sample_name, gr2co(which)))

    return(counts)

}

extractRangesFromFeatures <- function(x)
{

    ir <- ranges(x)
    mcols(ir) <- convertSlots(x)
    return(ir)

}

convertSlots <- function(x)
{

    slots <- GenomicRanges:::extraColumnSlotsAsDF(x)
    slots <- slots[c("type", "splice5p", "splice3p")]
    strand <- strand(x)

    i_pos <- which(strand == "+")
    i_neg <- which(strand == "-")

    type <- as.character(slots$type)
    type[i_pos] <- sub("D", "spliceR", type[i_pos], fixed = TRUE)
    type[i_pos] <- sub("A", "spliceL", type[i_pos], fixed = TRUE)
    type[i_neg] <- sub("D", "spliceL", type[i_neg], fixed = TRUE)
    type[i_neg] <- sub("A", "spliceR", type[i_neg], fixed = TRUE)
    slots$type <- type

    spliceL <- rep(NA, nrow(slots))
    spliceR <- rep(NA, nrow(slots))

    spliceL[i_pos] <- slots$splice5p[i_pos]
    spliceL[i_neg] <- slots$splice3p[i_neg]
    spliceR[i_pos] <- slots$splice3p[i_pos]
    spliceR[i_neg] <- slots$splice5p[i_neg]

    slots$splice5p <- NULL
    slots$splice3p <- NULL

    slots$spliceL <- spliceL
    slots$spliceR <- spliceR

    return(slots)

}

##' Create \code{SGFeatureCounts} object from rowRanges, colData and counts.
##'
##' @title Create \code{SGFeatureCounts} object
##' @param rowRanges \code{SGFeatures} object
##' @param colData Data frame with sample information
##' @param counts Integer matrix of counts
##' @param min_anchor Integer specifiying minimum anchor length
##' @return \code{SGFeatureCounts} object
##' @examples
##' sgfc <- makeSGFeatureCounts(sgf_pred, si,
##'   matrix(0L, length(sgf_pred), nrow(si)))
##' @author Leonard Goldstein

makeSGFeatureCounts <- function(rowRanges, colData, counts, min_anchor = 1)
{

    colData <- DataFrame(colData, row.names = colData$sample_name)
    colnames(counts) <- colData$sample_name
    x <- SummarizedExperiment(
        assays = list(counts = counts),
        rowRanges = rowRanges,
        colData = colData)
    assays(x)$FPKM <- getScaledCounts(x, min_anchor)
    x <- SGFeatureCounts(x)

    return(x)

}

getScaledCounts <- function(object, min_anchor,
    features = rowRanges(object),
    counts = assays(object)$counts,
    paired_end = colData(object)$paired_end,
    read_length = colData(object)$read_length,
    frag_length = colData(object)$frag_length,
    lib_size = colData(object)$lib_size)
{

    if (is(features, "SGFeatures")) {

        L <- getFeatureLengths(features, min_anchor)

    } else if (is(features, "SGVariants")) {

        L <- rep(-2 * min_anchor + 2, length(features))

    }

    E <- getEffectiveLengths(L, paired_end, read_length, frag_length)
    X <- counts / (E * 1e-3)
    X <- sweep(X, 2, lib_size * 1e-6, FUN = "/")

    return(X)

}

getFeatureLengths <- function(features, min_anchor)
{

    L <- rep(NA_integer_, length(features))
    i <- which(type(features) %in% c("J", "A", "D"))
    if (length(i) > 0) { L[i] <- -2 * min_anchor + 2 }
    i <- which(type(features) == "E")
    if (length(i) > 0) { L[i] <- width(features[i]) }

    return(L)

}

getEffectiveLengths <- function(L, paired_end, read_length, frag_length)
{

    ## R = read length
    ## F = fragment length
    ## I = distance between fragment ends (I = F - 2 * R)
    ## L = feature length
    ## M = min anchor length
    ## E = effective feature length

    ## Single-end data
    ## (1) if type equals J, D or A,
    ## then E = R - 2 * M + 1
    ## (2) if type equals I, F or L,
    ## then E = R + L - 1
    ## Note (1) is a special case of (2) for L = -2 * M + 2

    ## Paired-end data
    ## (1) if type equals J, D or A,
    ## if I >= -2 * M + 1 : E = 2 * (R - 2 * M + 1)
    ## if I <  -2 * M + 1 : E = F - 2 * M + 1
    ## (2) if type equals I, F or L,
    ## if I >= L - 1 : E = 2 * (R + L - 1)
    ## if I <  L - 1 : E = F + L - 1
    ## Note (1) is a special case of (2) for L = -2 * M + 2

    n_features <- length(L)
    n_samples <- length(paired_end)

    E <- matrix(NA_real_, nrow = n_features, ncol = n_samples)

    i_PE <- which(paired_end)
    n_PE <- length(i_PE)

    if (n_PE > 0) {

        L_PE <- rep(L, n_PE)
        R_PE <- rep(read_length[i_PE], rep(n_features, n_PE))
        F_PE <- rep(frag_length[i_PE], rep(n_features, n_PE))
        I_PE <- F_PE - 2 * R_PE
        E_PE <- L_PE + F_PE - 1 - pmax(I_PE - L_PE + 1, 0)
        E[, i_PE] <- matrix(E_PE, nrow = n_features, ncol = n_PE)

    }

    i_SE <- which(!paired_end)
    n_SE <- length(i_SE)

    if (n_SE > 0) {

        L_SE <- rep(L, n_SE)
        R_SE <- rep(read_length[i_SE], rep(n_features, n_SE))
        E_SE <- L_SE + R_SE - 1
        E[, i_SE] <- matrix(E_SE, nrow = n_features, ncol = n_SE)

    }

    return(E)

}

convertFpkmToCount <- function(fpkm, paired_end, read_length, frag_length,
    lib_size, feature_length = 0)
{

    if (paired_end) {

        I <- frag_length - 2 * read_length
        E <- feature_length + frag_length - 1 - max(I - feature_length + 1, 0)

    } else {

        E <- feature_length + read_length - 1

    }

    count <- fpkm * (lib_size * 1e-6) * (E * 1e-3)

    return(count)

}

getSGVariantCountsFromSGFeatureCounts <- function(variants, sgfc,
    min_denominator, cores)
{

    n_variants <- length(variants)
    n_sample <- ncol(sgfc)

    if (n_variants == 0) {

        return(SGVariantCounts())

    }

    v5p <- featureID5p(variants)
    v3p <- featureID3p(variants)
    e5p <- featureID5pEvent(variants)
    e3p <- featureID3pEvent(variants)

    sgfc <- sgfc[featureID(sgfc) %in% unique(c(unlist(e5p), unlist(e3p))), ]
    features <- rowRanges(sgfc)

    variant_i_5p <- relist(match(unlist(v5p), featureID(features)), v5p)
    variant_i_3p <- relist(match(unlist(v3p), featureID(features)), v3p)
    event_i_5p <- relist(match(unlist(e5p), featureID(features)), e5p)
    event_i_3p <- relist(match(unlist(e3p), featureID(features)), e3p)

    x <- assay(sgfc, "counts")

    variant_x_5p <- collapseRows(x, variant_i_5p, cores = cores)
    variant_x_3p <- collapseRows(x, variant_i_3p, cores = cores)
    event_x_5p <- collapseRows(x, event_i_5p, cores = cores)
    event_x_3p <- collapseRows(x, event_i_3p, cores = cores)

    assays <- list(
        "countsVariant5p" = variant_x_5p,
        "countsVariant3p" = variant_x_3p,
        "countsEvent5p" = event_x_5p,
        "countsEvent3p" = event_x_3p)

    object <- SummarizedExperiment(
        assays = assays,
        rowRanges = variants,
        colData = colData(sgfc))
    colnames(object) <- colnames(sgfc)
    object <- getVariantFreq(object, min_denominator)
    object <- SGVariantCounts(object)

    return(object)

}

collapseRows <- function(x, list_i, fun = sum, cores = 1)
{

    y <- matrix(NA_integer_, nrow = length(list_i), ncol = ncol(x))

    j <- which(elementNROWS(list_i) == 1)

    if (length(j) > 0) {

        i <- unlist(list_i[j])
        y[j, ] <- x[i, ]

    }

    j <- which(elementNROWS(list_i) > 1)

    if (length(j) > 0) {

        i <- unlist(list_i[j])
        f <- togroup0(list_i[j])
        y[j, ] <- do.call(cbind, mclapply(seq_len(ncol(x)),
            function(j) { tapply(x[i, j, drop = FALSE], f, fun) },
            mc.cores = cores))

    }

    colnames(y) <- colnames(x)

    return(y)

}

getVariantFreq <- function(object, min_denominator)
{

    U5p <- assay(object, "countsVariant5p")
    U3p <- assay(object, "countsVariant3p")
    V5p <- assay(object, "countsEvent5p")
    V3p <- assay(object, "countsEvent3p")

    variants <- rowRanges(object)

    v5p <- featureID5p(variants)
    v3p <- featureID3p(variants)

    e5p <- featureID5pEvent(variants)
    e3p <- featureID3pEvent(variants)

    d5p <- elementNROWS(v5p) > 0 & elementNROWS(e5p) > 0
    d3p <- elementNROWS(v3p) > 0 & elementNROWS(e3p) > 0

    ixp <- which(d5p & d3p)
    i5p <- which(d5p & !d3p)
    i3p <- which(!d5p & d3p)

    U <- matrix(NA_integer_, nrow = nrow(object), ncol = ncol(object))
    U[ixp, ] <- U5p[ixp, ] + U3p[ixp, ]
    U[i5p, ] <- U5p[i5p, ]
    U[i3p, ] <- U3p[i3p, ]

    V <- matrix(NA_integer_, nrow = nrow(object), ncol = ncol(object))
    V[ixp, ] <- V5p[ixp, ] + V3p[ixp, ]
    V[i5p, ] <- V5p[i5p, ]
    V[i3p, ] <- V3p[i3p, ]

    X <- U/V
    X[is.na(X)] <- NA_real_

    if (!is.na(min_denominator)) {

        I <- matrix(FALSE, nrow = nrow(object), ncol = ncol(object))
        I[ixp, ] <- V5p[ixp, ] < min_denominator & V3p[ixp, ] < min_denominator
        I[i5p, ] <- V5p[i5p, ] < min_denominator
        I[i3p, ] <- V3p[i3p, ] < min_denominator
        X[which(I)] <- NA_real_

    }

    assay(object, "variantFreq") <- X

    return(object)

}

getSGVariantCountsFromBamFiles <- function(variants, sample_info,
    min_denominator, min_anchor, counts_only = FALSE, verbose, cores)
{

    checkSampleInfo(sample_info)

    features <- createSGFeaturesFromSGVariants(variants)

    cores <- setCores(cores, nrow(sample_info))

    list_counts <- mcmapply(
        getSGVariantCountsPerSample,
        file_bam = sample_info$file_bam,
        paired_end = sample_info$paired_end,
        sample_name = sample_info$sample_name,
        MoreArgs = list(
            variants = variants,
            features = features,
            min_anchor = min_anchor,
            verbose = verbose,
            cores = cores$per_sample),
        SIMPLIFY = FALSE,
        USE.NAMES = FALSE,
        mc.preschedule = setPreschedule(cores$n_sample),
        mc.cores = cores$n_sample
    )

    checkApplyResultsForErrors(
        list_counts,
        "getSGVariantCountsPerSample",
        sample_info$sample_name,
        "character")

    if (counts_only) return(list_counts)

    sgvc <- makeSGVariantCounts(variants, sample_info, list_counts,
        min_denominator, cores$total)

    return(sgvc)

}

createSGFeaturesFromSGVariants <- function(variants)
{

    features <- unlist(variants)
    features <- features[!duplicated(featureID(features))]
    features <- addSpliceSites(features, variants, "D")
    features <- addSpliceSites(features, variants, "A")

    return(features)

}

addSpliceSites <- function(features, variants, type = c("D", "A"))
{

    type <- match.arg(type)

    if (type == "D") {

        fid <- featureID5pEvent(variants)

    } else if (type == "A") {

        fid <- featureID3pEvent(variants)

    }

    fid_unlisted <- unlist(fid)
    grp <- factor(togroup0(fid), levels = seq_along(variants))
    i <- which(!fid_unlisted %in% featureID(features))
    fid <- split(fid_unlisted[i], grp[i])

    if (any(elementNROWS(fid) > 1)) {

        stop("invalid variants")

    }

    i <- which(elementNROWS(fid) == 1)

    if (length(i) == 0) {

        return(features)

    }

    if (type == "D") {

        splicesites <- pos2gr(sub("^D:", "", from(variants)[i]))

    } else if (type == "A") {

        splicesites <- pos2gr(sub("^A:", "", to(variants)[i]))

    }

    mcols(splicesites)$type <- rep(type, length(splicesites))
    mcols(splicesites)$splice5p <- rep(NA, length(splicesites))
    mcols(splicesites)$splice3p <- rep(NA, length(splicesites))
    mcols(splicesites)$featureID <- unlist(fid[i], use.names = FALSE)
    mcols(splicesites)$geneID <- geneID(variants)[i]

    splicesites <- SGFeatures(splicesites)
    splicesites <- splicesites[!duplicated(featureID(splicesites))]

    new2old <- match(seqlevels(features), seqlevels(splicesites))
    seqinfo(splicesites, new2old = new2old) <- seqinfo(features)
    
    features <- c(features, splicesites)

    return(features)

}

getSGVariantCountsPerSample <- function(variants, features, file_bam,
    paired_end, sample_name, min_anchor, verbose, cores)
{

    variants_range <- unlist(range(variants))
    list_range <- range(variants_range) + 1

    hits_1 <- findOverlaps(variants_range, list_range)
    list_index <- split(queryHits(hits_1), subjectHits(hits_1))
    list_variants <- split(variants[queryHits(hits_1)], subjectHits(hits_1))

    hits_2 <- findOverlaps(features, list_range)
    list_features <- split(features[queryHits(hits_2)], subjectHits(hits_2))

    list_counts <- mcmapply(
        getSGVariantCountsPerStrand,
        variants = list_variants,
        features = list_features,
        MoreArgs = list(
            file_bam = file_bam,
            paired_end = paired_end,
            sample_name = sample_name,
            min_anchor = min_anchor,
            verbose = verbose),
        SIMPLIFY = FALSE,
        USE.NAMES = FALSE,
        mc.preschedule = setPreschedule(cores),
        mc.cores = cores)

    checkApplyResultsForErrors(
        list_counts,
        "getSGVariantCountsPerStrand",
        gr2co(list_range),
        "character")

    counts <- do.call(rbind, list_counts)
    counts <- counts[order(unlist(list_index)), ]

    generateCompleteMessage(sample_name)

    return(counts)

}

getSGVariantCountsPerStrand <- function(variants, features,
    file_bam, paired_end, sample_name, min_anchor, verbose)
{

    which <- range(unlist(variants))

    if (length(which) > 1) {

        stop("variants must be on the same seqlevel and strand")

    }

    vid <- factor(variantID(variants))
    v5p <- featureID5p(variants)
    v3p <- featureID3p(variants)
    vxp <- unique(pc(v5p, v3p))
    e5p <- featureID5pEvent(variants)
    e3p <- featureID3pEvent(variants)

    ## in the case of nested variants with representative features
    ## at both 5' and 3' end, consider features at the 3' end only
    ## (for nested variant, counts based on both ends can be affected not
    ## only by changes in the relative usage of the variant, but also by
    ## changes in the relative usage of the intra-variant nested variants)

    i <- which(elementNROWS(v5p) > 0 & elementNROWS(v3p) > 0 &
        grepl("(", featureID(variants), fixed = TRUE))

    if (length(i) > 0) {

        vxp[i] <- v3p[i]

    }

    ## extract ranges for relevant features
    fid <- sort(unique(c(unlist(e5p), unlist(e3p))))
    features <- features[match(fid, featureID(features))]
    ir <- extractRangesFromFeatures(features)

    ## obtain exons and introns
    gap <- readGap(file_bam, paired_end, which, sample_name, verbose)
    frag_exonic <- gap$frag_exonic
    frag_intron <- gap$frag_intron

    ## extract feature type and spliced boundaries
    type <- mcols(ir)$type
    spliceL <- mcols(ir)$spliceL
    spliceR <- mcols(ir)$spliceR

    ir_index <- IntegerList(vector("list", length(ir)))

    i_J <- which(type == "J")

    if (length(i_J) > 0) {

        ir_index[i_J] <- junctionCompatible(ir[i_J], frag_exonic,
            frag_intron, min_anchor, FALSE)

    }

    i_S <- which(type %in% c("spliceL", "spliceR"))

    if (length(i_S) > 0) {

        ir_index[i_S] <- splicesiteOverlap(ir[i_S],
            sub("splice", "", type[i_S], fixed = TRUE),
            frag_exonic, frag_intron, min_anchor, "unspliced", FALSE)

    }

    assay_names <- c(
        "countsVariant5p",
        "countsVariant3p",
        "countsEvent5p",
        "countsEvent3p",
        "countsVariant5pOr3p")

    counts <- matrix(NA_integer_, nrow = length(vid),
        ncol = length(assay_names))
    colnames(counts) <- assay_names

    list_opt <- c("Variant5p", "Variant3p", "Event5p", "Event3p",
        "Variant5pOr3p")

    for (opt in list_opt) {

        f <- switch(opt, Variant5p = v5p, Variant3p = v3p,
            Event5p = e5p, Event3p = e3p, Variant5pOr3p = vxp)
        i <- ir_index[match(unlist(f), featureID(features))]
        v <- vid[togroup0(f)]
        tmp <- IntegerList(tapply(unlist(i), v[togroup0(i)], unique,
            simplify = FALSE))
        x <- elementNROWS(tmp)[match(vid, names(tmp))]
        x[elementNROWS(f) == 0] <- NA_integer_
        counts[, paste0("counts", opt)] <- x

    }

    if (verbose) generateCompleteMessage(paste(sample_name, gr2co(which)))

    return(counts)

}

makeSGVariantCounts <- function(rowRanges, colData, counts, min_denominator,
    cores)
{

    assays <- list()

    for (k in colnames(counts[[1]])) {

        a <- do.call(cbind, mclapply(counts, function(x) { x[, k] },
            mc.cores = cores))
        colnames(a) <- colData$sample_name
        assays[[k]] <- a

    }

    colData <- DataFrame(colData, row.names = colData$sample_name)
    x <- SummarizedExperiment(
        assays = assays,
        rowRanges = rowRanges,
        colData = colData)
    x <- getVariantFreq(x, min_denominator)
    x <- SGVariantCounts(x)

    return(x)

}
