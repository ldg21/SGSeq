##' Obtain counts of compatible fragments for splice graph features.
##'
##' @title Compatible fragment counts for splice graph features
##' @inheritParams predictTxFeaturesPerSample
##' @param features \code{SGFeatures} object
##' @return Numeric vector of compatible fragment counts
##' @keywords internal
##' @author Leonard Goldstein

getSGFeatureCountsPerSample <- function(features, file_bam, paired_end,
    retain_coverage, sample_name, verbose, cores)
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
        retain_coverage = retain_coverage,
        sample_name = sample_name,
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
    
        counts <- unlist(list_counts)
        counts <- counts[order(unlist(list_index))]
        
    }    

    generateCompleteMessage(sample_name)
    
    return(counts)
    
}

getSGFeatureCountsPerStrand <- function(features, file_bam, paired_end,
    retain_coverage, sample_name, verbose)
{

    which <- range(features)

    if (length(which) > 1) {

        stop("features must be on the same seqlevel and strand")

    }

    seqlevel <- as.character(seqnames(which))
    strand <- as.character(strand(which))

    ir <- extractRangesFromFeatures(features)
    
    ## obtain GAlignmentPairs
    gap <- readGap(file_bam, paired_end, which)
    gap <- gap[mcols(gap)$strand %in% c(strand, "*")]
    
    ## get exons and introns
    frag_exonic <- reduce(ranges(grglist(gap, drop.D.ranges = TRUE)))
    frag_intron <- ranges(junctions(gap))
    
    ## extract feature type and spliced boundaries
    type <- mcols(ir)$type
    spliceL <- mcols(ir)$spliceL
    spliceR <- mcols(ir)$spliceR

    i_J <- which(type == "J")
    i_E <- which(type == "E")
    i_S <- which(type %in% c("spliceL", "spliceR"))

    N <- rep(NA_integer_, length(ir))
    
    if (length(i_J) > 0) {

        N[i_J] <- junctionCompatible(ir[i_J], frag_intron)

    }
    
    if (length(i_E) > 0) {
        
        E_index <- exonCompatible(ir[i_E], spliceL[i_E], spliceR[i_E],
            frag_exonic, frag_intron, FALSE)
        N[i_E] <- elementLengths(E_index)

    }
    
    if (length(i_S) > 0) {
        
        N[i_S] <- splicesiteOverlap(ir[i_S],
            sub("splice", "", type[i_S], fixed = TRUE),
            frag_exonic, frag_intron, "unspliced")

    }

    if (retain_coverage) {

        counts <- DataFrame(N = N)
        counts$N_splicesite <- IntegerList(as.integer())
        counts$coverage <- RleList(as.integer(), compress = TRUE)
        
        if (length(i_J) > 0) {
        
            counts$N_splicesite[i_J] <- splicesiteCounts(ir[i_J],
                frag_exonic, frag_intron, "junction")

        }

        if (length(i_E) > 0) {
            
            counts$N_splicesite[i_E] <- splicesiteCounts(ir[i_E],
                frag_exonic, frag_intron, "exon")
            counts$coverage[i_E] <- exonCoverage(ir[i_E],
                E_index, frag_exonic)
            
        }

        if (strand == "-") {

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
##' @title Create \code{SGFeatureCounts} object
##' @param rowRanges An \code{SGFeatures} object
##' @param colData Data frame with sample information 
##' @param counts Integer matrix of counts
##' @return An \code{SGFeatureCounts} object
##' @examples
##' sgfc <- makeSGFeatureCounts(sgf_pred, si,
##'   matrix(0L, length(sgf_pred), nrow(si)))
##' @author Leonard Goldstein

makeSGFeatureCounts <- function(rowRanges, colData, counts)
{

    colnames(counts) <- colData$sample_name    
    x <- SummarizedExperiment(
        assays = list(counts = counts),
        rowRanges = rowRanges,
        colData = DataFrame(colData))
    colnames(x) <- colData$sample_name
    x <- getScaledCounts(x)
    x <- SGFeatureCounts(x)
    
    return(x)
    
}

getScaledCounts <- function(SE,
    paired_end = colData(SE)$paired_end,
    read_length = colData(SE)$read_length,
    frag_length = colData(SE)$frag_length,
    lib_size = colData(SE)$lib_size)
{

    L <- getFeatureLengths(rowRanges(SE))
    E <- getEffectiveLengths(L, paired_end, read_length, frag_length)
    X <- assay(SE, "counts")
    X <- X / (E * 1e-3)
    X <- sweep(X, 2, lib_size * 1e-6, FUN = "/")    
    assays(SE)$FPKM <- X

    return(SE)
    
}

getFeatureLengths <- function(features)
{

    if (is(features, "SGFeatures")) {
        
        L <- rep(NA_integer_, length(features))
        i <- which(type(features) %in% c("J", "A", "D"))
        if (length(i) > 0) { L[i] <- 0 }        
        i <- which(type(features) == "E")
        if (length(i) > 0) { L[i] <- width(features[i]) }

    } else if (is(features, "SGSegments")) {

        features_unlisted <- unlist(features)
        i <- which(type(features_unlisted) == "E")
        exons <- split(features_unlisted[i], togroup(features)[i])
        L <- rep(0, length(features))
        L[as.integer(names(exons))] <- unlist(sum(width(exons)))
        
    }

    return(L)
    
}

getEffectiveLengths <- function(L, paired_end, read_length, frag_length)
{

    ## R = read length
    ## F = fragment length
    ## I = distance between fragment ends (I = F - 2 * R)
    ## L = feature length
    ## E = effective feature length

    ## Single-end data
    ## (1) if type equals J, D or A,
    ## then E = R - 1
    ## (2) if type equals I, F or L,
    ## then E = L + R - 1
    ## Thus (1) is a special case of (2) for L = 0
    
    ## Paired-end data
    ## (1) if type equals J, D or A,
    ## then E = 2 x (R - 1) + min(I + 1, 0)
    ## (2) if type equals I, F or L,
    ## then E = L + F - 1 - pmax(I - L + 1, 0)
    ## Thus (1) is a special case of (2) for L = 0

    n_feature <- length(L)
    n_sample <- length(paired_end)

    E <- matrix(NA_real_, nrow = n_feature, ncol = n_sample)

    i_PE <- which(paired_end)
    n_PE <- length(i_PE)
    
    if (n_PE > 0) {

        L_PE <- rep(L, n_PE)
        R_PE <- rep(read_length[i_PE], rep(n_feature, n_PE))
        F_PE <- rep(frag_length[i_PE], rep(n_feature, n_PE))
        I_PE <- F_PE - 2 * R_PE
        E_PE <- L_PE + F_PE - 1 - pmax(I_PE - L_PE + 1, 0)
        E[, i_PE] <- matrix(E_PE, nrow = n_feature, ncol = n_PE)

    }

    i_SE <- which(!paired_end)
    n_SE <- length(i_SE)

    if (n_SE > 0) {
    
        L_SE <- rep(L, n_SE)
        R_SE <- rep(read_length[i_SE], rep(n_feature, n_SE))
        E_SE <- L_SE + R_SE - 1
        E[, i_SE] <- matrix(E_SE, nrow = n_feature, ncol = n_SE)
        
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

getSGVariantCountsFromSGFeatureCounts <- function(variants, object,
    min_denominator, cores)
{

    n_variants <- length(variants)
    n_sample <- ncol(object)

    if (n_variants == 0) {

        return(SGVariantCounts())

    }

    featureIDs <- union(
        unlist(featureID5p(variants)), unlist(featureID3p(variants)))
    
    if (!all(featureIDs %in% featureID(object))) {

        stop("'object' is missing features included in 'variants'")

    }

    object <- object[featureID(object) %in% featureIDs, ]
    features <- rowRanges(object)
    
    variant_i_start <- relist(match(unlist(featureID5p(variants)),
        featureID(features)), featureID5p(variants))
    variant_i_end <- relist(match(unlist(featureID3p(variants)),
        featureID(features)), featureID3p(variants))

    event_i_start <- tapply(unlist(variant_i_start),
        eventID(variants)[togroup(variant_i_start)], unique, simplify = FALSE)
    event_i_end <- tapply(unlist(variant_i_end),
        eventID(variants)[togroup(variant_i_end)], unique, simplify = FALSE)

    event_uninformative_start <- names(which(tapply(
        elementLengths(variant_i_start) == 0, eventID(variants), any)))
    i <- which(names(event_i_start) %in% event_uninformative_start)
    event_i_start[i] <- vector("list", length(i))
    event_uninformative_end <- names(which(tapply(
        elementLengths(variant_i_end) == 0, eventID(variants), any)))
    i <- which(names(event_i_end) %in% event_uninformative_end)
    event_i_end[i] <- vector("list", length(i))
    
    x <- assay(object, "counts")

    variant_x_start <- collapseRows(x, variant_i_start, cores = cores)
    variant_x_end <- collapseRows(x, variant_i_end, cores = cores)
    
    event_x_start <- collapseRows(x, event_i_start, cores = cores)
    event_x_end <- collapseRows(x, event_i_end, cores = cores)
    
    i <- match(eventID(variants), names(event_i_start))
    variant_x_start_total <- event_x_start[i, , drop = FALSE]
    i <- match(eventID(variants), names(event_i_end))
    variant_x_end_total <- event_x_end[i, , drop = FALSE]

    assays <- list(
        "countsVariant5p" = variant_x_start,
        "countsTotal5p" = variant_x_start_total,
        "countsVariant3p" = variant_x_end,
        "countsTotal3p" = variant_x_end_total)

    SE <- SummarizedExperiment(
        assays = assays,
        rowRanges = variants,
        colData = colData(object))
    colnames(SE) <- colnames(object)
    SE <- getVariantFreq(SE, min_denominator)
    SE <- SGVariantCounts(SE)
    
    return(SE)
    
}

collapseRows <- function(x, list_i, fun = sum, cores = 1)
{

    y <- matrix(NA_integer_, nrow = length(list_i), ncol = ncol(x))
    
    j <- which(elementLengths(list_i) == 1)

    if (length(j) > 0) {

        i <- unlist(list_i[j])
        y[j, ] <- x[i, ]
        
    }

    j <- which(elementLengths(list_i) > 1)

    if (length(j) > 0) {

        i <- unlist(list_i[j])
        f <- togroup(list_i[j])
        y[j, ] <- do.call(cbind, mclapply(seq_len(ncol(x)),
            function(j) { tapply(x[i, j, drop = FALSE], f, fun) },
            mc.cores = cores))

    }

    colnames(y) <- colnames(x)
    
    return(y)
    
}

getVariantFreq <- function(SE, min_denominator)
{
    
    U5p <- assay(SE, "countsVariant5p")
    V5p <- assay(SE, "countsTotal5p")
    U3p <- assay(SE, "countsVariant3p")
    V3p <- assay(SE, "countsTotal3p")

    variants <- rowRanges(SE)
    eid <- eventID(variants)
    f5p <- featureID5p(variants)
    f3p <- featureID3p(variants)

    e5p <- names(which(tapply(elementLengths(f5p) > 0, eid, all)))
    e3p <- names(which(tapply(elementLengths(f3p) > 0, eid, all)))

    ixp <- which(eid %in% e5p & eid %in% e3p)
    i5p <- which(eid %in% e5p & !eid %in% e3p)
    i3p <- which(!eid %in% e5p & eid %in% e3p)

    U <- matrix(NA_integer_, nrow = nrow(SE), ncol = ncol(SE))
    U[ixp, ] <- U5p[ixp, ] + U3p[ixp, ]
    U[i5p, ] <- U5p[i5p, ]
    U[i3p, ] <- U3p[i3p, ]

    V <- matrix(NA_integer_, nrow = nrow(SE), ncol = ncol(SE))
    V[ixp, ] <- V5p[ixp, ] + V3p[ixp, ]
    V[i5p, ] <- V5p[i5p, ]
    V[i3p, ] <- V3p[i3p, ]

    X <- U/V
    X[is.na(X)] <- NA_real_

    if (!is.na(min_denominator)) {

        X[V < min_denominator] <- NA_real_

    }
    
    assay(SE, "variantFreq") <- X
        
    return(SE)
    
}

getSGVariantCountsFromBamFiles <- function(variants, features, sample_info,
    counts_only = FALSE, min_denominator, verbose, cores)
{

    checkSampleInfo(sample_info)

    featureIDs <- union(
        unlist(featureID5p(variants)), unlist(featureID3p(variants)))
    
    if (!all(featureIDs %in% featureID(features))) {

        stop("'features' is missing features included in 'variants'")

    }

    features <- features[featureID(features) %in% featureIDs]
    
    cores <- setCores(cores, nrow(sample_info))

    list_counts <- mcmapply(
        getSGVariantCountsPerSample,
        file_bam = sample_info$file_bam,
        paired_end = sample_info$paired_end,
        sample_name = sample_info$sample_name,
        MoreArgs = list(
            variants = variants,
            features = features,
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

getSGVariantCountsPerSample <- function(variants, features,
    file_bam, paired_end, sample_name, verbose, cores)
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
    file_bam, paired_end, sample_name, verbose)
{

    which <- range(unlist(variants))

    if (length(which) > 1) {

        stop("variants must be on the same seqlevel and strand")

    }

    seqlevel <- as.character(seqnames(which))
    strand <- as.character(strand(which))
    
    vid <- factor(variantID(variants))
    eid <- factor(eventID(variants))

    f5p <- featureID5p(variants)
    f3p <- featureID3p(variants)
    all <- unique(pc(f5p, f3p))
    
    ## in the case of nested variants with representative features
    ## at both 5' and 3' end, consider features at the 3' end only
    ## (for nested variant, counts based on both ends can be
    ## affected not only by changes in the relative usage of
    ## the variant under consideration, but also by changes in
    ## the relative usage of the intra-variant nested variants)
    i <- which(elementLengths(f5p) > 0 & elementLengths(f3p) > 0 &
        grep("(", featureID(variants), fixed = TRUE))
    if (length(i) > 0) all[i] <- f3p[i]

    ## extract ranges for relevant features
    fid <- sort(unique(c(unlist(f5p), unlist(f3p))))
    features <- features[match(fid, featureID(features))]
    ir <- extractRangesFromFeatures(features)
    
    ## obtain GAlignmentPairs
    gap <- readGap(file_bam, paired_end, which)
    gap <- gap[mcols(gap)$strand %in% c(strand, "*")]
    
    ## get exons and introns
    frag_exonic <- reduce(ranges(grglist(gap, drop.D.ranges = TRUE)))
    frag_intron <- ranges(junctions(gap))
    
    ## extract feature type and spliced boundaries
    type <- mcols(ir)$type
    spliceL <- mcols(ir)$spliceL
    spliceR <- mcols(ir)$spliceR

    i_J <- which(type == "J")
    i_S <- which(type %in% c("spliceL", "spliceR"))

    ir_index <- as(vector("list", length(ir)), "CompressedIntegerList")
    
    if (length(i_J) > 0) {

        ir_index[i_J] <- junctionCompatible(ir[i_J], frag_intron, FALSE)

    }
    
    if (length(i_S) > 0) {
        
        ir_index[i_S] <- splicesiteOverlap(ir[i_S],
            sub("splice", "", type[i_S], fixed = TRUE),
            frag_exonic, frag_intron, "unspliced", FALSE)

    }

    assay_names <- c(
        "countsVariant5p",
        "countsTotal5p",
        "countsVariant3p",
        "countsTotal3p",
        "countsVariant",
        "countsTotal")
    
    counts <- matrix(NA_integer_, nrow = length(vid), ncol = 6)
    colnames(counts) <- assay_names
    
    for (opt in c("5p", "3p", "all")) {

        f <- switch(opt, "5p" = f5p, "3p" = f3p, "all" = all)
        s <- switch(opt, "5p" = "5p", "3p" = "3p", "all" = "")
        
        i <- match(unlist(f), featureID(features))
        g <- togroup(f)

        tmp <- as(tapply(unlist(ir_index[i]),
            vid[g][togroup(ir_index[i])], unique,
            simplify = FALSE), "CompressedIntegerList")
        countsVariant <- elementLengths(tmp)[match(vid, names(tmp))]
        vid_uninformative <- vid[elementLengths(f) == 0]
        countsVariant[vid %in% vid_uninformative] <- NA_integer_
        counts[, paste0("countsVariant", s)] <- countsVariant
        
        tmp <- as(tapply(unlist(ir_index[i]),
            eid[g][togroup(ir_index[i])], unique,
            simplify = FALSE), "CompressedIntegerList")
        countsTotal <- elementLengths(tmp)[match(eid, names(tmp))]
        eid_uninformative <- names(which(tapply(
            elementLengths(f) == 0, eid, any)))
        countsTotal[eid %in% eid_uninformative] <- NA_integer_
        counts[, paste0("countsTotal", s)] <- countsTotal

    }
    
    if (verbose) generateCompleteMessage(paste(sample_name, gr2co(which)))
    
    return(counts)
    
}

makeSGVariantCounts <- function(rowRanges, colData, counts, min_denominator,
    cores)
{

    assays <- list()

    for (k in colnames(counts[[1]])) {

        assays[[k]] <- do.call(cbind,
            mclapply(counts, function(x) { x[, k] }, mc.cores = cores))
      
    }
        
    x <- SummarizedExperiment(
        assays = assays,
        rowRanges = rowRanges,
        colData = DataFrame(colData))
    colnames(x) <- colData$sample_name
    x <- getVariantFreq(x, min_denominator)
    x <- SGVariantCounts(x)

    return(x)
    
}
