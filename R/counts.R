##' Obtain counts of compatible fragments for splice graph features.
##'
##' @title Compatible fragment counts for splice graph features
##' @inheritParams predictTxFeaturesPerSample
##' @param features \code{SGFeatures} object
##' @return Numeric vector of compatible fragment counts
##' @keywords internal
##' @author Leonard Goldstein

getSGFeatureCountsPerSample <- function(features, file_bam, paired_end,
    retain_coverage = FALSE, cores = 1)
{

    if (!is(features, "SGFeatures")) {

        stop("argument features must be an SGFeatures object")

    }
    
    hits <- findOverlaps(features, range(features))
    list_index <- split(queryHits(hits), subjectHits(hits))
    list_features <- split(features[queryHits(hits)], subjectHits(hits))
            
    list_counts <- mclapply(
        list_features,
        getSGFeatureCountsRanges,
        file_bam = file_bam,
        paired_end = paired_end,
        retain_coverage = retain_coverage, 
        mc.preschedule = FALSE,
        mc.cores = cores)

    if (retain_coverage) {

        counts <- do.call(rbind, list_counts)
        counts <- counts[order(unlist(list_index)), ]
        
    } else {
    
        counts <- unlist(list_counts)
        counts <- counts[order(unlist(list_index))]
        
    }    

    return(counts)
    
}

getSGFeatureCountsRanges <- function(features, file_bam, paired_end,
    retain_coverage)
{

    which <- range(as(features, "GRanges"))

    if (length(which) > 1) {

        stop("features must be on the same seqlevel and strand")

    }

    seqlevel <- as.character(seqnames(which))
    strand <- as.character(strand(which))

    ir <- extractRangesFromFeatures(features)
    
    ## obtain GAlignmentPairs
    gap <- readGap(file_bam, paired_end, which)
    gap <- gap[strand(gap) %in% c(strand, "*")]
    
    ## get exons and introns
    frag_exonic <- ranges(grglist(gap, drop.D.ranges = TRUE))
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

makeSGFeatureCounts <- function(rowData, colData, counts)
{

    colnames(counts) <- colData$sample_name    
    x <- SummarizedExperiment(assays = list(counts = counts),
        rowData = rowData, colData = DataFrame(colData))
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

    E <- getEffectiveLengths(rowData(SE),
        paired_end, read_length, frag_length)    

    X <- assay(SE, "counts")
    X <- X / (E * 1e-3)
    X <- sweep(X, 2, lib_size * 1e-6, FUN = "/")
    
    assays(SE)$FPKM <- X

    return(SE)
    
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

getEffectiveLengths <- function(features, paired_end, read_length, frag_length)
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

    if (is(features, "SGFeatures")) {
                        
        L <- rep(NA, length(features))
        i <- which(type(features) %in% c("J", "A", "D"))
        if (length(i) > 0) { L[i] <- 0 }        
        i <- which(type(features) == "E")
        if (length(i) > 0) { L[i] <- width(features[i]) }

    }
    if (is(features, "TxSegments")) {

        features_unlisted <- unlist(features)
        i <- which(type(features_unlisted) == "E")
        exons <- split(features_unlisted[i], togroup(features)[i])
        L <- rep(0, length(features))
        L[as.integer(names(exons))] <- unlist(sum(width(exons)))
        
    }

    n_features <- length(features)
    n_samples <- length(paired_end)

    E <- matrix(NA, nrow = n_features, ncol = n_samples)

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

##' For transcript variants, obtain counts of compatible fragments
##' extending across the start and/or end of each variant.
##' Variant frequencies are estimated based on representive counts.
##' 
##' @title Representative counts and frequency estimates for
##'   transcript variants
##' @param object \code{SGFeatureCounts} object
##' @param variants \code{TxVariants} object
##' @return A \code{TxVariantCounts} object
##' @examples
##' txvc <- getTxVariantCounts(sgfc, txv)
##' @author Leonard Goldstein

getTxVariantCounts <- function(object, variants)
{
    
    n_variants <- length(variants)
    n_samples <- ncol(object)
    features <- rowData(object)

    if (n_variants == 0) {

        return(TxVariantCounts())

    }
    
    variant_i_start <- relist(match(unlist(featureID5p(variants)),
        featureID(features)), featureID5p(variants))
    variant_i_end <- relist(match(unlist(featureID3p(variants)),
        featureID(features)), featureID3p(variants))

    event_i_start <- tapply(unlist(variant_i_start),
        eventID(variants)[togroup(variant_i_start)], unique)
    event_i_end <- tapply(unlist(variant_i_end),
        eventID(variants)[togroup(variant_i_end)], unique)

    x <- assay(object, "counts")

    variant_x_start <- collapseRows(x, variant_i_start)
    variant_x_end <- collapseRows(x, variant_i_end)
    
    event_x_start <- collapseRows(x, event_i_start)
    event_x_end <- collapseRows(x, event_i_end)
    
    i <- match(eventID(variants), names(event_i_start))
    variant_x_start_total <- event_x_start[i, , drop = FALSE]
    i <- match(eventID(variants), names(event_i_end))
    variant_x_end_total <- event_x_end[i, , drop = FALSE]
    
    assays <- list(
        "countsVariant5p" = variant_x_start,
        "countsTotal5p" = variant_x_start_total,
        "countsVariant3p" = variant_x_end,
        "countsTotal3p" = variant_x_end_total)

    SE <- SummarizedExperiment(assays = assays, rowData = variants,
        colData = colData(object))
    colnames(SE) <- colnames(object)
    SE <- getVariantFreq(SE)
    SE <- TxVariantCounts(SE)
    
    return(SE)
    
}

representativeFeatures <- function(variant_info, features,
    node = c("start", "end"))
{

    node <- match.arg(node)

    if (node == "start") {

        variant_node <- variant_info$from
        variant_informative <- variant_info$closed3p
        pattern_incl <- "^[^R]"
        pattern_single_feature <- "^[^\\(]"
        pattern_multiple_features <- "^\\("
        pattern_remainder <- ",\\S+$"

    }
    if (node == "end") {

        variant_node <- variant_info$to
        variant_informative <- variant_info$closed5p
        pattern_incl <- "^[^L]"
        pattern_single_feature <- "[^\\)]$"
        pattern_multiple_features <- "\\)$"
        pattern_remainder <- "^\\S+,"

    }

    variant_all_id <- variant_info$featureID
    variant_rep_id <- vector("list", nrow(variant_info))

    ## for each variant identify first or last feature(s) in the variant
    
    i <- grep(pattern_single_feature, variant_all_id)
    i <- i[grep(pattern_incl, variant_node[i])]
    
    if (length(i) > 0) {
    
        variant_rep_id[i] <- as.list(as.integer(
            sub(pattern_remainder, "", variant_all_id[i])))

    }

    i <- grep(pattern_multiple_features, variant_all_id)
    i <- i[grep(pattern_incl, variant_node[i])]

    if (length(i) > 0) {
    
        variant_rep_id[i] <- lapply(expandPath(variant_all_id[i]),
            function(p) { unique(as.integer(sub(pattern_remainder, "", p))) })

    }

    ## replace exons with splice sites

    index <- which(elementLengths(variant_rep_id) > 0)
    tmp_id <- variant_rep_id[index]
    tmp_node <- variant_node[index]
    
    tmp_id_unlisted <- unlist(tmp_id)
    tmp_i_unlisted <- match(tmp_id_unlisted, featureID(features))

    i_E <- which(type(features)[tmp_i_unlisted] == "E")

    if (length(i_E) > 0) {
    
        tmp_i_unlisted[i_E] <- match(tmp_node[togroup(tmp_id)][i_E],
            feature2name(features))

    }
    
    tmp_id_unlisted <- featureID(features)[tmp_i_unlisted]
    tmp_id <- relist(tmp_id_unlisted, tmp_id)

    variant_rep_id[index] <- tmp_id
    variant_rep_id <- IntegerList(variant_rep_id)
    
    ## exclude variants due to open events or ambiguous features
    event_dup <- tapply(unlist(variant_rep_id),
        variant_info$eventID[togroup(variant_rep_id)],
        function(x) { any(duplicated(x)) })
    i <- which(!variant_informative |
        variant_info$eventID %in% names(which(event_dup)))
    variant_rep_id[i] <- vector("list", length(i))    

    return(variant_rep_id)
    
}

collapseRows <- function(x, list_i, fun = sum)
{

    y <- matrix(NA, nrow = length(list_i), ncol = ncol(x))
    
    j <- which(elementLengths(list_i) == 1)

    if (length(j) > 0) {

        i <- unlist(list_i[j])
        y[j, ] <- x[i, ]
        
    }

    j <- which(elementLengths(list_i) > 1)

    if (length(j) > 0) {

        i <- unlist(list_i[j])
        f <- togroup(list_i[j])
        y[j, ] <- apply(x[i, , drop = FALSE], 2,
            function(z) { tapply(z, f, fun) })

    }

    colnames(y) <- colnames(x)
    
    return(y)
    
}

getVariantFreq <- function(SE)
{
    
    U_start <- assay(SE, "countsVariant5p")
    V_start <- assay(SE, "countsTotal5p")
    U_end <- assay(SE, "countsVariant3p")
    V_end <- assay(SE, "countsTotal3p")

    variants <- rowData(SE)

    informative_start <- elementLengths(featureID5p(variants)) > 0
    informative_end <- elementLengths(featureID3p(variants)) > 0

    i_both <- which(informative_start & informative_end)
    i_start <- which(informative_start & !informative_end)
    i_end <- which(!informative_start & informative_end)

    U <- matrix(NA, nrow = nrow(SE), ncol = ncol(SE))
    U[i_both, ] <- U_start[i_both, ] + U_end[i_both, ]
    U[i_start, ] <- U_start[i_start, ]
    U[i_end, ] <- U_end[i_end, ]

    V <- matrix(NA, nrow = nrow(SE), ncol = ncol(SE))
    V[i_both, ] <- V_start[i_both, ] + V_end[i_both, ]
    V[i_start, ] <- V_start[i_start, ]
    V[i_end, ] <- V_end[i_end, ]

    X <- U/V
    X[is.na(X)] <- NA

    assay(SE, "variantFreq") <- X
        
    return(SE)
    
}
