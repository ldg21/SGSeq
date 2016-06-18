feature2name <- function(features, collapse_terminal = FALSE)
{

    if (is(features, "Features")) {

        features_type <- type(features)

    } else {

        features_type <- mcols(features)$type

    }

    name <- rep(NA, length(features))

    if (collapse_terminal) {

        i <- which(features_type %in% c("J", "I", "U", "E"))
        name[i] <- paste0(features_type[i], ":", gr2co(features[i]))
        i <- which(features_type %in% c("F", "L"))
        start <- ifelse(features_type[i] == "F", FALSE, TRUE)
        name[i] <- paste0(features_type[i], ":",
            gr2co(flank(features[i], -1, start)))

    } else {

        i <- which(features_type %in% c("J", "I", "F", "L", "U", "E"))
        name[i] <- paste0(features_type[i], ":", gr2co(features[i]))

    }

    i <- which(features_type %in% c("D", "A"))
    name[i] <- paste0(features_type[i], ":", gr2pos(features[i]))

    return(name)

}

co2str <- function(seqlevel, start, end, strand)
{

    paste0(seqlevel, ":", start, "-", end, ":", strand)

}

gr2co <- function(x)
{

    if (length(x) == 0) {

        return()

    } else {

        co2str(seqnames(x), start(x), end(x), strand(x))

    }

}

co2gr <- function(co)
{

    x <- strsplit(co, split = ":", fixed = TRUE)
    r <- strsplit(sapply(x, "[", 2), split = "-", fixed = TRUE)
    sn <- sapply(x, "[", 1)
    start <- as.integer(sapply(r, "[", 1))
    end <- as.integer(sapply(r, "[", 2))
    st <- sapply(x, "[", 3)
    GRanges(sn, IRanges(start, end), st)

}

pos2str <- function(seqlevel, position, strand)
{

    paste0(seqlevel, ":", position, ":", strand)

}

gr2pos <- function(x)
{

    if (length(x) == 0) {

        return()

    } else {

        pos2str(seqnames(x), start(x), strand(x))

    }

}

pos2gr <- function(x)
{

    x <- strsplit(x, split = ":", fixed = TRUE)
    sn <- sapply(x, "[", 1)
    pos <- as.integer(sapply(x, "[", 2))
    st <- sapply(x, "[", 3)
    GRanges(sn, IRanges(pos, pos), st)

}

readGap <- function(file, paired_end, which, sample_name, verbose)
{

    if (length(which) != 1) {

        stop("which argument must have length 1")

    }

    ## the following flags are set by functions
    ## readGAlignments and readGAlignmentPairs
    ## - isUnmappedQuery
    ## - isPaired
    ## - hasUnmappedMate

    flag <- scanBamFlag(isSecondaryAlignment = FALSE)
    param <- ScanBamParam(flag = flag, tag = "XS", which = which)

    if (paired_end) {

        gap <- suppressWarnings(readGAlignmentPairs(file = file,
            param = param))

        ## scanBam workaround start
        ## bamWhat(param) <- c("flag", "mrnm", "mpos")
        ## ga <- readGAlignments(file = file, use.names = TRUE, param = param)
        ## gap <- makeGAlignmentPairs(ga, use.names = TRUE, use.mcols = TRUE)
        ## names(gap) <- NULL
        ## scanBam workaround end

        gap <- propagateXS(gap)

    } else {

        gap <- suppressWarnings(readGAlignments(file = file, param = param))

    }

    gap <- filterGap(gap)

    mcols(gap)$strand <- XS2strand(mcols(gap)$XS)

    gap <- gap[mcols(gap)$strand %in% c(as.character(strand(which)), "*")]

    frag_exonic <- reduce(ranges(grglist(gap, drop.D.ranges = TRUE)))
    frag_intron <- ranges(junctions(gap))

    if (paired_end) {

        diff <- setdiff(frag_exonic, frag_intron)
        excl <- which(sum(width(frag_exonic)) > sum(width(diff)))

        if (length(excl) > 0) {

            if (verbose) {

                msg <- paste(
                    "filtered",
                    length(excl),
                    "inconsistent paired alignments in",
                    gr2co(which))

                generateWarningMessage(
                    "readGap",
                    sample_name,
                    msg)

            }

            frag_exonic <- frag_exonic[-excl]
            frag_intron <- frag_intron[-excl]

        }

    }

    gap <- list(frag_exonic = frag_exonic, frag_intron = frag_intron)

    return(gap)

}

propagateXS <- function(gap)
{

    first_xs <- mcols(first(gap))$XS
    last_xs <- mcols(last(gap))$XS
    xs <- first_xs
    xs[is.na(xs)] <- last_xs[is.na(xs)]
    mcols(gap)$XS <- xs
    return(gap)

}

XS2strand <- function(xs)
{

    s <- xs
    s[is.na(s)|s == "?"] <- "*"
    return(s)

}

filterGap <- function(gap)
{

    if (is(gap, "GAlignments")) {

        exclude <- filterGa(gap)

    }
    if (is(gap, "GAlignmentPairs")) {

        exclude <- filterGa(first(gap)) | filterGa(last(gap))

    }

    gap <- gap[!exclude]

    return(gap)

}

filterGa <- function(ga)
{

    grepl("(\\d+D\\d+N)|(\\d+N\\d+D)", cigar(ga))

}

dropMcols <- function(x)
{

    mcols(x) <- NULL
    return(x)

}

completeMcols <- function(x, retain_coverage)
{

    mcol <- c("type", "N")

    if (retain_coverage) {

        mcol <- c(mcol, "N_splicesite", "coverage")

    }

    for (m in setdiff(mcol, names(mcols(x)))) {

        if (m == "N") {

            mcols(x)[, m] <- NA_integer_

        } else if (m == "N_splicesite") {

            mcols(x)[, m] <- IntegerList(vector("list", length(x)))

        } else if (m == "coverage") {

            mcols(x)[, m] <- RleList(IntegerList(vector("list", length(x))))

        }

    }

    mcols(x) <- mcols(x)[, mcol, drop = FALSE]
    names(mcols(x)) <- mcol

    return(x)

}

getBamInfoPerSample <- function(file_bam, yieldSize, sample_name)
{

    if (is(file_bam, "BamFile")) {

        file_tmp <- file_bam

    } else {

        file_tmp <- BamFile(file_bam)

    }

    if (!is.null(yieldSize)) {

        yieldSize(file_tmp) <- yieldSize

    }

    flag <- scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE)
    what <- c("qname", "flag", "qwidth", "isize")
    param <- ScanBamParam(flag = flag, what = what, tag = "XS")
    bam <- scanBam(file = file_tmp, param = param)[[1]]

    XS <- !is.null(bam$tag$XS)
    paired_end <- any(bamFlagTest(bam$flag, "isPaired"))
    read_length <- median(bam$qwidth, na.rm = TRUE)

    if (paired_end) {

        isize <- bam$isize
        frag_length <- median(isize[which(isize > 0)], na.rm = TRUE)

    } else {

        frag_length <- NA_real_

    }

    x <- data.frame(
        XS = XS,
        paired_end = paired_end,
        read_length = read_length,
        frag_length = frag_length,
        stringsAsFactors = FALSE)

    if (is.null(yieldSize)) {

        x$lib_size <- length(unique(bam$qname))

    }

    generateCompleteMessage(sample_name)

    return(x)

}

expandUnstrandedRanges <- function(x)
{

    i <- which(strand(x) == "*")

    if (length(i) > 0) {

        additional <- x[i]
        strand(additional) <- "-"
        strand(x)[i] <- "+"
        x <- c(x, additional)

    }

    return(x)

}

uniqueFeatures <- function(features)
{

    i_duplicated <- vector()

    for (type in levels(type(features))) {

        i_type <- which(type(features) == type)
        i <- i_type[which(duplicated(features[i_type]))]
        i_duplicated <- c(i_duplicated, i)

    }

    if (length(i_duplicated) > 0) {

        features <- features[-i_duplicated]

    }

    return(features)

}

##' Export features to BED format. Splice sites are not included.
##'
##' @title Export to BED format
##' @param features \code{TxFeatures} or \code{SGFeatures} object
##' @param file Character string specifying output file
##' @return \code{NULL}
##' @examples
##' \dontrun{
##' exportFeatures(txf_pred, "txf.bed")
##' exportFeatures(sgf_pred, "sgf.bed")
##' }
##' NULL
##' @author Leonard Goldstein

exportFeatures <- function(features, file)
{

    if (!is(features, "Features")) {

        stop("features must be a TxFeatures or SGFeatures object")

    }

    features <- asGRanges(features)

    i_splicesite <- which(mcols(features)$type %in% c("D", "A"))

    if (length(i_splicesite) > 0) {

        features <- features[-i_splicesite]

    }

    i_junction <- which(mcols(features)$type == "J")
    color <- mcols(features)$color
    mcols(features) <- NULL

    bed <- split(features, seq_along(features))

    if (length(i_junction) > 0) {

        bed[i_junction] <- setdiff(
           split(features[i_junction], seq_along(i_junction)),
           split(features[i_junction] - 1, seq_along(i_junction)))

    }

    if (!is.null(color)) {

        itemRgb <- rgb(t(col2rgb(color)), maxColorValue = 255)
        mcols(bed)$itemRgb <- itemRgb

    }

    names(bed) <- feature2name(features)

    export(object = bed, con = file, format = "BED")

    return()

}

nextFrame <- function(f, w, prev = FALSE)
{

    if (is(f, "list") || is(f, "List")) {

        f_unlisted <- unlist(f)
        w_unlisted <- w[togroup0(f)]
        n_unlisted <- nextFrame(f_unlisted, w_unlisted, prev)
        n <- relist(n_unlisted, f)

    } else {

        if (prev) {

            n <- ifelse(f != -1, (f - w) %% 3, -1)

        } else {

            n <- ifelse(f != -1, (f + w) %% 3, -1)

        }

    }

    return(n)

}

splitCharacterList <- function(x, f)
{

    if (!is(f, "factor")) {

        stop("f must be a factor")

    }

    x_unlisted <- setNames(unlist(x), NULL)
    f_unlisted <- f[togroup0(x)]
    y <- CharacterList(split(x_unlisted, f_unlisted))
    y <- unique(y)

    return(y)

}

asGRanges <- function(from)
{

    granges(from, use.mcols = TRUE)

}

asGRangesList <- function(from)
{

    as(from, "GRangesList")

}

reorderFeatures <- function(x)
{

    x_names <- names(x)
    x_mc <- mcols(x)
    features <- unlist(x, use.names = FALSE)
    features_x <- togroup0(x)
    i_pos <- which(strand(features) == "+" | strand(features) == "*")
    i_neg <- which(strand(features) == "-")
    i_pos <- i_pos[order(features[i_pos])]
    i_neg <- i_neg[order(features[i_neg], decreasing = TRUE)]
    i_all <- c(i_pos, i_neg)
    x <- split(features[i_all], features_x[i_all])
    names(x) <- x_names
    mcols(x) <- x_mc

    return(x)

}

pfirst <- function(x, use_names = FALSE)
{

    unlist(phead(x, 1), use.names = use_names)

}

plast <- function(x, use_names = FALSE)
{

    unlist(ptail(x, 1), use.names = use_names)

}

rbindListOfDFs <- function(x, cores)
{

    k <- names(x[[1]])

    df <- vector("list", length(k))

    for (j in seq_along(k)) {

        df[[j]] <- do.call(c, mclapply(x, "[[", j, mc.cores = cores))

    }

    names(df) <- k

    df <- DataFrame(df, check.names = FALSE)

    return(df)

}

checkApplyResultsForErrors <- function(out, fun_name, items, error_class)
{

    failed <- sapply(out, is, error_class)

    if (any(failed)) {

        i <- which(failed)
        err <- makeErrorMessage(fun_name, items[i], out[i])
        err <- paste0("\n", err)
        stop(err, call. = FALSE)

    }

}

makeErrorMessage <- function(fun_name, items, msgs)
{

    msg <- paste0("Error in ", fun_name, " for ", items, ":", "\n", msgs)
    msg <- paste(msg, collapse = "\n")

    return(msg)

}

makeWarningMessage <- function(fun_name, items, msgs)
{

    msg <- paste0("Warning in ", fun_name, " for ", items, ":", "\n", msgs)
    msg <- paste(msg, collapse = "\n")

    return(msg)

}

makeCompleteMessage <- function(item) {

    paste(item, "complete.")

}

generateWarningMessage <- function(fun_name, item, msg)
{

    message(makeWarningMessage(fun_name, item, msg))

}

generateCompleteMessage <- function(item)
{

    message(makeCompleteMessage(item))

}

getCoverage <- function(sample_info, which, sizefactor, cores)
{

    if (!is(which, "GRanges") || length(which) > 1) {

        stop("which must be a GRanges object of length 1")

    }

    list_cov <- mcmapply(
        getCoveragePerSample,
        file_bam = sample_info$file_bam,
        paired_end = sample_info$paired_end,
        sample_name = sample_info$sample_name,
        sizefactor = sizefactor,
        MoreArgs = list(which = which),
        SIMPLIFY = FALSE,
        USE.NAMES = FALSE,
        mc.preschedule = setPreschedule(cores),
        mc.cores = cores)

    return(list_cov)

}

getCoveragePerSample <- function(file_bam, paired_end, sample_name,
    sizefactor, which)
{

    gap <- readGap(file_bam, paired_end, which, sample_name, FALSE)
    cov <- coverage(unlist(gap$frag_exonic), width = end(which))
    cov <- cov / sizefactor

    return(cov)

}

calculateSizeFactor <- function(sample_info)
{

    E <- rep(NA, nrow(sample_info))

    i_PE <- which(sample_info$paired_end)

    if (length(i_PE) > 0) {

        R_PE <- sample_info$read_length[i_PE]
        F_PE <- sample_info$frag_length[i_PE]
        I_PE <- F_PE - 2 * R_PE
        E[i_PE] <- F_PE - pmax(I_PE, 0)

    }

    i_SE <- which(!sample_info$paired_end)

    if (length(i_SE) > 0) {

        E[i_SE] <- sample_info$read_length[i_SE]

    }

    sizefactor <- sample_info$lib_size * E * 1e-9

    return(sizefactor)

}

checkIdenticalSummarizedExperiment <- function(target, current,
    check.colData = FALSE)
{

    checkTrue(is(target, "RangedSummarizedExperiment"))
    checkTrue(is(current, "RangedSummarizedExperiment"))

    slots <- c(
        "rowRanges",
        "colData",
        "assays",
        "NAMES",
        "elementMetadata",
        "metadata")

    checkIdentical(slots, slotNames(target))
    checkIdentical(slots, slotNames(current))

    slots <- slots[slots != "assays"]

    if (!check.colData) {

        slots <- slots[slots != "colData"]

    }

    for (s in slots) {

        checkIdentical(slot(target, s), slot(current, s))

    }

    assays_target <- names(assays(target))
    assays_current <- names(assays(current))

    checkIdentical(assays_target, assays_current)

    for (a in assays_target) {

        checkIdentical(assay(target, a), assay(current, a))

    }

    return(TRUE)

}

rbindDfsWithoutRowNames <- function(...)
{

    rbind(..., make.row.names = FALSE)

}

##' Import GFF file and generate a \code{GRangesList} of transcripts
##' suitable as input for functions \code{convertToTxFeatures} or
##' \code{predictVariantEffects}.
##'
##' @title Import transcripts from GFF file
##' @param file Character string specifying input GFF file
##' @param tag_tx GFF attribute tag for transcript identifier
##' @param tag_gene GFF attribute tag for gene identifier
##' @return \code{GRangesList} of exons grouped by transcipts with
##'   metadata columns txName, geneName, cdsStart, cdsEnd.
##' @examples
##' \dontrun{
##' tx <- importTranscripts(file)
##' }
##' NULL
##' @author Leonard Goldstein

importTranscripts <- function(file, tag_tx = "transcript_id",
    tag_gene = "gene_id")
{

  gff <- import(file)
  exons <- gff[mcols(gff)$type == "exon", ]
  df <- unique(data.frame(mcols(exons)[c(tag_tx, tag_gene)]))
  tx <- split(exons, mcols(exons)[[tag_tx]])
  cds <- gff[mcols(gff)$type == "CDS"]
  cds <- split(cds, mcols(cds)[[tag_tx]])
  cds <- unlist(range(cds))
  mcols(tx)$txName <- names(tx)
  mcols(tx)$geneName <- df[match(names(tx), df[, 1]), 2]
  mcols(tx)$cdsStart <- start(cds)[match(names(tx), names(cds))]
  mcols(tx)$cdsEnd <- end(cds)[match(names(tx), names(cds))]
  rownames(mcols(tx)) <- NULL

  return(tx)

}

convertToTranscripts <- function(txdb)
{

    tx <- exonsBy(txdb, "tx", use.names = TRUE)
    mcols(tx)$txName <- names(tx)
    df <- silent_select(txdb, names(tx), "GENEID", "TXNAME")
    mcols(tx)$geneName <- df$GENEID[match(names(tx), df$TXNAME)]
    cds <- unlist(range(cdsBy(txdb, "tx", use.names = TRUE)))
    cdsLeft(tx) <- start(cds)[match(names(tx), names(cds))]
    cdsRight(tx) <- end(cds)[match(names(tx), names(cds))]
    rownames(mcols(tx)) <- NULL

    return(tx)

}

checkTranscriptFormat <- function(x)
{

    if (!exonsOnSameChromAndStrand(x)) {

        msg <- "All exons for the same transcript must\n
            be on the same chromosome and strand"
        stop(msg, call. = FALSE)

    }

    mcol_type <- c(
        txName = "character",
        geneName = "character",
        cdsStart = "integer",
        cdsEnd = "integer")

    msg <- validMcols(x, mcol_type)

    if (!is.null(msg)) {

        stop(msg, call. = FALSE)

    }

    if (any(mcols(x)$cdsStart > mcols(x)$cdsEnd, na.rm = TRUE)) {

        msg <- "All coding transcripts must have cdsStart < cdsEnd"
        stop(msg, call. = FALSE)

    }

}

## Starting with IRanges 2.5.31, togroup() does not work on an arbitrary
## object anymore, only on a ManyToOneGrouping object.
## S4Vectors:::quick_togroup() is a replacement for the old togroup() that
## works on any object.

togroup0 <- S4Vectors:::quick_togroup

isOr <- function(object, class2)
{

    any(sapply(class2, is, object = object))

}

silent_select <- function(...)
{

    suppressMessages(select(...))

}
