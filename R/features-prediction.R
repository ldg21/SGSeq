##' Splice junctions and exons are predicted from genomic RNA-seq read
##' alignments in BAM format.
##'
##' For spliced alignments, the direction of transcription is inferred from
##' the XS tag in the BAM file and used to assign strand information to
##' the read, or fragment for paired-end data.
##'
##' Feature prediction is performed in two steps. First, splice junctions
##' are identified from spliced alignments. Second, exons
##' are identified based on regions that are flanked by splice
##' junctions and show sufficient coverage with compatible reads.
##'
##' Splice junctions implied by read alignments are filtered based on
##' fragment count and splice frequency. The splice frequency at the
##' splice donor (acceptor) is defined as x_J/x_D (x_J/x_A), where
##' x_J is the number of fragments containing the splice junction, and
##' x_D (x_A) is the number of fragments overlapping the exon/intron
##' (intron/exon) boundary. Fragments overlapping the spliced boundary
##' can be either spliced or extend into the intron. To be included in
##' predicted features, splice junctions must have fragment count at
##' least \code{min_junction_count} or FPKM at least \code{alpha}, and
##' splice frequency at both donor and acceptor at least \code{psi}.
##'
##' Regions between any pair of identified splice junctions with sufficient
##' compatible read coverage are considered candidate internal exons.
##' Read coverage for a candidate exon is computed based on compatible
##' fragments, i.e. fragments with matching (or missing) strand information
##' and introns consistent with the exon under consideration.
##' Candidate exons are included in predicted features if the minimum
##' coverage is at least \code{beta} * number of junction-containing
##' fragments for either flanking junctions.
##'
##' Terminal exons are regions downstream or upstream of splice junctions
##' with compatible fragment coverage at least \code{gamma} * number of
##' junction-containing fragments.
##'
##' @title Identification of splice junctions and exons from BAM file
##' @param file_bam BAM file with genomic RNA-seq read alignments
##' @param which \code{GRanges} of genomic regions to be considered for
##'   feature prediction, passed to \code{ScanBamParam}
##' @param paired_end Logical, \code{TRUE} for paired-end data,
##'   \code{FALSE} for single-end data
##' @param read_length Read length required for use with \code{alpha}
##' @param frag_length Fragment length for paired-end data required
##'   for use with \code{alpha}
##' @param lib_size Number of aligned fragments required for use with
##'   \code{alpha}
##' @param min_junction_count Minimum fragment count required for a splice
##'   junction to be included. If specified, argument \code{alpha} is ignored.
##' @param alpha Minimum FPKM required for a splice junction to be
##'   included. Internally, FPKMs are converted to counts, requiring arguments
##'   \code{read_length}, \code{frag_length} and \code{lib_size}.
##'   \code{alpha} is ignored if argument \code{min_junction_count}
##'   is specified.
##' @param psi Minimum splice frequency required for a splice junction
##'   to be included
##' @param beta Minimum relative coverage required for an internal exon
##'   to be included
##' @param gamma Minimum relative coverage required for a terminal exon
##'   to be included
##' @param min_anchor Integer specifiying minimum anchor length
##' @param include_counts Logical indicating whether counts of
##'   compatible fragments should be included in metadata column
##'   \dQuote{N}
##' @param retain_coverage Logical indicating whether coverage for each
##'   exon should be retained as an \code{RleList} in metadata
##'   column \dQuote{coverage}. This allows filtering of features
##'   using more stringent criteria after the initial prediction.
##' @param junctions_only Logical indicating whether predictions
##'   should be limited to identification of splice junctions only
##' @param max_complexity Maximum allowed complexity. If a locus exceeds
##'   this threshold, it is skipped, resulting in a warning.
##'   Complexity is defined as the maximum number of unique predicted
##'   splice junctions overlapping a given position.
##'   High complexity regions are often due to spurious read alignments
##'   and can slow down processing. To disable this filter, set to \code{NA}.
##' @param sample_name Sample name used in messages
##' @param verbose If \code{TRUE}, generate messages indicating progress
##' @param cores Number of cores available for parallel processing
##' @return \code{TxFeatures} object
##' @keywords internal
##' @author Leonard Goldstein

predictTxFeaturesPerSample <- function(file_bam, which, paired_end,
    read_length, frag_length, lib_size, min_junction_count,
    alpha, psi, beta, gamma, min_anchor, include_counts, retain_coverage,
    junctions_only, max_complexity, sample_name, verbose, cores)
{

    if (is.null(min_junction_count) && is.null(alpha)) {

        stop("Need to provide min_junction_count or alpha")

    }

    if (is.null(min_junction_count) &&
        (is.null(read_length) || is.null(frag_length) || is.null(lib_size))) {

        stop("For use with alpha, need to provide read_length,
            frag_length and lib_size")

    }

    if (!is.null(which) && !is(which, "GRanges")) {

        stop("argument which must be NULL or a GRanges object")

    }

    if (is.null(min_junction_count)) {

        min_junction_count <- convertFpkmToCount(alpha, paired_end,
            read_length, frag_length, lib_size)

    }

    if (is(file_bam, "BamFile")) {

        si <- seqinfo(file_bam)

    } else {

        si <- seqinfo(BamFile(file_bam))

    }

    if (is.null(which)) {

        sl <- rep(seqlevels(si), 2)
        st <- rep(c("+", "-"), rep(length(si), 2))
        which <- GRanges(sl, IRanges(1, seqlengths(si)[sl]), st)

    } else {

        which <- expandUnstrandedRanges(which)

    }

    list_which <- split(which, seq_along(which))

    list_features <- mclapply(
        list_which,
        predictTxFeaturesPerStrand,
        file_bam = file_bam,
        paired_end = paired_end,
        min_junction_count = min_junction_count,
        psi = psi,
        beta = beta,
        gamma = gamma,
        min_anchor = min_anchor,
        include_counts = include_counts,
        retain_coverage = retain_coverage,
        junctions_only = junctions_only,
        max_complexity = max_complexity,
        sample_name = sample_name,
        verbose = verbose,
        mc.preschedule = setPreschedule(cores),
        mc.cores = cores)

    checkApplyResultsForErrors(
        list_features,
        "predictTxFeaturesPerStrand",
        gr2co(unlist(range(list_which))),
        "try-error")

    list_features <- list_features[!vapply(list_features, is.null, logical(1))]

    if (length(list_features) == 0) {

        features <- TxFeatures()

    } else {

        features <- do.call(c, setNames(list_features, NULL))
        features <- sort(features)
        features <- TxFeatures(features)

    }

    generateCompleteMessage(sample_name)

    return(features)

}

##' Identification of splice junctions and exons for a given chromosome
##' and strand.
##'
##' @title Identification of splice junctions and exons for a given
##'   chromosome and strand
##' @inheritParams predictTxFeaturesPerSample
##' @return \code{GRanges} of predicted features
##' @keywords internal
##' @author Leonard Goldstein

predictTxFeaturesPerStrand <- function(file_bam, paired_end, which,
    min_junction_count, psi, beta, gamma, min_anchor, include_counts,
    retain_coverage, junctions_only, max_complexity, sample_name, verbose)
{

    seqlevel <- as.character(seqnames(which))
    strand <- as.character(strand(which))

    if (is(file_bam, "BamFile")) {

        si <- seqinfo(file_bam)

    } else {

        si <- seqinfo(BamFile(file_bam))

    }

    gap <- readGap(file_bam, paired_end, which, sample_name, verbose)
    frag_exonic <- gap$frag_exonic
    frag_intron <- gap$frag_intron

    if (length(frag_exonic) == 0) {

        gr <- NULL

    } else {

        ir <- predictSpliced(frag_exonic, frag_intron, min_junction_count,
            psi, beta, gamma, min_anchor, include_counts, retain_coverage,
            junctions_only, max_complexity, sample_name, seqlevel, strand)

        if (is.null(ir)) {

            gr <- NULL

        } else {

            gr <- constructGRangesFromRanges(ir, seqlevel, strand, si)

        }

    }

    if (verbose) generateCompleteMessage(paste(sample_name, gr2co(which)))

    return(gr)

}

##' Ranges-based identification of splice junctions and exons.
##'
##' @title Ranges-based identification of splice junctions and exons
##' @inheritParams predictTxFeaturesPerSample
##' @param frag_exonic \code{IRangesList} with exonic regions from alignments
##' @param frag_intron \code{IRangesList} with introns implied by spliced
##'   alignments
##' @param min_anchor Integer specifiying minimum anchor length
##' @param seqlevel \code{seqlevel} to be processed
##' @param strand \code{strand} to be processed
##' @return \code{IRanges} with predicted features
##' @keywords internal
##' @author Leonard Goldstein

predictSpliced <- function(frag_exonic, frag_intron, min_junction_count,
    psi, beta, gamma, min_anchor, include_counts, retain_coverage,
    junctions_only, max_complexity, sample_name, seqlevel, strand)
{

    junctions <- predictJunctions(frag_exonic, frag_intron,
        min_junction_count, psi, min_anchor, retain_coverage)
    if (is.null(junctions)) { return() }

    if (!junctions_only) {

        lower <- max(min_junction_count * beta, 1)
        frag_coverage <- coverage(unlist(frag_exonic))
        islands <- slice(frag_coverage, lower, rangesOnly = TRUE)

        ## skip problematic regions

        if (!is.na(max_complexity)) {

            ir <- as(slice(coverage(junctions), max_complexity), "IRanges")

            if (length(ir) > 0) {

                junctions_stripped <- junctions
                mcols(junctions_stripped) <- NULL
                loci <- reduce(c(junctions_stripped, islands))
                excl <- loci[loci %over% ir]

                junctions <- junctions[!junctions %over% excl]
                islands <- islands[!islands %over% excl]

                excl_str <- co2str(seqlevel, start(excl), end(excl), strand)

                generateWarningMessage(
                    "predictSpliced",
                    sample_name,
                    paste("skipping", excl_str))

                if (length(junctions) == 0) { return() }

            }

        }

        features <- junctions

        splicesites_L <- extractSplicesitesFromJunctions(junctions, "L")
        splicesites_R <- extractSplicesitesFromJunctions(junctions, "R")
        splicesites <- c(splicesites_L, splicesites_R)

        candidates <- predictCandidatesInternal(islands, splicesites,
            frag_coverage, beta)
        exons_I <- predictExonsInternal(candidates, frag_exonic,
            frag_intron, beta, min_anchor, include_counts, retain_coverage)
        if (!is.null(exons_I)) { features <- c(features, exons_I) }

        candidates <- predictCandidatesTerminal(islands, splicesites, "exon_L")
        exons_L <- predictExonsTerminal(candidates, frag_exonic, frag_intron,
            gamma, min_anchor, "exon_L", include_counts, retain_coverage)
        if (!is.null(exons_L)) { features <- c(features, exons_L) }

        candidates <- predictCandidatesTerminal(islands, splicesites, "exon_R")
        exons_R <- predictExonsTerminal(candidates, frag_exonic, frag_intron,
            gamma, min_anchor, "exon_R", include_counts, retain_coverage)
        if (!is.null(exons_R)) { features <- c(features, exons_R) }

    } else {

        features <- junctions

    }

    if (!include_counts) {

        mcols(features)$N <- NULL

    }

    features <- sort(features)

    return(features)

}

constructGRangesFromRanges <- function(x, seqname, strand, seqinfo)
{

    x_mcols <- mcols(x)
    mcols(x) <- NULL

    if (strand == "+") {

        x_mcols_type <- as.character(x_mcols$type)
        x_mcols_type <- sub("exon_L", "F", x_mcols_type, fixed = TRUE)
        x_mcols_type <- sub("exon_R", "L", x_mcols_type, fixed = TRUE)
        x_mcols$type <- factor(x_mcols_type,
            levels = c("J", "I", "F", "L", "U"))

    }
    if (strand == "-") {

        x_mcols_type <- as.character(x_mcols$type)
        x_mcols_type <- sub("exon_L", "L", x_mcols_type, fixed = TRUE)
        x_mcols_type <- sub("exon_R", "F", x_mcols_type, fixed = TRUE)
        x_mcols$type <- factor(x_mcols_type,
            levels = c("J", "I", "F", "L", "U"))

        if ("N_splicesite" %in% names(x_mcols)) {

            x_mcols$N_splicesite <- endoapply(x_mcols$N_splicesite, rev)

        }
        if ("coverage" %in% names(x_mcols)) {

            x_mcols$coverage <- endoapply(x_mcols$coverage, rev)

        }

    }

    gr <- GRanges(rep(seqname, length(x)), x, rep(strand, length(x)),
        x_mcols, seqinfo = seqinfo)

    return(gr)

}

##' Identify splice junctions from genomic RNA-seq read alignments.
##'
##' @title Identify splice junctions
##' @inheritParams predictTxFeaturesPerSample
##' @inheritParams predictSpliced
##' @return \code{IRanges} of splice junctions with metadata
##'   columns \dQuote{type} and \dQuote{N}, and optionally
##'   \dQuote{N_splicesite} for \code{retain_coverage = TRUE}
##' @keywords internal
##' @author Leonard Goldstein

predictJunctions <- function(frag_exonic, frag_intron,
    min_junction_count, psi, min_anchor, retain_coverage)
{

    ## extract all splice junctions

    junctions <- unique(unlist(frag_intron)) + 1
    if (length(junctions) == 0) { return() }
    mcols(junctions) <- DataFrame("type" = rep("J", length(junctions)))

    ## consider splice junctions with counts at least min_junction_count

    mcols(junctions)$N <- junctionCompatible(junctions, frag_exonic,
        frag_intron, min_anchor)
    junctions <- junctions[which(mcols(junctions)$N >= min_junction_count)]
    if (length(junctions) == 0) { return() }

    ## consider splice junctions and splicesites with
    ## counts >= psi * max(splicesite counts)

    ## Note left/right (L/R) nomenclature: for the LHS splice site,
    ## the spliced boundary is situated on the right, for the RHS
    ## splice site, the spliced boundary is situated on the left

    if (psi > 0 || retain_coverage) {

        mcols(junctions)$N_splicesite <- splicesiteCounts(junctions,
            frag_exonic, frag_intron, min_anchor, "junction", "all")
        junctions <- junctions[which(mcols(junctions)$N >=
            psi * max(mcols(junctions)$N_splicesite))]
        if (length(junctions) == 0) { return() }

    }

    if (!retain_coverage) {

        mcols(junctions)$N_splicesite <- NULL

    }

    junctions <- completeMcols(junctions, retain_coverage)

    return(junctions)

}

##' Identify candidate internal exons based on previously identified
##' splice sites and regions with sufficient read coverage.
##'
##' @title Identify candidate internal exons
##' @param islands \code{IRanges} of genomic regions with minimal read
##'   coverage required for internal exon prediction
##' @param splicesites \code{IRanges} of splice sites with metadata
##'   columns \dQuote{type} and \dQuote{N}
##' @param frag_coverage \code{Rle} object with fragment coverage
##' @param relCov Minimum relative coverage required for exon prediction
##' @return \code{IRanges} of candidate internal exons
##' @keywords internal
##' @author Leonard Goldstein

predictCandidatesInternal <- function(islands, splicesites, frag_coverage,
    relCov)
{

    ## for each island, identify overlapping splice sites

    island_splicesite <- as.list(findOverlaps(islands, splicesites))

    ## for each island, obtain all pairs of overlapping splice sites

    island_splicesite_pairs <- mapply(expand.grid,
        island_splicesite, island_splicesite, SIMPLIFY = FALSE)
    splicesite_pairs <- unique(do.call(rbindDfsWithoutRowNames,
        island_splicesite_pairs))

    ## retain pairs of splice sites that are consistent with
    ## flanking an internal exon

    N_1 <- mcols(splicesites)$N[splicesite_pairs[, 1]]
    N_2 <- mcols(splicesites)$N[splicesite_pairs[, 2]]
    type_1 <- mcols(splicesites)$type[splicesite_pairs[, 1]]
    type_2 <- mcols(splicesites)$type[splicesite_pairs[, 2]]
    pos_1 <- start(splicesites)[splicesite_pairs[, 1]]
    pos_2 <- start(splicesites)[splicesite_pairs[, 2]]
    i <- which(type_1 == "R" & type_2 == "L" & pos_1 <= pos_2)
    candidates <- IRanges(pos_1[i], pos_2[i])
    mcols(candidates) <- DataFrame(N = IntegerList(
        mapply(c, N_1[i], N_2[i], SIMPLIFY = FALSE)))

    if (length(candidates) > 0) {

        ## retain candidate internal exons with sufficient read coverage

        candidates_frag_coverage <- split(frag_coverage[candidates],
            togroup0(candidates))
        i <- which(min(candidates_frag_coverage) >=
            relCov * min(mcols(candidates)$N))
        candidates <- candidates[i]

    }

    return(candidates)

}

##' Identify internal exons based on candidate internal exons and
##' compatible read coverage.
##'
##' @title Identify internal exons
##' @inheritParams predictTxFeaturesPerSample
##' @inheritParams predictSpliced
##' @inheritParams predictCandidatesInternal
##' @param candidates \code{IRanges} of candidate internal exons
##' @return \code{IRanges} of internal exons with metadata column
##'   \dQuote{type} and optionally \dQuote{N} for
##'   \code{include_counts = TRUE}, \dQuote{N_splicesite},
##'   \dQuote{coverage} for \code{retain_coverage = TRUE}
##' @keywords internal
##' @author Leonard Goldstein

predictExonsInternal <- function(candidates, frag_exonic, frag_intron, relCov,
    min_anchor, include_counts, retain_coverage)
{

    if (length(candidates) == 0) { return() }

    candidate_index <- exonCompatible(candidates, TRUE, TRUE,
        frag_exonic, frag_intron, FALSE)
    candidate_coverage <- exonCoverage(candidates, candidate_index,
        frag_exonic)
    candidate_N_splicesite <- splicesiteCounts(candidates, frag_exonic,
        frag_intron, min_anchor, "exon", "spliced")

    index <- which(min(candidate_coverage) >=
        relCov * min(candidate_N_splicesite))
    if (length(index) == 0) { return() }

    exons <- candidates[index]
    mcols(exons) <- DataFrame("type" = rep("I", length(exons)))

    if (include_counts) {

        mcols(exons)$N <- exonCompatible(exons, TRUE, TRUE, frag_exonic,
            frag_intron)

    }

    if (retain_coverage) {

        mcols(exons)$N_splicesite <- candidate_N_splicesite[index]
        mcols(exons)$coverage <- candidate_coverage[index]

    }

    exons <- completeMcols(exons, retain_coverage)

    return(exons)

}

##' Identify candidate terminal exons based on previously identified
##' splice sites and regions with sufficient read coverage.
##'
##' @title Identify candidate terminal exons
##' @inheritParams predictCandidatesInternal
##' @param type Character string indicating whether terminal exons
##'   should be identified to the left (\dQuote{exon_L}) or right
##'   (\dQuote{exon_R}) of provided splice sites
##' @return \code{IRanges} of candidate terminal exons
##' @keywords internal
##' @author Leonard Goldstein

predictCandidatesTerminal <- function(islands, splicesites,
    type = c("exon_L", "exon_R"))
{

    type <- match.arg(type)

    splicesites <- splicesites[mcols(splicesites)$type ==
        switch(type, "exon_L" = "L", "exon_R" = "R")]
    hits <- findOverlaps(splicesites, islands)
    spliced_boundary <- splicesites[queryHits(hits)]
    island <- islands[subjectHits(hits)]

    if (type == "exon_L") {

        candidates <- IRanges(start(island), end(spliced_boundary))

    }
    if (type == "exon_R") {

        candidates <- IRanges(start(spliced_boundary), end(island))

    }

    mcols(candidates) <- DataFrame(N = mcols(splicesites)$N[queryHits(hits)])

    return(candidates)

}

##' Identify terminal exons based on candidate terminal exons and
##' compatible read coverage.
##'
##' @title Identify terminal exons
##' @inheritParams predictTxFeaturesPerSample
##' @inheritParams predictSpliced
##' @inheritParams predictExonsInternal
##' @inheritParams predictCandidatesTerminal
##' @param relCov Minimum relative coverage required for exon prediction
##' @return \code{IRanges} of terminal exons with metadata column
##'   \dQuote{type} and optionally \dQuote{N} for
##'   \code{include_counts = TRUE}, \dQuote{N_splicesite},
##'   \dQuote{coverage} for \code{retain_coverage = TRUE}
##' @keywords internal
##' @author Leonard Goldstein

predictExonsTerminal <- function(candidates, frag_exonic, frag_intron, relCov,
    min_anchor, type = c("exon_L", "exon_R"), include_counts, retain_coverage)
{

    type <- match.arg(type)

    if (length(candidates) == 0) { return() }

    spliceL <- switch(type, "exon_L" = FALSE, "exon_R" = TRUE)
    spliceR <- switch(type, "exon_L" = TRUE, "exon_R" = FALSE)

    index <- exonCompatible(candidates, spliceL, spliceR,
        frag_exonic, frag_intron, FALSE)
    coverage <- exonCoverage(candidates, index, frag_exonic)

    splicesite <- flank(candidates, -1,
        start = switch(type, "exon_L" = FALSE, "exon_R" = TRUE))
    N_splicesite <- splicesiteOverlap(splicesite,
        switch(type, "exon_L" = "R", "exon_R" = "L"),
        frag_exonic, frag_intron, min_anchor, "spliced")

    ranges <- (coverage >= relCov * N_splicesite)
    el <- elementNROWS(ranges)
    rl <- runLength(ranges)

    if (type == "exon_L") {

        ir <- IRanges(end = el, width = plast(rl))

    }
    if (type == "exon_R") {

        ir <- IRanges(start = 1, width = pfirst(rl))

    }

    if (length(ir) ==    0) { return() }

    exons <- shift(ir, start(candidates) - 1)
    mcols(exons) <- DataFrame(type = rep(type, length(exons)))

    if (include_counts) {

        mcols(exons)$N <- exonCompatible(exons, spliceL, spliceR,
            frag_exonic, frag_intron)

    }

    if (retain_coverage) {

        mcols(exons)$N_splicesite <- as(N_splicesite, "IntegerList")
        mcols(exons)$coverage <- coverage[
            setNames(split(ir, seq_along(ir)), NULL)]

    }

    exons <- completeMcols(exons, retain_coverage)

    return(exons)

}

extractSplicesitesFromJunctions <- function(junctions, type = c("L", "R"))
{

    type <- match.arg(type)
    S <- flank(junctions, -1, switch(type, "L" = TRUE, "R" = FALSE))
    S_pos <- as.character(start(S))
    pos_N <- tapply(mcols(junctions)$N, S_pos, sum)
    pos_N <- setNames(as.integer(pos_N), names(pos_N))
    i <- which(!duplicated(S_pos))
    S <- S[i]
    S_pos <- S_pos[i]
    mcols(S) <- DataFrame(type = rep(type, length(S)),
        N = pos_N[match(S_pos, names(pos_N))])
    return(S)

}
