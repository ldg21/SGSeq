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
##' @param read_length Read length (required for use with \code{alpha})
##' @param frag_length Fragment length for paired-end data (required
##'   for use with \code{alpha})
##' @param lib_size Number of aligned fragments (required for use with
##'   \code{alpha})
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
##' @param include_counts Logical indicating whether counts of
##'   compatible fragments should be included in elementMetadata column
##'   \dQuote{N}
##' @param retain_coverage Logical indicating whether coverages for each
##'   exon should be retained as an \code{IntegerList} in elementMetadata
##'   column \dQuote{coverage}. This allows filtering of features
##'   using more stringent criteria after the initial prediction.
##' @param junctions_only Logical indicating whether predictions
##'   should be limited to identification of splice junctions only
##' @param cores Number of cores available for parallel processing
##' @return A \code{TxFeatures} object
##' @keywords internal
##' @author Leonard Goldstein

predictTxFeaturesPerSample <- function(file_bam, which = NULL,
    paired_end, read_length = NULL, frag_length = NULL, lib_size = NULL,
    min_junction_count = NULL, alpha = NULL, psi, beta, gamma,
    include_counts = TRUE, retain_coverage = FALSE,
    junctions_only = FALSE, cores = 1)
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

    si <- seqinfo(BamFile(file_bam))

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
        predictTxFeaturesRanges,
        file_bam = file_bam,
        paired_end = paired_end,
        min_junction_count = min_junction_count,
        psi = psi,
        beta = beta,
        gamma = gamma,
        include_counts = include_counts,
        retain_coverage = retain_coverage,
        junctions_only = junctions_only,
        mc.preschedule = FALSE,
        mc.cores = cores)
    
    list_features <- list_features[!sapply(list_features, is.null)]

    if (length(list_features) == 0) { return(TxFeatures()) }
    
    features <- do.call(c, setNames(list_features, NULL))
    features <- sort(features)
    features <- TxFeatures(features)
    
    return(features)
    
}

##' @title Identification of splice junctions and exons for a given
##'   chromosome and strand
##' @inheritParams predictTxFeaturesPerSample
##' @return \code{GRanges} of predicted features
##' @keywords internal
##' @author Leonard Goldstein

predictTxFeaturesRanges <- function(file_bam, paired_end, which,
    min_junction_count, psi, beta, gamma, include_counts, retain_coverage,
    junctions_only)
{

    seqlevel <- as.character(seqnames(which))
    strand <- as.character(strand(which))
    
    si <- seqinfo(BamFile(file_bam))
    
    gap <- readGap(file_bam, paired_end, which)
    gap <- gap[strand(gap) %in% c(strand, "*")]
    if (length(gap) == 0) { return() }

    frag_exonic <- ranges(grglist(gap, drop.D.ranges = TRUE))
    frag_intron <- ranges(junctions(gap))
    
    ir <- predictSpliced(frag_exonic, frag_intron, min_junction_count,
        psi, beta, gamma, include_counts, retain_coverage, junctions_only)
    if (is.null(ir)) { return() }
        
    gr <- constructGRangesFromRanges(ir, seqlevel, strand, si)
    gr <- gr[gr %over% which]
    
    return(gr)
    
}

##' @title \code{IRanges}-based identification of splice junctions and exons
##' @inheritParams predictTxFeaturesPerSample
##' @param frag_exonic \code{IRangesList} with exonic regions from alignments
##' @param frag_intron \code{IRangesList} with introns implied by spliced
##'   alignments
##' @return \code{IRanges} with predicted features
##' @keywords internal
##' @author Leonard Goldstein

predictSpliced <- function(frag_exonic, frag_intron, min_junction_count,
    psi, beta, gamma, include_counts, retain_coverage, junctions_only)
{
    
    junctions <- predictJunctions(frag_exonic, frag_intron, min_junction_count,
        psi, include_counts, retain_coverage)
    if (is.null(junctions)) { return() }

    features <- junctions

    if (!junctions_only) {
    
        lower <- max(min_junction_count * beta, 1)
        islands <- slice(coverage(unlist(frag_exonic)), lower,
            rangesOnly = TRUE)

        candidates <- predictCandidatesInternal(junctions, islands)
        exons_I <- predictExonsInternal(candidates, frag_exonic,
            frag_intron, beta, include_counts, retain_coverage)
        if (!is.null(exons_I)) { features <- c(features, exons_I) }
        
        candidates <- predictCandidatesTerminal(junctions, islands, "exon_L")
        exons_L <- predictExonsTerminal(candidates, frag_exonic, frag_intron,
            gamma, "exon_L", include_counts, retain_coverage)
        if (!is.null(exons_L)) { features <- c(features, exons_L) }

        candidates <- predictCandidatesTerminal(junctions, islands, "exon_R")
        exons_R <- predictExonsTerminal(candidates, frag_exonic, frag_intron,
            gamma, "exon_R", include_counts, retain_coverage)
        if (!is.null(exons_R)) { features <- c(features, exons_R) }

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
##' @return \code{IRanges} of splice junctions with elementMetadata
##'   column \dQuote{type} and optionally \dQuote{N} (for
##'   \code{include_counts = TRUE}), \dQuote{N_splicesite} (for
##'   \code{retain_coverage = TRUE})
##' @keywords internal
##' @author Leonard Goldstein

predictJunctions <- function(frag_exonic, frag_intron, min_junction_count,
    psi, include_counts, retain_coverage)
{

    ## extract all splice junctions
    
    junctions <- unique(unlist(frag_intron)) + 1
    if (length(junctions) == 0) { return() }
    mcols(junctions) <- DataFrame("type" = rep("J", length(junctions)))
    mcols(junctions)$N <- junctionCompatible(junctions, frag_intron)

    ## consider splice junctions with counts at least min_junction_count

    junctions <- junctions[which(mcols(junctions)$N >= min_junction_count)]
    if (length(junctions) == 0) { return() }
        
    ## consider splice junctions and splicesites with
    ## counts >= psi * max(splicesite counts)
    
    ## Note left/right (L/R) nomenclature: for the LHS splice site,
    ## the spliced boundary is situated on the right, for the RHS
    ## splice site, the spliced boundary is situated on the left

    splicesite_N <- splicesiteCounts(junctions, frag_exonic, frag_intron,
        "junction")
    index <- which(mcols(junctions)$N >= psi * max(splicesite_N))
    if (length(index) == 0) { return() }

    junctions <- junctions[index]

    if (!include_counts) {

        mcols(junctions)$N <- NULL

    }

    if (retain_coverage) {

        mcols(junctions)$N_splicesite <- splicesite_N[index]

    }
    
    junctions <- completeMcols(junctions, include_counts, retain_coverage)
    
    return(junctions)

}

##' Identify candidate internal exons based on previously identified
##' splice junctions and regions with minimal read coverage required
##' for internal exon prediction.
##' 
##' @title Identify candidate internal exons
##' @param junctions \code{IRanges} of splice junctions
##' @param islands \code{IRanges} of genomic regions with minimal read
##'   coverage required for internal exon prediction
##' @return \code{IRanges} of candidate internal exons 
##' @keywords internal
##' @author Leonard Goldstein

predictCandidatesInternal <- function(junctions, islands)
{

    ## for each island, identify overlapping splice junctions

    island_junction <- as.list(findOverlaps(islands, junctions))

    ## for each island, obtain all pairs of overlapping splice junctions

    island_junction_pairs <- mapply(expand.grid,
        island_junction, island_junction, SIMPLIFY = FALSE)
    junction_pairs <- unique(do.call(rbind, island_junction_pairs))

    ## retain pairs of splice junctions that are consistent with
    ## flanking an internal exon

    candidate_start <- end(junctions)[junction_pairs[, 1]]
    candidate_end <- start(junctions)[junction_pairs[, 2]]
    i <- which(candidate_start <= candidate_end)
    candidates <- unique(IRanges(candidate_start[i], candidate_end[i]))

    return(candidates)

}

##' Identify internal exons based on candidate internal exons and
##' compatible read coverage.
##' 
##' @title Identify internal exons
##' @inheritParams predictTxFeaturesPerSample
##' @inheritParams predictSpliced
##' @param candidates \code{IRanges} of candidate internal exons
##' @param relCov Minimum relative coverage required for exon prediction
##' @return \code{IRanges} of internal exons with elementMetadata column
##'   \dQuote{type} and optionally \dQuote{N} (for
##'   \code{include_counts = TRUE}), \dQuote{N_splicesite},
##'   \dQuote{coverage} (for \code{retain_coverage = TRUE})
##' @keywords internal
##' @author Leonard Goldstein

predictExonsInternal <- function(candidates, frag_exonic, frag_intron, relCov,
    include_counts, retain_coverage)
{

    if (length(candidates) == 0) { return() }
    
    candidate_index <- exonCompatible(candidates, TRUE, TRUE,
        frag_exonic, frag_intron, FALSE)
    candidate_coverage <- exonCoverage(candidates, candidate_index,
        frag_exonic)
    candidate_N_splicesite <- splicesiteCounts(candidates, frag_exonic,
        frag_intron, "exon")
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

    exons <- completeMcols(exons, include_counts, retain_coverage)
    
    return(exons)
    
}

##' Identify candidate terminal exons based on previously identified
##' splice junctions and regions with minimal read coverage required
##' for terminal exon prediction.
##' 
##' @title Identify candidate terminal exons
##' @inheritParams predictCandidatesInternal
##' @param type Character string indicating whether terminal exons
##'   should be identified to the left (\dQuote{exon_L}) or right
##'   (\dQuote{exon_R}) of provided splice junctions 
##' @return \code{IRanges} of candidate terminal exons
##' @keywords internal
##' @author Leonard Goldstein

predictCandidatesTerminal <- function(junctions, islands,
    type = c("exon_L", "exon_R"))
{

    type <- match.arg(type)
    
    splicesite <- unique(flank(junctions, -1,
        start = switch(type, "exon_L" = TRUE, "exon_R" = FALSE)))
    hits <- findOverlaps(splicesite, islands)
    spliced_boundary <- splicesite[queryHits(hits)]
    island <- islands[subjectHits(hits)]

    if (type == "exon_L") {

        candidates <- IRanges(start(island), end(spliced_boundary))

    }
    if (type == "exon_R") {

        candidates <- IRanges(start(spliced_boundary), end(island))

    }
    
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
##' @return \code{IRanges} of terminal exons with elementMetadata column
##'   \dQuote{type} and optionally \dQuote{N} (for
##'   \code{include_counts = TRUE}), \dQuote{N_splicesite},
##'   \dQuote{coverage} (for \code{retain_coverage = TRUE})
##' @keywords internal
##' @author Leonard Goldstein

predictExonsTerminal <- function(candidates, frag_exonic, frag_intron, relCov,
    type = c("exon_L", "exon_R"), include_counts, retain_coverage)
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
        frag_exonic, frag_intron, "spliced")
    ranges <- (coverage >= relCov * N_splicesite)
    el <- elementLengths(ranges)
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
            
        mcols(exons)$N_splicesite <- as(N_splicesite, "CompressedIntegerList")
        mcols(exons)$coverage <- coverage[setNames(split(ir, seq_along(ir)),
            NULL)]
            
    }
        
    exons <- completeMcols(exons, include_counts, retain_coverage)

    return(exons)
    
}
