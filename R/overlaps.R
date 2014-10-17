##' Identify fragments compatible with splice junctions.
##'
##' @title Compatible fragment counts for splice junctions
##' @inheritParams exonCompatible
##' @param junctions \code{IRanges} of splice junctions
##' @return Counts or list of indices of compatible fragments
##' @keywords internal
##' @author Leonard Goldstein

junctionCompatible <- function(junctions, frag_intron, counts = TRUE)
{

    introns <- junctions - 1
    hits <- findOverlapsRanges(introns, frag_intron, "equal")
    junction_index <- as.list(hits)
    
    if (counts) elementLengths(junction_index)
    else junction_index

}

##' Identify fragments compatible with exons.
##' 
##' @title Compatible fragment counts for exons
##' @param exons \code{IRanges} of exons
##' @param spliceL Logical vector indicating whether LHS boundary is spliced
##' @param spliceR Logical vector indicating whether RHS boundary is spliced
##' @param frag_exonic \code{IRangesList} of exonic regions, one entry
##'   per fragment
##' @param frag_intron \code{IRangesList} of introns, one entry per fragment
##' @param counts Logical indicating whether counts or indices of
##'   compatible fragments should be returned
##' @return Counts or list of indices of compatible fragments
##' @keywords internal
##' @author Leonard Goldstein

exonCompatible <- function(exons, spliceL, spliceR, frag_exonic,
    frag_intron, counts = TRUE)
{

    if (length(spliceL) == 1) { spliceL <- rep(spliceL, length(exons)) }
    if (length(spliceR) == 1) { spliceR <- rep(spliceR, length(exons)) }

    ## for each exon, identify alignments with overlapping exonic regions
    ## and no overlapping introns

    hits_exonic <- findOverlapsRanges(exons, frag_exonic)
    hits_introns <- findOverlapsRanges(exons, frag_intron)
    hits <- setdiff(hits_exonic, hits_introns)
    
    ## for exons with a LHS spliced boundary, exclude alignments with
    ## exonic regions overlapping the flanking LHS intronic position

    excl <- findOverlapsRanges(flank(exons, 1, TRUE), frag_exonic)
    excl <- excl[which(spliceL[queryHits(excl)])]
    hits <- setdiff(hits, excl)

    ## for exons with a RHS spliced boundary, exclude alignments with
    ## exonic regions overlapping the flanking RHS intronic position

    excl <- findOverlapsRanges(flank(exons, 1, FALSE), frag_exonic)
    excl <- excl[which(spliceR[queryHits(excl)])]
    hits <- setdiff(hits, excl)
    
    exon_index <- as.list(hits)
    
    if (counts) elementLengths(exon_index)
    else exon_index

}

##' Identify fragments with alignments extending across exon/intron boundaries.
##' 
##' @title Compatible fragment counts for splice sites
##' @inheritParams exonCompatible
##' @param splicesites \code{IRanges} of splice sites
##' @param side Character vector indicating whether the spliced boundary
##'   is to the left (\dQuote{L}) or right (\dQuote{R}) of the splice site
##' @param include Character string indicating whether considered fragments
##'   should be all that overlap the splice site (\dQuote{all}), those
##'   that are spliced at the site (\dQuote{spliced}) or those that are
##'   not spliced, i.e. extend into the adjacent intron (\dQuote{unspliced})
##' @return Counts or list of indices of compatible fragments
##' @keywords internal
##' @author Leonard Goldstein

splicesiteOverlap <- function(splicesites, side, frag_exonic,
    frag_intron, include = c("all", "spliced", "unspliced"), counts = TRUE)
{

    include <- match.arg(include)

    if (length(side) == 1) { side <- rep(side, length(splicesites)) }

    ## for each splice site, identify alignments with overlapping exonic region
    ## and no overlapping introns (the second condition is required in case of
    ## inconsistent paired end alignments)

    hits_exonic <- findOverlapsRanges(splicesites, frag_exonic)
    hits_introns <- findOverlapsRanges(splicesites, frag_intron)
    hits <- setdiff(hits_exonic, hits_introns)

    ## to obtain splice site counts that are comparable with junction
    ## counts, require the flanking intronic position to overlap an intron
    ## (include = "spliced"), an exonic region (include = "unspliced")
    ## or either (include = "all")

    ## Note that a fragment can be both spliced and unspliced with respect
    ## to a splice site if paired ends have inconsistent alignments.
    ## Here we consider a fragment spliced if at least one end is spliced,
    ## and unspliced otherwise.

    start <- setNames(c(L = TRUE, R = FALSE)[side], NULL)
    flanking <- flank(splicesites, 1, start = start)
    
    if (include == "all") {
    
        hits_flanking_intron <- findOverlapsRanges(flanking, frag_intron)
        hits_flanking_exonic <- findOverlapsRanges(flanking, frag_exonic)
        hits <- intersect(hits, union(hits_flanking_intron,
            hits_flanking_exonic))

    } else if (include == "spliced") {
        
        hits_flanking_intron <- findOverlapsRanges(flanking, frag_intron)
        hits <- intersect(hits, hits_flanking_intron)

    } else if (include == "unspliced") {

        hits_flanking_intron <- findOverlapsRanges(flanking, frag_intron)
        hits_flanking_exonic <- findOverlapsRanges(flanking, frag_exonic)
        hits <- intersect(hits, hits_flanking_exonic)
        hits <- setdiff(hits, hits_flanking_intron)

    }
        
    splicesite_index <- as.list(hits)

    if (counts) elementLengths(splicesite_index)
    else splicesite_index

}

##' Modified \code{findOverlaps} function for \code{IRanges},
##' \code{IRangesList} objects that behaves analogous to
##' \code{findOverlaps} for \code{GRanges}, \code{GRangesList} objects.
##' 
##' @title Modified \code{findOverlaps} function for \code{IRanges},
##'   \code{IRangesList} objects
##' @param query \code{IRanges} or \code{IRangesList} object
##' @param subject \code{IRanges} or \code{IRangesList} object
##' @param type Passed to \code{findOverlaps}
##' @return \code{Hits} object
##' @keywords internal
##' @author Leonard Goldstein

findOverlapsRanges <- function(query, subject, type = "any")
{

    if (is(query, "IRangesList")) {
    
        query_unlisted <- unlist(query)
        query_togroup <- togroup(query)

    } else {

        query_unlisted <- query
        query_togroup <- seq_along(query)

    }

    if (is(subject, "IRangesList")) {

        subject_unlisted <- unlist(subject)
        subject_togroup <- togroup(subject)

    } else {

        subject_unlisted <- subject
        subject_togroup <- seq_along(subject)
        
    }
    
    if (type == "equal") {

        hits_unlisted <- findMatches(query_unlisted, subject_unlisted)
        
    } else {

        hits_unlisted <- findOverlaps(query_unlisted, subject_unlisted,
            type = type)

    }

    qH <- query_togroup[queryHits(hits_unlisted)]
    sH <- subject_togroup[subjectHits(hits_unlisted)]
    hits <- unique(cbind(qH, sH))

    hits <- new2("Hits",
        queryHits = hits[, 1],
        subjectHits = hits[, 2],
        queryLength = length(query),
        subjectLength = length(subject),
        check = FALSE)

    return(hits)
    
}

splicesiteCounts <- function(x, frag_exonic, frag_intron,
    option = c("junction", "exon"))
{

    option <- match.arg(option)

    include <- switch(option, junction = "all", exon = "spliced")
    
    N_L <- splicesiteOverlap(flank(x, -1, TRUE),
        switch(option, junction = "R", exon = "L"),
        frag_exonic, frag_intron, include)
    N_R <- splicesiteOverlap(flank(x, -1, FALSE),
        switch(option, junction = "L", exon = "R"),
        frag_exonic, frag_intron, include)
    N <- IntegerList(mapply(c, N_L, N_R, SIMPLIFY = FALSE))

    return(N)

}

exonCoverage <- function(exons, exons_i_frag, frag_exonic)
{

    frag_exonic <- reduce(frag_exonic)

    expanded_exon <- factor(togroup(exons_i_frag), seq_along(exons))
    expanded_frag_exonic <- frag_exonic[unlist(exons_i_frag)]

    expanded_exon <- expanded_exon[togroup(expanded_frag_exonic)]
    expanded_frag_exonic <- unlist(expanded_frag_exonic)
    
    irl <- split(expanded_frag_exonic, expanded_exon)
    
    coverage <- coverage(irl, shift = -start(exons) + 1, width = width(exons))
    coverage <- RleList(coverage, compress = TRUE)
    
    return(coverage)
    
}
