##' Identify fragments compatible with splice junctions.
##'
##' @title Compatible fragment counts for splice junctions
##' @inheritParams exonCompatible
##' @param junctions \code{IRanges} of splice junctions
##' @param min_anchor Integer specifying minimum anchor length
##' @return Counts or list of indices of compatible fragments
##' @keywords internal
##' @author Leonard Goldstein

junctionCompatible <- function(junctions, frag_exonic, frag_intron,
    min_anchor, counts = TRUE)
{

    frag_intron <- filterIntrons(frag_intron, frag_exonic, min_anchor)

    introns <- junctions - 1

    hits <- findOverlapsRanges(introns, frag_intron, "equal")
    junction_index <- as.list(hits)

    if (counts) elementNROWS(junction_index)
    else junction_index

}

filterIntrons <- function(frag_intron, frag_exonic, min_anchor)
{

    unlisted_intron <- unlist(frag_intron)

    f_L <- as(flank(unlisted_intron, min_anchor, TRUE), "RangesList")
    f_L <- intersect(f_L, frag_exonic[togroup0(frag_intron)])
    w_L <- sum(width(f_L))

    f_R <- as(flank(unlisted_intron, min_anchor, FALSE), "RangesList")
    f_R <- intersect(f_R, frag_exonic[togroup0(frag_intron)])
    w_R <- sum(width(f_R))

    i <- which(w_L == min_anchor & w_R == min_anchor)

    filtered <- setNames(split(unlisted_intron[i],
        factor(togroup0(frag_intron)[i], seq_along(frag_intron))), NULL)

    return(filtered)

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

    if (length(spliceL) == 1) spliceL <- rep(spliceL, length(exons))
    if (length(spliceR) == 1) spliceR <- rep(spliceR, length(exons))

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

    if (counts) elementNROWS(exon_index)
    else exon_index

}

##' Identify fragments with alignments extending across exon/intron boundaries.
##'
##' @title Compatible fragment counts for splice sites
##' @inheritParams exonCompatible
##' @param splicesites \code{IRanges} of splice sites
##' @param side Character vector indicating whether the spliced boundary
##'   is to the left (\dQuote{L}) or right (\dQuote{R}) of the splice site
##' @param min_anchor Integer specifiying minimum anchor length
##' @param include Character string indicating whether considered fragments
##'   should be all that overlap the splice site (\dQuote{all}), those
##'   that are spliced at the site (\dQuote{spliced}) or those that are
##'   not spliced, i.e. extend into the adjacent intron (\dQuote{unspliced})
##' @return Counts or list of indices of compatible fragments
##' @keywords internal
##' @author Leonard Goldstein

splicesiteOverlap <- function(splicesites, side, frag_exonic, frag_intron,
    min_anchor, include = c("all", "spliced", "unspliced"), counts = TRUE)
{

    include <- match.arg(include)

    if (length(side) == 1) side <- rep(side, length(splicesites))

    start <- setNames(c(L = TRUE, R = FALSE)[side], NULL)
    flanking_intron <- flank(splicesites, 1, start = start)
    hits <- findOverlapsRanges(splicesites, frag_exonic)

    if (include == "spliced" || include == "all") {

        frag_intron <- filterIntrons(frag_intron, frag_exonic, min_anchor)
        tmp <- findOverlapsRanges(flanking_intron, frag_intron)
        hits_spliced <- intersect(hits, tmp)

    }

    if (include == "unspliced" || include == "all") {

        tmp <- findOverlapsRanges(flanking_intron, frag_exonic)
        hits_unspliced <- intersect(hits, tmp)
        qH <- queryHits(hits_unspliced)
        sH <- subjectHits(hits_unspliced)
        intron_anchor <- as(flank(splicesites, min_anchor, start = start),
            "RangesList")
        exonic_anchor <- as(flank(splicesites, -min_anchor, start = start),
            "RangesList")
        w_E <- sum(width(intersect(exonic_anchor[qH], frag_exonic[sH])))
        w_I <- sum(width(intersect(intron_anchor[qH], frag_exonic[sH])))
        hits_unspliced <- hits_unspliced[w_E == min_anchor & w_I == min_anchor]

    }

    if (include == "spliced") {

        hits <- hits_spliced

    } else if (include == "unspliced") {

        hits <- hits_unspliced

    } else if (include == "all") {

        hits <- union(hits_spliced, hits_unspliced)

    }

    splicesite_index <- as.list(hits)

    if (counts) elementNROWS(splicesite_index)
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
        query_togroup <- togroup0(query)

    } else {

        query_unlisted <- query
        query_togroup <- seq_along(query)

    }

    if (is(subject, "IRangesList")) {

        subject_unlisted <- unlist(subject)
        subject_togroup <- togroup0(subject)

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
    hits <- Hits(as.integer(hits[, 1]), as.integer(hits[, 2]),
        length(query), length(subject), sort.by.query = TRUE)

    return(hits)

}

splicesiteCounts <- function(x, frag_exonic, frag_intron, min_anchor,
    option = c("junction", "exon"), include)
{

    option <- match.arg(option)

    N_L <- splicesiteOverlap(flank(x, -1, TRUE),
        switch(option, junction = "R", exon = "L"),
        frag_exonic, frag_intron, min_anchor, include)
    N_R <- splicesiteOverlap(flank(x, -1, FALSE),
        switch(option, junction = "L", exon = "R"),
        frag_exonic, frag_intron, min_anchor, include)
    N <- IntegerList(mapply(c, N_L, N_R, SIMPLIFY = FALSE))

    return(N)

}

exonCoverage <- function(exons, exons_i_frag, frag_exonic)
{

    expanded_exon <- factor(togroup0(exons_i_frag), seq_along(exons))
    expanded_frag_exonic <- frag_exonic[unlist(exons_i_frag)]

    expanded_exon <- expanded_exon[togroup0(expanded_frag_exonic)]
    expanded_frag_exonic <- unlist(expanded_frag_exonic)

    irl <- split(expanded_frag_exonic, expanded_exon)
    coverage <- coverage(irl, shift = -start(exons) + 1, width = width(exons))
    ## coverage() returns a SimpleRleList
    ## need RleList() to obtain a CompressedRleList
    coverage <- RleList(coverage)

    return(coverage)

}
