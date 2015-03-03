##' Features in \code{query} are annotated with respect to transcript
##' features in \code{subject}. 
##'
##' Annotation happens at two levels: For feature-centric annotation,
##' query features are assigned all transcript names associated with any
##' matching subject features. For gene-centric annotation, query features
##' are assigned all gene names associated with subject features that are
##' part of the same gene (connected component in the splice graph)
##' as any matching query features. 
##'
##' Feature matching is performed as follows:
##' Query splice junctions are matched with identical subject splice
##' junctions. Query splice sites are matched with splice sites
##' implied by subject splice junctions. Query exon bins are matched with
##' overlapping subject exons. Spliced boundaries of query exon bins
##' must match spliced subject exon boundaries. Query exon bins cannot
##' extend across spliced subject exon boundaries.
##' 
##' @title Annotation with respect to transcript features
##' @param query \code{SGFeatures}, \code{TxVariants},
##'   \code{SGFeatureCounts} or \code{TxVariantCounts} object
##' @param subject \code{TxFeatures} object
##' @return \code{query} with updated \code{txName}, \code{geneName}
##'   column slots
##' @examples
##' sgf_annotated <- annotate(sgf, txf)
##' txv_annotated <- annotate(txv, txf)
##' @author Leonard Goldstein

annotate <- function(query, subject)
{

    if (!is(subject, "TxFeatures")) {

        stop("subject must be a TxFeatures object")

    }

    if (is(query, "SGFeatures")) {

        query <- annotateFeatures(query, subject)
        
    } else if (is(query, "TxVariants")) {

        query_class <- class(query)
        query_mcols <- mcols(query)
        query_unlisted <- unlist(query, use.names = FALSE)
        query_unlisted <- annotate(query_unlisted, subject)
        query <- relist(query_unlisted, query)
        mcols(query) <- query_mcols
        query <- new(query_class, query)
        query <- annotatePaths(query)
        
    } else if (is(query, "Counts")) {

        rd <- rowRanges(query)
        rd <- annotate(rd, subject)
        rowRanges(query) <- rd

    }
        
    return(query)

}

annotateFeatures <- function(query, subject)
{

    i_Z <- which(type(subject) %in% c("F", "L"))

    if (length(i_Z) > 0) {

        subject <- c(subject[-i_Z], mergeExonsTerminal(subject[i_Z], 1))

    }

    if (is(query, "TxFeatures")) {

        hits <- matchTxFeatures(query, subject)

    } else if (is(query, "SGFeatures")) {

        hits <- matchSGFeatures(query, subject)

    }

    qH <- queryHits(hits)
    sH <- subjectHits(hits)

    for (option in c("feature", "gene")) {

        slot_id <- switch(option, feature = "featureID", gene = "geneID")
        slot_ann <- switch(option, feature = "txName", gene = "geneName")
        query_id <- factor(slot(query, slot_id))
        subject_ann <- slot(subject, slot_ann)
        id_ann <- collapseCharacterList(subject_ann[sH], query_id[qH])
        query_ann <- id_ann[match(query_id, names(id_ann))]
        slot(query, slot_ann) <- setNames(query_ann, NULL)
        
    }

    if (!is.null(mcols(subject)$frame)) {

        query <- annotateFrame(query, subject)

    }
    
    return(query)
    
}

annotateSegments <- function(ids, features)
{

    ## ids is a character vector with comma-separatated lists of
    ## feature IDs, or txNames in the format {tx1;...;txn}
  
    out <- vector("list", length(ids))

    segment_ids <- strsplit(ids, ",")
    segment_n <- elementLengths(segment_ids)

    ids <- unlist(segment_ids)
    ids_segment <- togroup(segment_ids)
    ids_ann <- vector("list", length(ids))
    ids_ann <- as(ids_ann, "CompressedCharacterList")
    
    i <- grep("^\\d+$", ids)

    if (length(i) > 0) {

        ids_ann[i] <- txName(features)[match(ids[i], featureID(features))]

    }

    i <- grep("^\\{\\S*\\}$", ids)

    if (length(i) > 0) {

        tmp <- ids[i]
        tmp <- gsub("\\{|\\}", "", tmp)
        tmp <- strsplit(tmp, ";", fixed = TRUE)
        tmp <- as(tmp, "CompressedCharacterList")
        ids_ann[i] <- tmp

    }

    ann <- unlist(ids_ann)
    ann_id <- togroup(ids_ann)
    ann_segment <- ids_segment[ann_id]

    if (length(ann) == 0) {

        return(out)

    }
    
    segment_ann_n <- table(paste0(ann_segment, ":", ann))
    x_segment <- sub(":\\S+$", "", names(segment_ann_n))
    x_ann <- sub("^\\S+:", "", names(segment_ann_n))

    i <- which(segment_ann_n == segment_n[as.integer(x_segment)])

    if (length(i) > 0) {

        x_segment <- x_segment[i]
        x_ann <- x_ann[i]
        segment_ann <- split(x_ann, x_segment)
        out[as.integer(names(segment_ann))] <- segment_ann

    }
        
    return(out)

}

annotatePaths <- function(paths)
{

    features <- unlist(paths)
    
    x <- featureID(paths)
    i <- grep("(", x, fixed = TRUE)
    
    while (length(i) > 0) {

        ## find inner event
        m <- regexpr("\\([^\\(\\)]+\\)", x[i])
        l <- attr(m, "match.length")
        
        ## split at event
        u <- substr(x[i], 1, m - 1)
        b <- substr(x[i], m + 1, m + l - 2)
        v <- substr(x[i], m + l, nchar(x)[i])
        
        ## annotate event
        list_ids <- strsplit(b, "|", fixed = TRUE)
        ids <- unlist(list_ids)
        ids_path <- i[togroup(list_ids)]
        ids_ann <- annotateSegments(ids, features)
        ann <- unlist(ids_ann)
        ann_id <- togroup(ids_ann)
        ann_path <- ids_path[ann_id]

        path_ann <- vector("list", length(i))

        if (!is.null(ann)) {
          
            tmp <- tapply(ann, ann_path, unique, simplify = FALSE)
            path_ann[match(names(tmp), i)] <- tmp

        }

        path_ann <- paste0("{", unstrsplit(path_ann, ";"), "}")

        x[i] <- paste0(u, path_ann, v)
        i <- grep("(", x, fixed = TRUE)

    }

    out <- annotateSegments(x, features)
    out <- as(out, "CompressedCharacterList")
    txName(paths) <- out
    geneName(paths) <- geneName(features)[match(seq_along(paths),
        togroup(paths))]

    return(paths)
    
}

plintersect <- function(x, y)
{
  
    n <- length(x)
    ix <- paste0(togroup(x), ":", unlist(x))
    iy <- paste0(togroup(y), ":", unlist(y))
    ix <- ix[ix %in% iy]
    i <- sub(":\\S+$", "", ix)
    x <- sub("^\\S+:", "", ix)
    z <- split(x, i)
    z <- z[match(seq_len(n), names(z))]
    names(z) <- NULL

    return(z)

}

annotateFrame <- function(query, subject)
{
    
    hits <- matchExon(query, subject)
    qH <- queryHits(hits)
    sH <- subjectHits(hits)

    i_pos <- which(strand(query[qH]) == "+")
    i_neg <- which(strand(query[qH]) == "-")

    offset <- rep(NA_integer_, length(hits))
    offset[i_pos] <- start(query)[qH][i_pos] - start(subject)[sH][i_pos]
    offset[i_neg] <- end(subject)[sH][i_neg] - end(query)[qH][i_neg]
    
    sH_frame <- mcols(subject)$frame[sH]
    qH_frame <- nextFrame(sH_frame, offset)

    qH_frameStart <- qH_frame
    qH_frameEnd <- nextFrame(qH_frame, width(query)[qH])

    sH_cdsStart <- mcols(subject)$cdsStart[sH]
    sH_cdsEnd <- mcols(subject)$cdsEnd[sH]
    
    qH_cdsStart_unlisted <- unlist(sH_cdsStart) - offset[togroup(sH_cdsStart)]
    i_na <- which(qH_cdsStart_unlisted < 0 |
        qH_cdsStart_unlisted > width(query)[qH][togroup(sH_cdsStart)])
    qH_cdsStart_unlisted[i_na] <- NA_integer_
    
    qH_cdsEnd_unlisted <- unlist(sH_cdsEnd) - offset[togroup(sH_cdsEnd)]
    i_na <- which(qH_cdsEnd_unlisted < 0 |
        qH_cdsEnd_unlisted > width(query)[qH][togroup(sH_cdsEnd)])
    qH_cdsEnd_unlisted[i_na] <- NA_integer_
    
    qH_frameStart_unlisted <- unlist(qH_frameStart)
    qH_frameStart_unlisted[!is.na(qH_cdsStart_unlisted)] <- -1
    
    qH_frameEnd_unlisted <- unlist(qH_frameEnd)
    qH_frameEnd_unlisted[!is.na(qH_cdsEnd_unlisted)] <- -1
    
    x <- paste(qH_frameStart_unlisted, qH_frameEnd_unlisted, sep = ",")
    x <- tapply(x, qH[togroup(sH_cdsStart)], setdiff, "NA,NA",
        simplify = FALSE)
    y <- strsplit(unlist(x), ",")

    frameStart <- split(suppressWarnings(as.integer(sapply(y, "[", 1))),
        names(x)[togroup(x)])
    frameStart <- IntegerList(frameStart)

    frameEnd <- split(suppressWarnings(as.integer(sapply(y, "[", 2))),
        names(x)[togroup(x)])
    frameEnd <- IntegerList(frameEnd)

    query_frameStart <- IntegerList(vector("list", length(query)))
    query_frameStart[as.integer(names(frameStart))] <- frameStart
    mcols(query)$frameStart <- setNames(query_frameStart, NULL)

    query_frameEnd <- IntegerList(vector("list", length(query)))
    query_frameEnd[as.integer(names(frameEnd))] <- frameEnd
    mcols(query)$frameEnd <- setNames(query_frameEnd, NULL)

    return(query)
    
}

matchTxFeatures <- function(query, subject)
{

    q <- feature2name(query, TRUE)
    s <- feature2name(subject, TRUE)
    s2i <- split(seq_along(s), s)

    qH <- which(q %in% names(s2i))
    sH <- match(q[qH], names(s2i))

    qH <- qH[togroup(sH)]
    sH <- unlist(sH)
    
    new2("Hits",
         queryHits = qH,
         subjectHits = sH,
         queryLength = length(query),
         subjectLength = length(subject),
         check = FALSE)

}

matchSGFeatures <- function(query, subject) {

    Reduce(union, list(
        matchJunction(query, subject),
        matchExon(query, subject),
        matchSplice(query, subject, "D"),
        matchSplice(query, subject, "A")
    ))
    
}

matchJunction <- function(query, subject) {

    i_q <- which(type(query) == "J")
    i_s <- which(type(subject) == "J")

    hits <- findMatches(query[i_q], subject[i_s])

    new2("Hits",
        queryHits = i_q[queryHits(hits)],
        subjectHits = i_s[subjectHits(hits)],
        queryLength = length(query),
        subjectLength = length(subject),
        check = FALSE)
    
}

matchExon <- function(query, subject) {

    i_q <- which(type(query) == "E")    
    i_s <- which(type(subject) %in% c("I", "F", "L", "U"))

    q2 <- query[i_q]
    s2 <- subject[i_s]
    
    hits <- findOverlaps(q2, s2)
    
    qH <- queryHits(hits)
    sH <- subjectHits(hits)

    ## exclude hits with inconsistent 5' spliced boundary
    
    i <- which(splice5p(q2)[qH])

    if (length(i) > 0) {

        q2_5p <- start(flank(q2[qH][i], -1, TRUE))
        s2_5p <- start(flank(s2[sH][i], -1, TRUE))
        s2_type <- type(s2)[sH][i]
        excl_q_5p <- i[which(q2_5p != s2_5p | s2_type %in% c("F", "U"))]

    } else { excl_q_5p <- integer() }

    ## exclude hits with inconsistent 3' spliced boundary

    i <- which(splice3p(q2)[qH])

    if (length(i) > 0) {
    
        q2_3p <- start(flank(q2[qH][i], -1, FALSE))
        s2_3p <- start(flank(s2[sH][i], -1, FALSE))
        s2_type <- type(s2)[sH][i]
        excl_q_3p <- i[which(q2_3p != s2_3p | s2_type %in% c("L", "U"))]

    } else { excl_q_3p <- integer() }

    ## exclude hits with query overlapping 5' flanking intron
    
    i <- which(type(s2)[sH] %in% c("I", "L"))

    if (length(i) > 0) {

        hits_2 <- findOverlaps(q2[qH][i], flank(s2[sH][i], 1, TRUE))
        hits_2 <- hits_2[queryHits(hits_2) == subjectHits(hits_2)]
        excl_s_5p <- i[queryHits(hits_2)]

    } else { excl_s_5p <- integer() }

    ## exclude hits with query overlapping 3' flanking intron

    i <- which(type(s2)[sH] %in% c("F", "I"))

    if (length(i) > 0) {

        hits_2 <- findOverlaps(q2[qH][i], flank(s2[sH][i], 1, FALSE))
        hits_2 <- hits_2[queryHits(hits_2) == subjectHits(hits_2)]
        excl_s_3p <- i[queryHits(hits_2)]

    } else { excl_s_3p <- integer() }
    
    excl <- unique(c(excl_q_5p, excl_q_3p, excl_s_5p, excl_s_3p))

    if (length(excl) > 0) { hits <- hits[-excl] }
    
    new2("Hits",
        queryHits = i_q[queryHits(hits)],
        subjectHits = i_s[subjectHits(hits)],
        queryLength = length(query),
        subjectLength = length(subject),
        check = FALSE)

}

matchSplice <- function(query, subject, type = c("D", "A")) {

    type <- match.arg(type)
    
    i_q <- which(type(query) == type)    
    q <- query[i_q]

    if (type == "D") {

        i_s_J <- which(type(subject) == "J")
        s_J <- flank(subject[i_s_J], -1, TRUE)

        i_s_E <- which(type(subject) %in% c("I", "F"))
        s_E <- flank(subject[i_s_E], -1, FALSE)
      
    } else if (type == "A") {

        i_s_J <- which(type(subject) == "J")
        s_J <- flank(subject[i_s_J], -1, FALSE)

        i_s_E <- which(type(subject) %in% c("I", "L"))
        s_E <- flank(subject[i_s_E], -1, TRUE)

    }

    i_s <- c(i_s_J, i_s_E)
    s <- c(s_J, s_E)
    
    hits <- findMatches(q, s)

    new2("Hits",
        queryHits = i_q[queryHits(hits)],
        subjectHits = i_s[subjectHits(hits)],
        queryLength = length(query),
        subjectLength = length(subject),
        check = FALSE)

}
