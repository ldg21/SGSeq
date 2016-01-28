##' Features in \code{query} are annotated with respect to transcript
##' features in \code{subject}.
##'
##' Annotation is performed at the gene and transcript level.
##' For transcript-level annotation, query features are assigned all
##' transcript names associated with matching subject features.
##' For gene-level annotation, query features are assigned all gene
##' names associated with subject features that belong to the same
##' gene (connected component of the splice graph) as matching
##' query features.
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
##' @param query \code{SGFeatures}, \code{SGVariants},
##'   \code{SGFeatureCounts} or \code{SGVariantCounts} object
##' @param subject \code{TxFeatures} object
##' @return \code{query} with updated \code{txName}, \code{geneName}
##'   column slots
##' @examples
##' sgf_annotated <- annotate(sgf_pred, txf_ann)
##' sgv_annotated <- annotate(sgv_pred, txf_ann)
##' @author Leonard Goldstein

annotate <- function(query, subject)
{

    if (!is(subject, "TxFeatures")) {

        stop("subject must be a TxFeatures object")

    }

    if (is(query, "SGFeatures")) {

        query <- annotateFeatures(query, subject)

    } else if (is(query, "SGVariants")) {

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

    i <- which(type(subject) %in% c("F", "L"))

    if (length(i) > 0) {

        subject <- c(subject[-i], mergeExonsTerminal(subject[i], 1))

    }

    if (is(query, "TxFeatures")) {

        hits <- matchTxFeatures(query, subject)

    } else if (is(query, "SGFeatures")) {

        hits <- matchSGFeatures(query, subject)

    }

    qH <- queryHits(hits)
    sH <- subjectHits(hits)

    for (option in c("tx", "gene")) {

        q_id <- factor(slot(query, "featureID"))
        s_ann <- slot(subject, paste0(option, "Name"))
        id_ann <- splitCharacterList(s_ann[sH], q_id[qH])
        q_ann <- setNames(id_ann[match(q_id, names(id_ann))], NULL)
        slot(query, paste0(option, "Name")) <- q_ann

    }

    if (is(query, "SGFeatures")) {

        query <- propagateAnnotation(query)

    }

    return(query)

}

propagateAnnotation <- function(query)
{

    g <- spliceGraph(query)

    ## Remove annotated parts of the splice graph. Specifically,
    ## remove annotated edges with annotated 'to' and 'from' vertices.
    ## Then remove edge-free vertices (vertices that are annotated and with
    ## all incoming and outgoing edges annotated in the original graph)

    gv <- nodes(g)
    gd <- edges(g)

    gd_ann <- elementNROWS(gd$geneName) > 0
    i <- match(gv$featureID[match(gd$from, gv$name)], featureID(query))
    from_ann <- elementNROWS(geneName(query))[i] > 0
    i <- match(gv$featureID[match(gd$to, gv$name)], featureID(query))
    to_ann <- elementNROWS(geneName(query))[i] > 0

    excl <- which(
        (is.na(from_ann) | is.na(to_ann)) |
        (gd_ann & from_ann & to_ann))

    if (length(excl) > 0)  g <- delete.edges(g, excl)

    gd <- edges(g)

    excl <- which(!gv$name %in% c(gd$from, gd$to))

    if (length(excl) > 0) g <- delete.vertices(g, excl)

    ## Unnannotated parts of the splice graph are assigned gene names
    ## of connected annotated vertices or edges.

    gv <- nodes(g)
    gd <- edges(g)

    gv_geneName <- geneName(query)[match(gv$featureID, featureID(query))]
    gd_geneName <- geneName(query)[match(gd$featureID, featureID(query))]

    gv_cluster <- as.character(clusters(g)$membership)
    gd_cluster <- gv_cluster[match(gd$from, gv$name)]

    ann <- DataFrame(
        featureID = c(gv$featureID, gd$featureID),
        geneName = c(gv_geneName, gd_geneName),
        cluster = c(gv_cluster, gd_cluster))

    cluster_geneName <- splitCharacterList(ann$geneName, factor(ann$cluster))
    ann <- ann[elementNROWS(ann$geneName) == 0, ]
    i <- match(ann$cluster, names(cluster_geneName))
    ann$geneName <- cluster_geneName[i]
    i <- match(ann$featureID, featureID(query))
    geneName(query)[i] <- ann$geneName

    return(query)

}

annotateSGSegments <- function(ids, features)
{

    ## ids is a character vector with comma-separatated lists of
    ## feature IDs, or txNames in the format {tx1;...;txn}

    out <- vector("list", length(ids))

    segment_ids <- strsplit(ids, ",")
    segment_n <- elementNROWS(segment_ids)

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
    x_segment <- sapply(strsplit(names(segment_ann_n), ":"), "[", 1)
    x_ann <- sub("^\\d+:", "", names(segment_ann_n))

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

    ## encode txName
    txName_unique <- unique(unlist(txName(features)))
    txName(features) <- relist(as.character(match(unlist(txName(features)),
        txName_unique)), txName(features))

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
        ids_ann <- annotateSGSegments(ids, features)
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

    out <- annotateSGSegments(x, features)
    out <- as(out, "CompressedCharacterList")

    ## decode txName
    out <- relist(txName_unique[as.integer(unlist(out))], out)

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

matchExon <- function(query, subject, stringent = TRUE) {

    i_q <- which(type(query) == "E")
    i_s <- which(type(subject) %in% c("I", "F", "L", "U"))

    q2 <- query[i_q]
    s2 <- subject[i_s]

    hits <- findOverlaps(q2, s2)

    qH <- queryHits(hits)
    sH <- subjectHits(hits)

    hits <- new2("Hits",
        queryHits = i_q[queryHits(hits)],
        subjectHits = i_s[subjectHits(hits)],
        queryLength = length(query),
        subjectLength = length(subject),
        check = FALSE)

    if (!stringent) { return(hits) }

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

    return(hits)

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
