##' Convert a \code{TxDb} object or a \code{GRangesList} of exons
##' grouped by transcripts to a \code{TxFeatures} object.
##'
##' In the returned \code{TxFeatures} object, column slot \code{txName}
##' is based on either \code{tx_name} in the \code{TxDb} object
##' or \code{names(x)} if \code{x} is a \code{GRangesList}. 
##' 
##' @title Convert to TxFeatures object
##' @param x \code{TxDb} object, or \code{GRangesList} of exons
##'   grouped by transcripts
##' @return A \code{TxFeatures} object
##' @examples
##' gr <- GRanges(c(1, 1), IRanges(c(1, 201), c(100, 300)), c("+", "+"))
##' grl <- split(gr, 1)
##' txf <- convertToTxFeatures(grl)
##' @author Leonard Goldstein

convertToTxFeatures <- function(x)
{
    
    if (is(x, "TxDb")) {
        
        tx <- exonsBy(x, "tx", use.names = TRUE)
        
    } else if (is(x, "GRangesList")) {

        if (!exonsOnSameChromAndStrand(x)) {
            
            stop("All ranges in the same element of x must be on the same
                chromosome and strand")

        }
        
        tx <- x
        
        if (is.null(names(tx))) { names(tx) <- seq_along(tx) }

    } else {

        stop("x must be a TxDb or GRangesList object")

    }
    
    ## reduce ranges to ensure there are no adjacent exons

    tx <- reduce(tx)

    ## remove tx with duplicate names
    
    dup <- names(tx)[duplicated(names(tx))]
    
    if (length(dup) > 0) {

        remove <- which(names(tx) %in% dup)
        tx <- tx[-remove]
        warning(paste("Removed", length(remove),
            "transcripts with duplicate names"))

    }

    exons <- extractExons(tx)
    junctions <- extractJunctions(tx)
    features <- c(exons, junctions)
    features <- collapseFeatures(features)
    features <- TxFeatures(features)

    if (is(x, "TxDb")) {

        tx_name <- txName(features)
        tx_name_unlisted <- unlist(tx_name)
        gene_id_unlisted <- select(x, tx_name_unlisted, "GENEID",
            "TXNAME")$GENEID
        i_not_na <- which(!is.na(gene_id_unlisted))
        i_gene_id <- split(gene_id_unlisted[i_not_na],
            togroup(tx_name)[i_not_na])
        i_gene_id <- CharacterList(i_gene_id)
        gene_id <- CharacterList(vector("list", length(features)))
        gene_id[as.integer(names(i_gene_id))] <- i_gene_id
        gene_id <- unique(gene_id)
        geneName(features) <- gene_id
        
    }
    
    return(features)
        
}

exonsOnSameChromAndStrand <- function(x)
{

    sn <- unique(CharacterList(seqnames(x)))
    st <- unique(CharacterList(strand(x)))
                
    if (any(elementLengths(sn) > 1) || any(elementLengths(st) > 1)) {

        valid <- FALSE
        
    } else {

        valid <- TRUE

    }

    return(valid)

}

extractExons <- function(tx)
{

    tx <- reorderFeatures(tx)
    exons <- setNames(unlist(tx), NULL)
    exons_tx <- togroup(tx)

    tx_i_exon <- IntegerList(split(seq_along(exons), exons_tx))
    i_X <- pfirst(tx_i_exon)
    i_Y <- plast(tx_i_exon)
    
    type <- rep("I", length(exons))
    type[i_X] <- "F"
    type[i_Y] <- "L"
    type[intersect(i_X, i_Y)] <- "U"
    
    mcols(exons) <- DataFrame(
        type = type,
        txName = names(tx)[exons_tx])

    if (!is.null(mcols(tx)$cdsStart) && !is.null(mcols(tx)$cdsEnd) &&
        !is.null(mcols(tx)$status)) {

        exons <- exonFrame(exons, exons_tx, tx)

    }
    
    return(exons)
    
}

exonFrame <- function(ex, ex_tx, tx)
{
    
    tx_chrom <- as.character(seqnames(unlist(range(tx))))
    tx_strand <- as.character(strand(unlist(range(tx))))
    tx_cdsStart <- mcols(tx)$cdsStart
    tx_cdsEnd <- mcols(tx)$cdsEnd
    tx_status <- mcols(tx)$status
    
    i_tx_coding <- which(!is.na(tx_cdsStart) & !is.na(tx_cdsEnd)) 
    
    tx_coding_cds <- GRanges(
        tx_chrom[i_tx_coding],
        IRanges(tx_cdsStart[i_tx_coding], tx_cdsEnd[i_tx_coding]),
        tx_strand[i_tx_coding]
    )

    hits <- findOverlaps(ex, tx_coding_cds)
    hits <- hits[ex_tx[queryHits(hits)] == i_tx_coding[subjectHits(hits)]]

    i_ex_coding <- queryHits(hits)
    ex_coding <- ex[i_ex_coding]
    ex_coding_tx <- ex_tx[i_ex_coding]
    
    tx_i_ex_coding <- IntegerList(split(seq_along(ex_coding), ex_coding_tx))
    i_5p <- pfirst(tx_i_ex_coding, use_names = TRUE)
    i_3p <- plast(tx_i_ex_coding, use_names = TRUE)
    
    i_5p_pos <- i_5p[which(strand(ex_coding)[i_5p] == "+")]
    i_5p_neg <- i_5p[which(strand(ex_coding)[i_5p] == "-")]

    i_3p_pos <- i_3p[which(strand(ex_coding)[i_3p] == "+")]
    i_3p_neg <- i_3p[which(strand(ex_coding)[i_3p] == "-")]
        
    ## set frame for coding exons
    
    ex_coding_cds <- restrict(ex_coding, start = tx_cdsStart[ex_coding_tx],
        end = tx_cdsEnd[ex_coding_tx])
    ex_coding_cds_length <- unlist(tapply(width(ex_coding_cds),
        ex_coding_tx, cumsum), use.names = FALSE)
    ex_coding_cds_offset <- c(NA_integer_,
        ex_coding_cds_length[-length(ex_coding_cds_length)])
    ex_coding_cds_offset[i_5p_pos] <- start(ex_coding)[i_5p_pos] -
        start(ex_coding_cds)[i_5p_pos]
    ex_coding_cds_offset[i_5p_neg] <- end(ex_coding_cds)[i_5p_neg] -
        end(ex_coding)[i_5p_neg]
    ex_coding_frame <- ex_coding_cds_offset %% 3
    
    ## exclude transcripts with invalid coding sequence
    
    tx_invalid <- as.integer(names(i_3p))[which(
        ex_coding_cds_length[i_3p] %% 3 != 0)]
    ex_coding_frame[ex_coding_tx %in% tx_invalid] <- NA_integer_

    ## set frame for all exons
    
    ex_frame <- ifelse(is.na(tx_status[ex_tx]), NA_integer_, -1)
    ex_frame[i_ex_coding] <- ex_coding_frame
    mcols(ex)$frame <- as.integer(ex_frame)

    ## set cdsStart / cdsEnd for all exons

    ex_cdsStart <- rep(NA_integer_, length(ex))
    ex_cdsStart[i_ex_coding][i_5p_pos] <- tx_cdsStart[ex_coding_tx][i_5p_pos] -
        start(ex_coding)[i_5p_pos]
    ex_cdsStart[i_ex_coding][i_5p_neg] <- end(ex_coding)[i_5p_neg] -
        tx_cdsEnd[ex_coding_tx][i_5p_neg]
    ex_cdsStart[ex_tx %in% tx_invalid] <- NA_integer_
    mcols(ex)$cdsStart <- as.integer(ex_cdsStart)
    
    ex_cdsEnd <- rep(NA_integer_, length(ex))
    ex_cdsEnd[i_ex_coding][i_3p_pos] <- tx_cdsEnd[ex_coding_tx][i_3p_pos] -
        start(ex_coding)[i_3p_pos] + 1
    ex_cdsEnd[i_ex_coding][i_3p_neg] <- end(ex_coding)[i_3p_neg] -
        tx_cdsStart[ex_coding_tx][i_3p_neg] + 1
    ex_cdsEnd[ex_tx %in% tx_invalid] <- NA_integer_
    mcols(ex)$cdsEnd <- as.integer(ex_cdsEnd)

    return(ex)
    
}

extractJunctions <- function(tx)
{
    
    introns <- psetdiff(range(tx), tx)
    junctions <- setNames(unlist(introns), NULL) + 1
    mcols(junctions) <- DataFrame(
        type = rep("J", length(junctions)),
        txName = names(introns)[togroup(introns)])

    if (!is.null(mcols(tx)$cdsStart) && !is.null(mcols(tx)$cdsEnd) &&
        !is.null(mcols(tx)$status)) {

        for (k in c("frame", "cdsStart", "cdsEnd")) {
            
            mcols(junctions)[[k]] <- rep(NA_integer_, length(junctions))

        }

    }

    return(junctions)
    
}

collapseFeatures <- function(features)
{
        
    features_name <- factor(feature2name(features))
    collapsed <- features[match(levels(features_name), features_name)]
    mcols(collapsed)$txName <- setNames(collapseCharacterList(
        mcols(features)$txName, features_name), NULL)

    if (!is.null(mcols(features)$frame)) {

        cols <- c("frame", "cdsStart", "cdsEnd")

        x <- paste(
            mcols(features)$frame,
            mcols(features)$cdsStart,
            mcols(features)$cdsEnd,
            sep = ",")

        x <- tapply(x, features_name, setdiff, "NA,NA,NA", simplify = FALSE)
        y <- strsplit(unlist(x), ",")

        for (j in seq_along(cols)) {
        
            z <- relist(suppressWarnings(as.integer(sapply(y, "[", j))), x)
            z <- IntegerList(z)
            mcols(collapsed)[[cols[j]]] <- setNames(z, NULL)

        }
        
    }

    collapsed <- sort(collapsed)

    return(collapsed)

}

##' Convert transcript features, predicted from RNA-seq data or
##' extracted from transcript annotation, to splice graph features. 
##'
##' Splice junctions are unaltered. Exons are disjoined into
##' non-overlapping exon bins. Adjacent exon bins without a splice site
##' at the shared boundary are merged. All exon bins are assigned type
##' \dQuote{E}.
##' 
##' Entries for splice donor and acceptor sites (positions immediately
##' upstream and downstream of introns, respectively) are added.
##'
##' In the returned \code{SGFeatures} object, column slots \code{splice5p}
##' and \code{splice3p} indicate whether compatibility with an exon bin
##' requires a fragment to be spliced at the 5' or 3' boundary, respectively.
##' \code{splice5p} (\code{splice3p}) is \code{TRUE} if the first
##' (last) position of the exon coincides with a splice acceptor (donor),
##' and it is not adjacent to a neighboring exon bin.
##'
##' Each feature is assigned a unique feature and gene identifier,
##' stored in column slots \code{featureID} and \code{geneID},
##' respectively. The latter indicates features that belong to the
##' same gene, represented by a connected component in the splice graph.
##' 
##' @title Convert transcript features to splice graph features
##' @param x \code{TxFeatures} object
##' @param coerce Logical indicating whether transcript features
##'   should be coerced to splice graph features without disjoining
##'   exons and omitting splice donor and acceptor sites
##' @return An \code{SGFeatures} object
##' @examples
##' sgf <- convertToSGFeatures(txf)
##' @author Leonard Goldstein

convertToSGFeatures <- function(x, coerce = FALSE)
{

    if (!is(x, "TxFeatures")) {

        stop("x must be a TxFeatures object")

    }

    features <- granges(x)
    mcols(features)$type <- as.character(type(x))

    if (coerce) {

        splice5p <- mcols(features)$type %in% c("I", "L")
        splice3p <- mcols(features)$type %in% c("I", "F")

        splice5p[mcols(features)$type == "J"] <- NA
        splice3p[mcols(features)$type == "J"] <- NA
        
        mcols(features)$type[mcols(features)$type != "J"] <- "E"
        mcols(features)$splice5p <- splice5p
        mcols(features)$splice3p <- splice3p

    } else {

        junctions <- features[mcols(features)$type == "J"]
        exons <- features[mcols(features)$type %in% c("I", "F", "L", "U")]
        features <- processFeatures(exons, junctions)

    }

    features <- addFeatureID(features)
    features <- addGeneID(features)
    
    features <- SGFeatures(features)
    features <- annotate(features, x)
    
    return(features)
    
}

processFeatures <- function(exons, junctions)
{

    D <- unique(flank(junctions, -1, TRUE))
    mcols(D)$type <- "D"
    A <- unique(flank(junctions, -1, FALSE))
    mcols(A)$type <- "A"
    splicesites <- c(D, A)
    other <- c(junctions, splicesites)
    
    ## Disjoin exons into non-overlapping exon bins,
    ## merge adjacent exon bins that do not have a
    ## splice site at the shared boundary
    
    exons <- disjoin(exons)
    exons_start <- flank(exons, -1, TRUE)
    exons_end <- flank(exons, -1, FALSE)    
    
    i_q <- which(!exons_end %over% splicesites)
    i_s <- which(!exons_start %over% splicesites)    
    
    ol <- findOverlaps(suppressWarnings(flank(exons[i_q], 1, FALSE)),
        exons_start[i_s])
    
    if (length(ol) > 0) {
        
        qH <- i_q[queryHits(ol)]
        sH <- i_s[subjectHits(ol)]
        i_to_be_merged <- union(qH, sH)
        d <- data.frame(from = qH, to = sH)
        v <- data.frame(name = i_to_be_merged)
        g <- graph.data.frame(d = d, directed = TRUE, vertices = v)
        k <- clusters(g)$membership
        exons_to_be_merged <- split(exons[i_to_be_merged], k)
        exons_merged <- unlist(reduce(exons_to_be_merged))
        
        if (length(exons_to_be_merged) != length(exons_merged)) {
            
            stop("cannot merge non-adjacent exons")
            
        }
        
        exons <- c(exons[-i_to_be_merged], exons_merged)
        
    }
    
    exons_start <- flank(exons, -1, TRUE)
    exons_end <- flank(exons, -1, FALSE)    
    
    splice5p <- rep(FALSE, length(exons))
    i_spliced <- unique(queryHits(findOverlaps(exons_start, A)))
    i_adjacent <- unique(queryHits(findOverlaps(
        suppressWarnings(flank(exons, 1, TRUE)), exons)))
    splice5p[setdiff(i_spliced, i_adjacent)] <- TRUE
    
    splice3p <- rep(FALSE, length(exons))
    i_spliced <- unique(queryHits(findOverlaps(exons_end, D)))
    i_adjacent <- unique(queryHits(findOverlaps(
        suppressWarnings(flank(exons, 1, FALSE)), exons)))
    splice3p[setdiff(i_spliced, i_adjacent)] <- TRUE
    
    mcols(exons)$type <- "E"
    mcols(exons)$splice5p <- splice5p
    mcols(exons)$splice3p <- splice3p
    
    mcols(other)$splice5p <- NA
    mcols(other)$splice3p <- NA

    ## combine exons and other features
    
    features <- setNames(c(exons, other), NULL)
    features <- sort(features)

    return(features)
    
}

addFeatureID <- function(features)
{
    
    mcols(features)$featureID <- seq_along(features)

    return(features)
    
}

addGeneID <- function(features)
{
    
    g <- spliceGraph(features)
    gd <- edges(g)
    gv <- nodes(g)

    geneID <- rep(NA, length(features))
    geneID[as.integer(gd$featureID)] <- gd$geneID
    i <- which(!is.na(gv$featureID))
    geneID[as.integer(gv$featureID)[i]] <- gv$geneID[i]
    mcols(features)$geneID <- geneID

    return(features)

}
    
convertToTxSegments <- function(x, cores = 1)
{

    x_segmentID <- findSegments(x, cores)

    i <- which(!is.na(x_segmentID))
    segments <- split(x[i], x_segmentID[i])
    segments <- reorderFeatures(segments)

    list_type <- split(as.character(type(unlist(segments))),
        togroup(segments))
    list_featureID <- CharacterList(split(featureID(unlist(segments)),
        togroup(segments)))
    
    gd <- edges(spliceGraph(x))

    i_first <- match(pfirst(list_featureID), gd$featureID)
    i_last <- match(plast(list_featureID), gd$featureID)
    
    mcols(segments) <- DataFrame(
        from = gd$from[i_first],
        to = gd$to[i_last],
        type = unstrsplit(list_type, ""),
        splice5p = gd$splice5p[i_first],
        splice3p = gd$splice3p[i_last],
        featureID = unstrsplit(list_featureID, ","),
        segmentID = seq_along(segments),
        geneID = gd$geneID[i_first])

    segments <- TxSegments(segments)
    segments <- annotatePaths(segments)
    
    return(segments)

}
