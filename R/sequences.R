getContextSeq <- function(query = NULL, features, genome, translate = FALSE)
{

    if (length(query) > 1) {

        stop("query must have length 1")

    }
    
    ## subset features

    featureID <- as.integer(unlist(strsplit(as.character(query), ",")))
    i_query <- which(featureID(features) %in% featureID)
    geneID <- unique(geneID(features)[i_query])
    features <- features[geneID(features) %in% geneID]

    ## commplete frames

    if (translate) {

        inference <- "upstream"        
        features <- completeFrames(features, inference)

        i <- which(featureID(features) %in% featureID)
        f <- mcols(features)$frameStart[i][[1]]
        
        if (length(f) == 0) {

            inference <- "downstream"        
            features <- completeFrames(features, inference)

        }
        
    }

    ## obtain nucleotide sequence context

    segments <- convertToTxSegments(features)
    context <- getContextFeatures(query, segments)
    exons_us <- getExons(context$upstream, features)
    exons_qu <- getExons(context$query, features)
    exons_ds <- getExons(context$downstream, features)
    exons_us_incl <- c(exons_us, exons_qu)
    exons_context <- c(exons_us, exons_qu, exons_ds)
    dna <- getDNA(exons_context, genome)

    l <- cumsum(width(exons_context))
        
    if (length(exons_us) > 0) { start <- l[length(exons_us)] + 1
    } else { start <- 1 }
    if (length(exons_ds) > 0) { end <- l[length(exons_us_incl)]
    } else { end <- nchar(dna) }

    results <- data.frame(
        nt_seq = dna,
        nt_start = start,
        nt_end = end,
        stringsAsFactors = FALSE)

    ## translation

    if (translate) {

        frame <- getFrameStart(exons_qu, sum(width(exons_us)))
        frame_ds <- getFrameStart(exons_ds)

        if (length(frame) == 0 || length(frame_ds) == 0) {

            in_frame <- NA

        } else {

            if (all(frame == -1) && all(frame_ds == -1)) {

                in_frame <- NA

            } else {
            
                frame_next <- nextFrame(frame, sum(width(exons_us_incl)))
                in_frame <- any(setdiff(frame_next, -1) %in%
                    setdiff(frame_ds, -1))

            }

        }

        results$frame_inferred_from <- inference
        results$frame <- paste(frame, collapse = ",")
        results$frame_ds <- paste(frame_ds, collapse = ",")
        results$in_frame <- in_frame

        if (length(frame) != 1 ||
            (length(frame) == 1 && frame == -1)) {

            results$aa_seq <- NA_character_
            results$aa_start <- NA_integer_
            results$aa_end <- NA_integer_
            results$aa_start_codon <- NA_integer_
            results$aa_stop_codon <- NA_integer_
            results$nt_start_codon <- NA_integer_
            results$nt_stop_codon <- NA_integer_

            return(results)

        }

        offset <- 3 - frame
        aa <- getAA(dna, offset)
        l <- cumsum(width(exons_context)) - offset
        
        if (length(exons_us) > 0) {
          
            start <- floor(l[length(exons_us)] / 3) + 1
            
        } else {

            start <- 1

        }

        if (length(exons_ds) > 0) {

            end <- ceiling(l[length(exons_us_incl)] / 3)
            
        } else {

            end <- nchar(aa)

        }

        codons <- startStopCodons(aa, start, end)
        
        results$aa_seq <- aa
        results$aa_start <- start
        results$aa_end <- end
        results$aa_start_codon <- codons$start
        results$aa_stop_codon <- codons$stop
        results$nt_start_codon <- (codons$start - 1) * 3 + 1 + offset
        results$nt_stop_codon <- (codons$stop - 1) * 3 + 1 + offset
        
    }
    
    return(results)
    
}

completeFrames <- function(features, inference = c("upstream", "downstream"))
{

    inference <- match.arg(inference)
    
    s <- as.character(strand(unlist(range(asGRanges(features)))))
    
    i <- which(type(features) == "E" &
        elementLengths(mcols(features)$frameStart) == 0)
    i <- i[order(features[i],
        decreasing = (s == switch(inference,
            upstream = "-", downstream = "+")))]

    for (j in i) {

        features <- completeFramePerExon(j, features,
            switch(inference, upstream = TRUE, downstream = FALSE))

    }
    
    return(features)
    
}

completeFramePerExon <- function(i, features, start)
{

    Q <- features[i]
    S <- flank(Q, -1, start)
    J <- features[type(features) == "J"]
    J_prox <- flank(J, -1, !start)
    J_dist <- flank(J, -1, start)
    E <- features[type(features) == "E"]
    E_S <- flank(E, -1, !start)
    E2J <- findMatches(S, J_prox)
    J2E <- findMatches(J_dist[subjectHits(E2J)], E_S)
    
    if (length(J2E) == 0) { return(features) }
    
    if (start) {
        
        frameStart <- mcols(E)$frameEnd[subjectHits(J2E)]
        frameEnd <- nextFrame(frameStart, width(Q))
        
    } else {

        frameEnd <- mcols(E)$frameStart[subjectHits(J2E)]
        frameStart <- nextFrame(frameEnd, width(Q), TRUE)
        
    }

    x <- setdiff(paste(unlist(frameStart), unlist(frameEnd),
        sep = ","), ",")

    if (length(x) == 0) { return(features) }
    
    y <- strsplit(x, ",")

    mcols(features)$frameStart[[i]] <-
        suppressWarnings(as.integer(sapply(y, "[", 1)))
    mcols(features)$frameEnd[[i]] <-
        suppressWarnings(as.integer(sapply(y, "[", 2)))
    
    return(features)
    
}

getContextFeatures <- function(query, segments)
{

    ids_query <- as.integer(unlist(
        strsplit(as.character(query), ",", fixed = TRUE)))

    i <- grep(paste0("(^|,)", query, "(,|$)"), featureID(segments))
    from <- from(segments)[i]
    to <- to(segments)[i]
    ids_seg_query <- as.integer(unlist(
        strsplit(featureID(segments)[i], ",", fixed = TRUE)))
    
    i_segs_upstream <- integer()
    i <- which(to(segments) == from)
    
    while (length(i) == 1) {

        i_segs_upstream <- c(i, i_segs_upstream)
        i <- which(to(segments) == from(segments)[i])

    }

    i_segs_downstream <- integer()
    i <- which(from(segments) == to)
    
    while (length(i) == 1) {

        i_segs_downstream <- c(i_segs_downstream, i)
        i <- which(from(segments) == to(segments)[i])

    }
    
    ids_segs_upstream <- as.integer(unlist(strsplit(
        featureID(segments)[i_segs_upstream], ",", fixed = TRUE)))
    ids_segs_downstream <- as.integer(unlist(strsplit(
        featureID(segments)[i_segs_downstream], ",", fixed = TRUE)))

    ids_context <- c(ids_segs_upstream, ids_seg_query, ids_segs_downstream)

    i_first <- which(ids_context == head(ids_query, 1))
    i_last <- which(ids_context == tail(ids_query, 1))

    if (i_first > 1) {

        ids_upstream <- ids_context[seq_len(i_first - 1)]

    } else {

        ids_upstream <- NULL

    }

    if (i_last < length(ids_context)) {

        ids_downstream <- ids_context[(i_last + 1):length(ids_context)]

    } else {

        ids_downstream <- NULL

    }
    
    list_ids <- list(
        upstream = ids_upstream,
        query = ids_query,
        downstream = ids_downstream
    )

    return(list_ids)
    
}

getExons <- function(ids, features)
{
    
    ids_unlisted <- unlist(ids)
    i <- match(ids_unlisted, featureID(features))
    j <- which(type(features)[i] == "E")
    
    if (is(ids, "list")) {
        
        exons <- split(features[i][j], togroup(ids)[j])

    } else {

        exons <- features[i][j]

    }
    
    return(exons)
    
}

getDNA <- function(exons, genome)
{

    exons_unlisted <- unlist(exons)
    dna <- as.character(getSeq(genome, exons_unlisted))

    if (is(exons, "GRangesList")) {

        dna <- tapply(dna, togroup(exons), paste, collapse = "")

    } else {

        dna <- paste(dna, collapse = "")

    }
    
    return(dna)

}

getAA <- function(dna, offset)
{

    dna <- substr(dna, start = offset + 1, stop = nchar(dna))
    dna <- substr(dna, start = 1, stop = nchar(dna) - nchar(dna) %% 3)
    aa <- as.character(translate(DNAString(dna)))
    
    return(aa)
    
}

getFrameStart <- function(query, offset = 0)
{

    if (length(query) == 0) {

        return()

    }
    
    f <- unique(mcols(query)$frameStart[[1]])

    if (length(f) == 1 && f == -1) {

        return(f)

    }

    f <- f[f != -1]
    f <- nextFrame(f, offset, TRUE)

    return(f)
    
}

startStopCodons <- function(seq, start, end)
{

    pos_start <- gregexpr("M", seq, fixed = TRUE)[[1]]
    pos_stop <- gregexpr("*", seq, fixed = TRUE)[[1]]

    pos_start <- pos_start[pos_start >= start & pos_start <= end]
    pos_stop <- pos_stop[pos_stop >= start & pos_stop <= end]
    
    if (any(pos_start > max(pos_stop))) {
        
        pos_start <- min(pos_start[pos_start > max(pos_stop)])
        
    } else {
        
        pos_start <- NA_integer_
        
    }
    
    if (length(pos_stop) > 0) {
        
        pos_stop <- min(pos_stop)
        
    } else {
        
        pos_stop <- NA_integer_
        
    }

    out <- list(start = pos_start, stop = pos_stop)
    
    return(out)

}
