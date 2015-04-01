spliceGraph <- function(features)
{

    if (is(features, "SGFeatures")) {
        
        features <- asGRanges(features)

    } else if (is(features, "SGSegments")) {
        
        features <- asGRangesList(features)

    }

    if (is(features, "GRanges")) {
        
        option <- "features"

    } else if (is(features, "GRangesList")) {
        
        option <- "segments"

    }
    
    if (option == "features") {

        edges <- features[mcols(features)$type %in% c("J", "E")]
        d <- cbind(from = NA, to = NA, as.data.frame(mcols(edges)))

        nodes <- features[mcols(features)$type %in% c("D" ,"A")]
        nodes_pos <- gr2pos(nodes)
        nodes_name <- feature2name(nodes)

        i_D <- which(mcols(nodes)$type == "D")
        i_A <- which(mcols(nodes)$type == "A")

        ## identify 'from' nodes
        
        i <- which(mcols(edges)$type == "J")
        pos <- gr2pos(flank(edges[i], -1, start = TRUE))
        d$from[i] <- paste0("D:", pos)
        i <- which(mcols(edges)$type == "E" & mcols(edges)$splice5p)
        pos <- gr2pos(flank(edges[i], -1, start = TRUE))
        d$from[i] <- paste0("A:", pos)
        i <- which(mcols(edges)$type == "E" & !mcols(edges)$splice5p)
        pos <- gr2pos(flank(edges[i], -1, start = TRUE))
        d$from[i] <- nodes_name[i_A][match(pos, nodes_pos[i_A])]
        i <- i[which(is.na(d$from[i]))]
        pos <- gr2pos(suppressWarnings(flank(edges[i], 1, start = TRUE)))
        d$from[i] <- nodes_name[i_D][match(pos, nodes_pos[i_D])]
        i <- i[which(is.na(d$from[i]))]
        pos <- gr2pos(flank(edges[i], -1, start = TRUE))
        d$from[i] <- paste0("S:", pos)

        ## identify 'to' nodes
        
        i <- which(mcols(edges)$type == "J")
        pos <- gr2pos(flank(edges[i], -1, start = FALSE))
        d$to[i] <- paste0("A:", pos)
        i <- which(mcols(edges)$type == "E" & mcols(edges)$splice3p)
        pos <- gr2pos(flank(edges[i], -1, start = FALSE))
        d$to[i] <- paste0("D:", pos)
        i <- which(mcols(edges)$type == "E" & !mcols(edges)$splice3p)
        pos <- gr2pos(flank(edges[i], -1, start = FALSE))
        d$to[i] <- nodes_name[i_D][match(pos, nodes_pos[i_D])]
        i <- i[which(is.na(d$to[i]))]
        pos <- gr2pos(suppressWarnings(flank(edges[i], 1, start = FALSE)))
        d$to[i] <- nodes_name[i_A][match(pos, nodes_pos[i_A])]
        i <- i[which(is.na(d$to[i]))]
        pos <- gr2pos(flank(edges[i], -1, start = FALSE))
        d$to[i] <- paste0("E:", pos)

    } else {
        
        d <- as.data.frame(mcols(features))
        
    }
                    
    g <- graph.data.frame(d = d, directed = TRUE)
    gd <- edges(g)
    gv <- nodes(g)
    gv$type <- substr(gv$name, 1, 1)

    if (option == "features") {
    
        gv$featureID <- mcols(nodes)$featureID[match(gv$name, nodes_name)]

    }
    
    gv_cluster <- as.integer(clusters(g)$membership)
    gd_cluster <- gv_cluster[match(gd$from, gv$name)]
    
    if (is.null(mcols(features)$geneID)) {
            
        gv$geneID <- gv_cluster
        gd$geneID <- gd_cluster
            
    } else {

        cluster2geneID <- tapply(gd$geneID, gd_cluster, unique,
            simplify = FALSE)

        if (any(elementLengths(cluster2geneID) > 1)) {

            stop("splice graph inconsistent with geneIDs")

        } else {

            cluster2geneID <- unlist(cluster2geneID)

        }

        gv$geneID <- cluster2geneID[match(gv_cluster, names(cluster2geneID))]

    }
        
    ## For each gene reorder nodes in 5' to 3' direction 
    ## and by type (S -> A -> D -> E)

    tmp <- strsplit(gv$name, split = ":", fixed = TRUE)
    type <- c(S = 0, A = 1, D = 2, E = 3)[sapply(tmp, "[", 1)]
    pos <- as.integer(sapply(tmp, "[", 3))
    strand <- sapply(tmp, "[", 4)
    i_neg <- which(strand == "-")
    pos[i_neg] <- -1 * pos[i_neg]
    gv <- gv[order(gv$geneID, pos, type), ]    

    ## For each gene reorder edges in 5' to 3' direction 

    split_from <- strsplit(gd$from, split = ":", fixed = TRUE)
    split_to <- strsplit(gd$to, split = ":", fixed = TRUE)
    start <- as.integer(sapply(split_from, "[", 3))
    end <- as.integer(sapply(split_to, "[", 3))
    strand <- sapply(split_from, "[", 4)
    i_neg <- which(strand == "-")
    tmp_start_neg <- start[i_neg]
    tmp_end_neg <- end[i_neg]
    start[i_neg] <- -1 * tmp_end_neg
    end[i_neg] <- -1 * tmp_start_neg
    gd <- gd[order(gd$geneID, start, end), ]    

    ## Create splice graph

    g <- graph.data.frame(d = gd, directed = TRUE, vertices = gv)

    return(g)
        
}

nodes <- function(g)
{

    get.data.frame(g, "vertices")

}

edges <- function(g)
{

    get.data.frame(g, "edges")

}

neighborhood2 <- function(graph, order, nodes, mode)
{

    n <- neighborhood(graph, order, nodes, mode)
    n <- mapply(setdiff, n, nodes, SIMPLIFY = FALSE)
    
    return(n)

}

addRootAndLeafNodes <- function(g)
{

    gv <- nodes(g)
    gd <- edges(g)

    strand <- strsplit(gv$name[1], ":", fixed = TRUE)[[1]][4]
    
    ## Add unique root node and edges
    
    i <- which(gv$type == "S")

    if (length(i) > 0) {

        pos <- as.integer(sapply(strsplit(gv$name[i], ":", fixed = TRUE),
            "[", 3))
        i_R <- i[which.min(switch(strand, "+" = 1, "-" = -1) * pos)]
        R <- sub("^S", "R", gv$name[i_R])
        
        gd_R <- data.frame(matrix(NA, nrow = length(i), ncol = ncol(gd)))
        names(gd_R) <- names(gd)
        gd_R$from <- R
        gd_R$to <- gv$name[i]
        gd <- rbind(gd_R, gd)

        gv_R <- data.frame(matrix(NA, nrow = 1, ncol = ncol(gv)))
        names(gv_R) <- names(gv)
        gv_R$name <- R
        gv <- rbind(gv_R, gv)

    }

    ## Add unique leaf nodes and edges

    i <- which(gv$type == "E")

    if (length(i) > 0) {

        pos <- as.integer(sapply(strsplit(gv$name[i], ":", fixed = TRUE),
            "[", 3))
        i_L <- i[which.max(switch(strand, "+" = 1, "-" = -1) * pos)]
        L <- sub("^E", "L", gv$name[i_L])

        gd_L <- data.frame(matrix(NA, nrow = length(i), ncol = ncol(gd)))
        names(gd_L) <- names(gd)
        gd_L$from <- gv$name[i]
        gd_L$to <- L
        gd <- rbind(gd, gd_L)

        gv_L <- data.frame(matrix(NA, nrow = 1, ncol = ncol(gv)))
        names(gv_L) <- names(gv)
        gv_L$name <- L
        gv <- rbind(gv, gv_L)

    }
    
    ## Create extended splice graph
    
    g <- graph.data.frame(d = gd, directed = TRUE, vertices = gv)

    return(g)
    
}

subgraph <- function(g, geneIDs)
{

    gv <- nodes(g)
    i <- which(gv$geneID %in% geneIDs)
    g <- induced.subgraph(g, i)

    return(g)
    
}

findSegments <- function(features, cores = 1)
{

    g <- spliceGraph(features)

    gv <- nodes(g)
    gd <- edges(g)

    i_branch <- which(gv$type == "S" | gv$type == "E" |
        degree(g, mode = "out") > 1 | degree(g, mode = "in") > 1)
    geneID_n_branch <- table(gv$geneID[i_branch])

    ## Collapse graph for geneIDs with single isoform
    ## NOTE edges must be ordered in 5' to 3' direction

    geneIDs_1 <- as.integer(names(which(geneID_n_branch == 2)))
    i <- which(gd$geneID %in% geneIDs_1)
    segments_1 <- IntegerList(split(gd$featureID[i], gd$geneID[i]))

    ## Collapse graph for geneID with multiple isoforms
    ## NOTE nodes must be ordered in 5' to 3' direction

    geneIDs_2 <- as.integer(names(which(geneID_n_branch > 2)))
    list_segments_2 <- mclapply(geneIDs_2, findSegmentsPerGene, g = g,
        mc.cores = cores)
    segments_2 <- IntegerList(do.call(c, list_segments_2))

    segments <- c(segments_1, segments_2) 
        
    segmentID <- togroup(segments)[match(featureID(features),
        unlist(segments))]

    return(segmentID)

}

findSegmentsPerGene <- function(g, geneID)
{

    h <- subgraph(g, geneID)
    
    ## Extract nodes and edges
    
    hv <- nodes(h)
    hd <- edges(h)
    
    ## Find source and target nodes
    
    sources <- which(hv$type != "E" & (hv$type == "S" |
        degree(h, mode = "out") > 1 | degree(h, mode = "in") > 1))
    targets <- which(hv$type != "S" & (hv$type == "E" |
        degree(h, mode = "out") > 1 | degree(h, mode = "in") > 1))

    names(sources) <- NULL
    names(targets) <- NULL
    
    list_neighbors <- neighborhood2(h, 1, sources, "out")
    
    s <- sources[togroup(list_neighbors)]
    n <- unlist(list_neighbors)
    
    r <- is.finite(shortest.paths(h, n, targets, "out"))
    t <- targets[apply(r, 1, function(x) { min(which(x)) })]
    
    fun <- function(from, to)
    {

        if (from == to) {

            p <- integer()

        } else if (edge.connectivity(h, from, to) == 1) {

            p <- get.shortest.paths(h, from, to, output = "epath")$epath[[1]]
            p <- as.integer(p)

        } else {

            p <- which(hd$from == hv$name[from] & hd$to == hv$name[to])
            p <- as.list(p)

        }

        return(p)

    }
    
    segments_1 <- mapply(fun, s, n, SIMPLIFY = FALSE)
    segments_2 <- mapply(fun, n, t, SIMPLIFY = FALSE)

    i <- which(sapply(segments_1, is.list))

    if (length(i) > 0) {

        segments_2 <- c(segments_2[-i],
            segments_2[rep(i, elementLengths(segments_1[i]))])
        segments_1 <- c(segments_1[-i],
            unlist(segments_1[i], recursive = FALSE))

    }

    segments <- pc(IntegerList(segments_1), IntegerList(segments_2))
    segments <- relist(hd$featureID[unlist(segments)], segments)

    return(segments)
    
}

##' @title Find transcript variants from splice graph
##' @param features \code{SGFeatures} object
##' @param maxnvariant If more than \code{maxnvariant} variants are
##'   identified in an event, the gene is skipped, resulting in a warning.
##'   Set to \code{NA} to include all genes.
##' @param annotate_events Logical indicating whether identified
##'   transcript variants should be annotated in terms of canonical events.
##'   For details see help page for \code{\link{annotateSGVariants}}.
##' @param cores Number of cores available for parallel processing
##' @return A \code{SGVariants} object
##' @examples
##' sgv <- findSGVariants(sgf)
##' @author Leonard Goldstein

findSGVariants <- function(features, maxnvariant = 20, annotate_events = TRUE,
    cores = 1)
{

    if (!is(features, "SGFeatures")) {

        stop("features must be an SGFeatures object")

    } 

    variants <- findSGVariantsFromSGFeatures(features, maxnvariant, cores)

    if (length(variants) == 0) {

        warning("features do not include any alternative transcript
            events, no transcript variants identified")
        return(variants)

    }
    
    if (annotate_events) {

        message("Annotate variants...")
        variants <- annotateSGVariants(variants)

    }
    
    variantName(variants) <- makeVariantNames(variants)

    return(variants)

}

findSGVariantsFromSGFeatures <- function(features, maxnvariant, cores = 1)
{

    message("Find segments...")
    segments <- convertToSGSegments(features, cores)

    message("Find variants...")
    g <- spliceGraph(segments)
    i <- which(degree(g, mode = "out") > 1 | degree(g, mode = "in") > 1)
    geneIDs <- unique(nodes(g)$geneID[i])

    if (length(geneIDs) == 0) {

        return(SGVariants())

    }
    
    list_variant_info <- mclapply(geneIDs, findSGVariantsPerGene,
        g = g, maxnvariant = maxnvariant, mc.cores = cores)
    variant_info <- rbindListOfDFs(list_variant_info, cores)

    if (!is.na(maxnvariant) && nrow(variant_info) == 0) {

        return(SGVariants())

    }

    variant_info$eventID <- eventIDs(variant_info)
    variant_info <- variant_info[order(variant_info$eventID), ]
    variant_info$variantID <- seq_len(nrow(variant_info))
    variant_info$featureID5p <- getRepresentativeFeatureIDs(
        variant_info, features, TRUE)
    variant_info$featureID3p <- getRepresentativeFeatureIDs(
        variant_info, features, FALSE)
    
    variant_featureID <- variant_info$featureID
    variant_featureID <- gsub("(", "", variant_featureID, fixed = TRUE)
    variant_featureID <- gsub(")", "", variant_featureID, fixed = TRUE)
    variant_featureID <- gsub("|", ",", variant_featureID, fixed = TRUE)
    variant_featureID <- strsplit(variant_featureID, ",", fixed = TRUE)

    variants <- split(features[match(unlist(variant_featureID),
        featureID(features))], togroup(variant_featureID))
    mcols(variants) <- variant_info
    variants <- SGVariants(variants)
    variants <- annotatePaths(variants)

    return(variants)

}

eventIDs <- function(variant_info)
{

    from_to <- paste(variant_info$from, variant_info$to)

    unique_from_to <- unique(from_to)
    unique_from <- sapply(strsplit(unique_from_to, " ", fixed = TRUE), "[", 1)
    unique_to <- sapply(strsplit(unique_from_to, " ", fixed = TRUE), "[", 2)
    
    sl <- sapply(strsplit(unique_from, ":", fixed = TRUE), "[", 2)
    st <- sapply(strsplit(unique_from, ":", fixed = TRUE), "[", 4)
    start <- as.integer(sapply(strsplit(unique_from, ":", fixed = TRUE),
        "[", 3))
    end <- as.integer(sapply(strsplit(unique_to, ":", fixed = TRUE), "[", 3))
    tmp_start <- start
    tmp_end <- end
    i_neg <- which(st == "-")
    start[i_neg] <- -tmp_end[i_neg]
    end[i_neg] <- -tmp_start[i_neg]    
    o <- order(sl, start, end)
    
    eventIDs <- match(from_to, unique_from_to[o])

    return(eventIDs)
    
}

getRepresentativeFeatureIDs <- function(variant_info, features, start = TRUE)
{

    if (start) {

        variant_node <- variant_info$from
        variant_informative <- variant_info$closed3p

    } else {

        variant_node <- variant_info$to
        variant_informative <- variant_info$closed5p

    }

    variant_rep_id <- getTerminalFeatureIDs(variant_info$featureID, start)

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

findSGVariantsPerGene <- function(g, geneID, maxnvariant)
{

    ## Extract subgraph corresponding to geneID
    h <- subgraph(g, geneID)
        
    ## Add unique root and leaf nodes    
    h <- addRootAndLeafNodes(h)

    ## Extract data frames of nodes and edges    
    hv <- nodes(h)
    hd <- edges(h)

    ## Find events
    b <- findEvents(h)

    ## Sort events
    b <- sortEvents(b, h)

    ## Initialize list of alternative paths
    list_path_info <- vector("list", nrow(b))
    
    ## Initialize data frame of recursively defined paths
    ref <- hd[, c("from", "to", "type", "featureID", "segmentID")]
    ref$segmentID <- as.character(ref$segmentID)
    
    for (k in seq_len(nrow(b))) {

        from <- hv$name[b$source[k]]
        to <- hv$name[b$target[k]]

        paths_index_ref <- findAllPaths(from, to, NULL, ref, hv$name)

        i <- unlist(paths_index_ref)
        f <- togroup(paths_index_ref)

        paths_type <- unstrsplit(split(ref$type[i], f), "")
        paths_featureID <- unstrsplit(split(ref$featureID[i], f), ",")
        paths_segmentID <- unstrsplit(split(ref$segmentID[i], f), ",")

        if (!is.na(maxnvariant) && length(paths_index_ref) > maxnvariant) {

            warning(paste("number of variants exceeds maxnvariant in gene",
                geneID), call. = FALSE, immediate. = TRUE)
            return()

        }

        o <- order(nchar(paths_type), paths_type)

        paths_type <- paths_type[o]
        paths_featureID <- paths_featureID[o]
        paths_segmentID <- paths_segmentID[o]

        ## nodes within event (including source and target nodes)
        bv <- intersect(subcomponent(h, b$source[k], "out"),
            subcomponent(h, b$target[k], "in"))

        ## nodes outside of event
        ov <- setdiff(seq_len(nrow(hv)), bv)

        ## nodes within event, excluding source and target nodes
        bv <- setdiff(bv, c(b$source[k], b$target[k]))

        ## find adjacent nodes for all nodes within event
        bn_out <- unique(unlist(neighborhood(h, 1, bv, "out")))
        bn_in <- unique(unlist(neighborhood(h, 1, bv, "in")))

        closed3p <- length(intersect(bn_out, ov)) == 0
        closed5p <- length(intersect(bn_in, ov)) == 0
        
        if (k < nrow(b) && closed3p && closed5p) {
          
            ## Include event in data frame of recursively defined paths
            ref <- rbind(ref, data.frame(
                from = from,
                to = to,
                type = paste0("(", paste(paths_type, collapse = "|"), ")"),
                featureID = paste0("(", paste(paths_featureID,
                    collapse = "|"), ")"),
                segmentID = paste0("(", paste(paths_segmentID,
                    collapse = "|"), ")"),
                stringsAsFactors = FALSE))

            ref <- ref[-unique(unlist(paths_index_ref)), ]

        }
        
        list_path_info[[k]] <- data.frame(
            from = rep(from, length(o)),
            to = rep(to, length(o)),
            type = paths_type,
            featureID = paths_featureID,
            segmentID = paths_segmentID,
            closed3p = closed3p,
            closed5p = closed5p,
            stringsAsFactors = FALSE)
         
    }
    
    path_info <- do.call(rbind, list_path_info)

    ## remove references to artifical root and leaf nodes
    path_info$type <- gsub("NA", "", path_info$type, fixed = TRUE)
    path_info$featureID <- gsub("NA,", "", path_info$featureID, fixed = TRUE)
    path_info$featureID <- gsub(",NA", "", path_info$featureID, fixed = TRUE)
    path_info$segmentID <- gsub("NA,", "", path_info$segmentID, fixed = TRUE)
    path_info$segmentID <- gsub(",NA", "", path_info$segmentID, fixed = TRUE)

    ## include geneID
    path_info$geneID <- geneID

    return(path_info)
    
}

findEvents <- function(g)
{

    gv <- nodes(g)

    i_source <- setNames(which(degree(g, mode = "out") > 1), NULL)
    
    list_source <- vector()
    list_target <- vector()
    
    for (s in i_source) {

        ## Initialize vector of target nodes completing events from s
        t <- vector()
        
        ## Consider alternative nodes immediately downstream (proximal) from s
        alt_prox <- neighborhood2(g, 1, s, "out")[[1]]

        ## Proximal nodes with connectivity > 1 complete events
        i <- which(sapply(alt_prox, edge.connectivity, graph = g, source = s)
            > 1)
        t <- c(t, alt_prox[i])        
        
        ## Consider branchpoints reachable from s (excluding s)
        branchpts <- which(is.finite(shortest.paths(g, s, mode = "out")))[-1]

        ## Keep track of alternative paths in terms of both proximal and
        ## downstream (distal) nodes when proceeding through branchpoints
        alt_dist <- alt_prox

        j <- 1
        reachable <- is.finite(shortest.paths(g, alt_dist, branchpts[j],
            "out"))

        while (!all(reachable)) {

            i <- which(reachable)

            prox <- unique(alt_prox[i])
            dist <- neighborhood2(g, 1, branchpts[j], "out")[[1]]
            
            if (length(prox) > 1) {

                t <- c(t, branchpts[j])
                prox <- branchpts[j]

            }

            alt_prox <- alt_prox[-i]
            alt_dist <- alt_dist[-i]

            n_prox <- length(prox)
            n_dist <- length(dist)
            
            alt_prox <- c(alt_prox, rep(prox, rep(n_dist, n_prox)))
            alt_dist <- c(alt_dist, rep(dist, n_prox))
            
            j <- j + 1
            reachable <- is.finite(shortest.paths(g, alt_dist,
                branchpts[j], "out"))

        }

        t <- c(t, branchpts[j])
        t <- unique(t)
        
        list_source <- c(list_source, rep(s, length(t)))
        list_target <- c(list_target, t)

    }
    
    ## exclude special case when source and target node are
    ## root and leaf node, respectively (cf. PLEKHA5)

    exclude <- intersect(grep("^R", gv$name[list_source]),
        grep("^L", gv$name[list_target]))
    
    if (length(exclude) > 0) {
        
        list_source <- list_source[-exclude]
        list_target <- list_target[-exclude]
        
    }
    
    events <- data.frame(source = list_source, target = list_target)
    
    return(events)
    
}

sortEvents <- function(events, g)
{

    name_split <- strsplit(nodes(g)$name, ":", fixed = TRUE)
    gv_type <- sapply(name_split, "[", 1)
    gv_pos <- as.integer(sapply(name_split, "[", 3))
    st <- setdiff(sapply(name_split, "[", 4), NA)

    source_type <- gv_type[events$source]
    source_pos <- gv_pos[events$source]
    
    target_type <- gv_type[events$target]
    target_pos <- gv_pos[events$target]

    i_I <- which(source_type != "R" & target_type != "L")
    i_I <- i_I[order(target_pos[i_I] - source_pos[i_I],
        decreasing = (st == "-"))]

    i_R <- which(source_type == "R")
    i_R <- i_R[order(target_pos[i_R], decreasing = (st == "-"))]

    i_L <- which(target_type == "L")
    i_L <- i_L[order(source_pos[i_L], decreasing = (st == "+"))]
    
    i <- c(i_I, i_R, i_L)

    events <- events[i, ]
    
    return(events)
    
}

findAllPaths <- function(from, to, path, ref, nodes)
{

    if (from == to && !is.null(path)) {

        return(path)
        
    } else {

        i <- which(ref$from == from)
        i <- i[match(ref$to[i], nodes) <= match(to, nodes)]

        if (length(i) > 0) {
        
            paths <- lapply(i, append, path, 0)
            
            paths <- mapply(findAllPaths,
                from = ref$to[i],
                path = paths,
                MoreArgs = list(to = to, ref = ref, nodes = nodes),
                SIMPLIFY = FALSE,
                USE.NAMES = FALSE)
            
            i <- which(sapply(paths, is.list))

            if (length(i) > 0) {
                
                paths <- c(paths[-i], unlist(paths[i], recursive = FALSE))

            }
            
            i <- which(elementLengths(paths) == 0)
            
            if (length(i) > 0) { paths <- paths[-i] }
            
            return(paths)

        }
        
    }
        
}

##' Annotate transcript variants in terms of canonical events.
##'
##' The following events are considered:
##' \dQuote{SE} (skipped exon),
##' \dQuote{S2E} (two consecutive exons skipped),
##' \dQuote{RI} (retained intron),
##' \dQuote{MXE} (mutually exclusive exons),
##' \dQuote{A5SS} (alternative 5' splice site),
##' \dQuote{A3SS} (alternative 3' splice site),
##' \dQuote{AFE} (alternative first exon),
##' \dQuote{ALE} (alternative last exon),
##' \dQuote{AS} (alternative start other than \dQuote{AFE}) and
##' \dQuote{AE} (alternative end other than \dQuote{ALE}).
##'
##' These are binary events, defined by two alternative variants.
##' A variant is annotated as a canonical event if it coincides with one
##' of the two variants in the canonical event, and there is at least one
##' variant in the same event that coincides with the second variant of the
##' canonical event.
##'
##' @title Annotate transcript variants in terms of canonical events
##' @param variants \code{SGVariants} object
##' @return \code{variants} with added elementMetadata column
##'   \dQuote{variantType} indicating canonical event(s)
##' @keywords internal
##' @author Leonard Goldstein

annotateSGVariants <- function(variants)
{

    ## asymmetric events
  
    list_ae_event_1 <- c("SE:I", "S2E:I", "RI:E", "A5SS:P", "A3SS:P")
    list_ae_event_2 <- c("SE:S", "S2E:S", "RI:R", "A5SS:D", "A3SS:D")

    list_ae_type_1 <- c("^JE+J$", "^JE+JE+J$", "^J$", "^E+J$", "^JE+$")
    list_ae_type_2 <- c("^J$", "^J$", "^E+$", "^J$", "^J$")

    ## symmetric events
    
    list_se_event <- c("AFE", "ALE", "MXE")
    list_se_type <- c("^E+J$", "^JE+$", "^JE+J$")

    ## maximum number of Js in type patterns

    t <- c(list_ae_type_1, list_ae_type_2, list_se_type)
    max_n_J <- max(elementLengths(gregexpr("J", t)))
    
    ## preliminaries
    
    path_from_to <- paste(from(variants), to(variants))
    path_event <- CharacterList(vector("list", length(variants)))
    path_type <- expandType(type(variants), max_n_J)
    
    ## asymmetric events

    for (k in seq_along(list_ae_type_1)) {

        event_1 <- list_ae_event_1[k]
        event_2 <- list_ae_event_2[k]

        type_1 <- list_ae_type_1[k]
        type_2 <- list_ae_type_2[k]

        i_1 <- unique(togroup(path_type)[grep(type_1, unlist(path_type))])
        i_2 <- unique(togroup(path_type)[grep(type_2, unlist(path_type))])
        from_to <- intersect(path_from_to[i_1], path_from_to[i_2])
        i_1 <- i_1[path_from_to[i_1] %in% from_to]
        i_2 <- i_2[path_from_to[i_2] %in% from_to]

        path_event[i_1] <- pc(path_event[i_1], rep(event_1, length(i_1)))
        path_event[i_2] <- pc(path_event[i_2], rep(event_2, length(i_2)))

    }

    ## symmetric events
    
    for (k in seq_along(list_se_type)) {

        event <- list_se_event[k]
        type <- list_se_type[k]

        i <- unique(togroup(path_type)[grep(type, unlist(path_type))])

        if (event == "AFE") { i <- i[grep("^R", from(variants)[i])] }
        if (event == "ALE") { i <- i[grep("^L", to(variants)[i])] }
        
        if (length(i) > 0) {

            from_to_n <- table(path_from_to[i])
            from_to <- names(which(from_to_n >= 2))
            i <- i[path_from_to[i] %in% from_to]
            path_event[i] <- pc(path_event[i], rep(event, length(i)))
            
        }
        
    }

    ## alternative start/end

    i <- intersect(grep("^R", from(variants)),
        which(!any(path_event == "AFE")))
    path_event[i] <- pc(path_event[i], rep("AS", length(i)))

    i <- intersect(grep("^L", to(variants)),
        which(!any(path_event == "ALE")))
    path_event[i] <- pc(path_event[i], rep("AE", length(i)))

    variantType(variants) <- path_event

    return(variants)

}

expandString <- function(x, return_full = FALSE)
{

    if (length(grep("(\\[|\\]|:)", x)) > 0) {

        stop("x contains characters '[', ']' or ':'")
      
    }
  
    x_f <- seq_along(x)
    x_d <- as(x, "CompressedCharacterList")    

    i <- grep("(", x, fixed = TRUE)
    
    while (length(i) > 0) {

        z <- maskInnerEvents(x[i])
      
        ## find event
        m <- regexpr("\\(\\S+\\)", z)
        l <- attr(m, "match.length")
        
        ## split at event
        u <- substr(z, 1, m - 1)
        b <- substr(z, m + 1, m + l - 2)
        v <- substr(z, m + l, nchar(z))

        ## expand event
        b <- strsplit(b, "|", fixed = TRUE)
        n <- elementLengths(b)
        y <- paste0(rep(u, n), unlist(b), rep(v, n))
        y <- unmaskEvents(y)
        y_f <- x_f[i][togroup(b)]
        y_d <- pc(x_d[i][togroup(b)], unlist(b))
        y_d <- as(y_d, "CompressedCharacterList")
        
        ## update x, x_f
        x <- c(x[-i], y)
        x_f <- c(x_f[-i], y_f)
        x_d <- c(x_d[-i], y_d)

        i <- grep("(", x, fixed = TRUE)
        
    }

    if (return_full) {

        out <- DataFrame(f = x_f, x = x,
            d = as(x_d, "CompressedCharacterList"))

    } else {
    
        out <- split(x, x_f)

    }
    
    return(out)
    
}

expandType <- function(x, max_n_J = NA)
{

    x_f <- seq_along(x)    

    i <- grep("(", x, fixed = TRUE)
    
    while (length(i) > 0) {

        z <- maskInnerEvents(x[i])
      
        ## find event
        m <- regexpr("\\(\\S+\\)", z)
        l <- attr(m, "match.length")
        
        ## split at event
        u <- substr(z, 1, m - 1)
        b <- substr(z, m + 1, m + l - 2)
        v <- substr(z, m + l, nchar(z))

        if (!is.na(max_n_J)) {

            u2 <- sub("\\[\\S+$", "", u)
            v2 <- sub("^\\S+\\]", "", v)
            min_n_J <- elementLengths(gregexpr("J", paste0(u2, v2)))
            excl <- which(min_n_J > max_n_J)
            
            if (length(excl) > 0) {

                x[i][excl] <- ""
                i <- i[-excl]
                u <- u[-excl]
                b <- b[-excl]
                v <- v[-excl]
                
            }
            
        }

        if (length(i) > 0) {
          
            ## expand event
            b <- strsplit(b, "|", fixed = TRUE)
            n <- elementLengths(b)
            y <- paste0(rep(u, n), unlist(b), rep(v, n))
            y <- unmaskEvents(y)
            y_f <- x_f[i][togroup(b)]
            
            ## update x, x_f
            x <- c(x[-i], y)
            x_f <- c(x_f[-i], y_f)

        }

        i <- grep("(", x, fixed = TRUE)
        
    }

    out <- split(x, x_f)
    
    return(out)
    
}

getTerminalFeatureIDs <- function(x, start = TRUE)
{

    x_f <- seq_along(x)

    i <- grep("(", x, fixed = TRUE)
    
    while (length(i) > 0) {
        
        z <- maskInnerEvents(x[i])
      
        ## find event
        m <- regexpr("\\(\\S+\\)", z)
        l <- attr(m, "match.length")
      
        ## split at event
        u <- substr(z, 1, m - 1)
        b <- substr(z, m + 1, m + l - 2)
        v <- substr(z, m + l, nchar(z))

        if (start) {

            terminal <- u

        } else {

            terminal <- v

        }

        done <- which(nchar(terminal) > 0)

        if (length(done) > 0) {

            x[i][done] <- unmaskEvents(terminal[done])
            i <- i[-done]
            u <- u[-done]
            b <- b[-done]
            v <- v[-done]

        }

        if (length(i) > 0) {
        
            ## expand event
            b <- strsplit(b, "|", fixed = TRUE)
            n <- elementLengths(b)
            y <- paste0(rep(u, n), unlist(b), rep(v, n))
            y <- unmaskEvents(y)
            y_f <- x_f[i][togroup(b)]
            
            ## update x, x_f
            x <- c(x[-i], y)
            x_f <- c(x_f[-i], y_f)

        }

        i <- grep("(", x, fixed = TRUE)
        
    }

    if (start) {

        x <- sub(",\\S*$", "", x)

    } else {

        x <- sub("^\\S*,", "", x)

    }

    x <- suppressWarnings(as.integer(x))
    out <- as(tapply(x, x_f, setdiff, NA, simplify = FALSE), "list")
    
    return(out)
    
}

maskInnerEvents <- function(x)
{

    pattern <- "\\([^\\(\\)]+\\)"
    m <- regexpr(pattern, x)
    done <- rep(FALSE, length(x))
    done[m == -1] <- TRUE
    
    while (!all(done)) {

        i <- which(!done)
        m <- regexpr(pattern, x[i])
        l <- attr(m, "match.length")

        p1 <- m
        p2 <- m + l - 1
        
        substr(x[i], p1, p1) <- "["
        substr(x[i], p2, p2) <- "]"

        m2 <- regexpr(pattern, x[i])
    
        i_done <- which(m2 == -1)
        i_mask <- which(m2 != -1)

        if (length(i_done) > 0) {

            done[i][i_done] <- TRUE
            substr(x[i][i_done], p1[i_done], p1[i_done]) <- "("
            substr(x[i][i_done], p2[i_done], p2[i_done]) <- ")"

        }

        if (length(i_mask) > 0) {

            substr(x[i][i_mask], p1[i_mask], p2[i_mask]) <- gsub("|", ":",
                substr(x[i][i_mask], p1[i_mask], p2[i_mask]), fixed = TRUE)
      
        }
    
    }

    return(x)
  
}

unmaskEvents <- function(x)
{

    x <- gsub("[", "(", x, fixed = TRUE)
    x <- gsub("]", ")", x, fixed = TRUE)
    x <- gsub(":", "|", x, fixed = TRUE)

    return(x)
  
}

expandSGVariantCounts <- function(object, eventID = NULL, cores = 1)
{
    
    variants <- rowRanges(object)
    features <- unlist(variants)

    if (is.null(eventID)) {

        variants_selected <- variants

    } else {

        variants_selected <- variants[eventID(variants) %in% eventID]

    }
    
    expanded <- expandString(featureID(variants_selected), TRUE)
    expandedIDs <- strsplit(expanded$x, ",", fixed = TRUE)    
    expanded_i <- relist(match(unlist(expandedIDs),
        featureID(features)), expandedIDs)
    
    paths_expanded <- variants_selected[expanded$f]
    type(paths_expanded) <- unlist(expandString(type(variants_selected)))
    featureID(paths_expanded) <- expanded$x

    f <- togroup(expanded$d)
    i <- match(unlist(expanded$d), featureID(variants))

    X <- variantFreq(object)[i, , drop = FALSE]
    X <- do.call(cbind, mclapply(seq_len(ncol(X)),
        function(j) { tapply(X[, j], f, prod) }, mc.cores = cores))

    rd <- split(features[unlist(expanded_i)], togroup(expanded_i))
    mcols(rd) <- mcols(paths_expanded)
    rd <- SGVariants(rd)
    
    object_expanded <- SummarizedExperiment(assays = list("variantFreq" = X),
        rowRanges = rd, colData = colData(object))

    colnames(object_expanded) <- colnames(object)
    rownames(object_expanded) <- NULL

    return(object_expanded)
    
}

##' This function creates interpretable transcript variant names
##' taking the format GENE_EVENT_VARIANT/ORDER_TYPE.
##' GENE is based on geneName if available, and geneID otherwise.
##' EVENT and VARIANT enumerate events and variants for the same gene
##' and event, respectively. ORDER indicates the total number of
##' variants in the same event (e.g. 1/2 refers to the first out of two
##' transcript variants in the event). TYPE is based on variantType.
##' @title Create interpretable transcript variant names 
##' @param variants \code{SGVariants} object
##' @return Character vector with transcript variant names
##' @examples
##' makeVariantNames(sgv)
##' @author Leonard Goldstein

makeVariantNames <- function(variants)
{

    if (all(elementLengths(geneName(variants)) == 0)) {
        
        GENE <- geneID(variants)

    } else {

        GENE <- unstrsplit(geneName(variants), ",")
        GENE[GENE == ""] <- "NA"

    }

    original_eventID <- as.character(eventID(variants))    
    tmp <- unique(cbind(GENE, original_eventID))
    tmp_index <- reindex(tmp[, 1])
    EVENT <- tmp_index[match(original_eventID, tmp[, 2])]
    VARIANT <- reindex(original_eventID)
    eventID_n <- table(original_eventID)
    ORDER <- eventID_n[match(original_eventID, names(eventID_n))]

    if (all(elementLengths(variantType(variants)) == 0)) {

        TYPE <- NA

    } else {
    
        TYPE <- unstrsplit(variantType(variants), ",")
        TYPE <- gsub(":\\S", "", TYPE)
        TYPE[TYPE == ""] <- "OTHER"

    }
    
    variantName <- paste(GENE, EVENT,
        paste0(VARIANT, "/", ORDER), TYPE, sep = "_")

    return(variantName)
    
}

reindex <- function(f)
{

    f <- as.character(f)
    f_n <- table(f)
    f_n <- f_n[order(names(f_n))]
    o <- order(f)
    i <- as.integer(IRanges(1, f_n))
    i[match(seq_along(f), o)]
        
}
