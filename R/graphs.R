spliceGraph <- function(features)
{

    if (is(features, "SGFeatures")) {

        features <- asGRanges(features)

    } else if (is(features, "SGSegments")) {
    
        features <- asGRangesList(features)

    }
  
    if (is(features, "GRanges")) {

        edges <- features[mcols(features)$type %in% c("J", "E")]
        d <- cbind(from = NA, to = NA, as.data.frame(mcols(edges)))

        nodes <- features[mcols(features)$type %in% c("D" ,"A")]
        nodes_pos <- gr2pos(nodes)
        nodes_name <- feature2name(nodes)

        i_D <- which(mcols(nodes)$type == "D")
        i_A <- which(mcols(nodes)$type == "A")

        ## determine 'from' nodes

        ## splice junctions -> D
        i <- which(mcols(edges)$type == "J")
        pos <- gr2pos(flank(edges[i], -1, start = TRUE))
        d$from[i] <- paste0("D:", pos)
        
        ## exons with 5' splice -> A
        i <- which(mcols(edges)$type == "E" & mcols(edges)$splice5p)
        pos <- gr2pos(flank(edges[i], -1, start = TRUE))
        d$from[i] <- paste0("A:", pos)

        ## exons without 5' splice -> A (if included)
        i <- which(mcols(edges)$type == "E" & !mcols(edges)$splice5p)
        pos <- gr2pos(flank(edges[i], -1, start = TRUE))
        d$from[i] <- nodes_name[i_A][match(pos, nodes_pos[i_A])]

        ## remaining exons without 5' splice -> D (if included)
        i <- i[which(is.na(d$from[i]))]
        pos <- gr2pos(suppressWarnings(flank(edges[i], 1, start = TRUE)))
        d$from[i] <- nodes_name[i_D][match(pos, nodes_pos[i_D])]

        ## remaining edges -> S
        i <- i[which(is.na(d$from[i]))]
        pos <- gr2pos(flank(edges[i], -1, start = TRUE))
        d$from[i] <- paste0("S:", pos)

        ## determine 'to' nodes

        ## splice junctions -> A
        i <- which(mcols(edges)$type == "J")
        pos <- gr2pos(flank(edges[i], -1, start = FALSE))
        d$to[i] <- paste0("A:", pos)

        ## exons with 3' splice -> D
        i <- which(mcols(edges)$type == "E" & mcols(edges)$splice3p)
        pos <- gr2pos(flank(edges[i], -1, start = FALSE))
        d$to[i] <- paste0("D:", pos)

        ## exons without 3' splice -> D (if included)
        i <- which(mcols(edges)$type == "E" & !mcols(edges)$splice3p)
        pos <- gr2pos(flank(edges[i], -1, start = FALSE))
        d$to[i] <- nodes_name[i_D][match(pos, nodes_pos[i_D])]

        ## remaining exons without 3' splice -> A (if included)
        i <- i[which(is.na(d$to[i]))]
        pos <- gr2pos(suppressWarnings(flank(edges[i], 1, start = FALSE)))
        d$to[i] <- nodes_name[i_A][match(pos, nodes_pos[i_A])]

        ## remaining edges -> E
        i <- i[which(is.na(d$to[i]))]
        pos <- gr2pos(flank(edges[i], -1, start = FALSE))
        d$to[i] <- paste0("E:", pos)

    } else if (is(features, "GRangesList")) {
        
        d <- as.data.frame(mcols(features))
        
    }
                    
    g <- graph.data.frame(d = d, directed = TRUE)
    gd <- edges(g)
    gv <- nodes(g)
    gv$type <- substr(gv$name, 1, 1)

    if (is(features, "GRanges")) {
    
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
        
    ## For each gene, reorder nodes 5' to 3' and by type (S -> A -> D -> E)

    tmp <- strsplit(gv$name, split = ":", fixed = TRUE)
    type <- c(S = 0, A = 1, D = 2, E = 3)[sapply(tmp, "[", 1)]
    pos <- as.integer(sapply(tmp, "[", 3))
    strand <- sapply(tmp, "[", 4)
    i_neg <- which(strand == "-")
    pos[i_neg] <- -1 * pos[i_neg]
    gv <- gv[order(gv$geneID, pos, type), ]    

    ## For each gene, reorder edges in 5' to 3' direction 

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

addSourceAndSinkNodes <- function(g)
{

    gv <- nodes(g)
    gd <- edges(g)

    strand <- strsplit(gv$name[1], ":", fixed = TRUE)[[1]][4]
    
    ## Add unique source node and edges
    
    i <- which(gv$type == "S")

    if (length(i) > 0) {

        gd_R <- data.frame(matrix(NA, nrow = length(i), ncol = ncol(gd)))
        names(gd_R) <- names(gd)
        gd_R$from <- "R"
        gd_R$to <- gv$name[i]
        gd <- rbind(gd_R, gd)

        gv_R <- data.frame(matrix(NA, nrow = 1, ncol = ncol(gv)))
        names(gv_R) <- names(gv)
        gv_R$name <- "R"
        gv <- rbind(gv_R, gv)

    }

    ## Add unique sink nodes and edges

    i <- which(gv$type == "E")

    if (length(i) > 0) {

        gd_K <- data.frame(matrix(NA, nrow = length(i), ncol = ncol(gd)))
        names(gd_K) <- names(gd)
        gd_K$from <- gv$name[i]
        gd_K$to <- "K"
        gd <- rbind(gd, gd_K)

        gv_K <- data.frame(matrix(NA, nrow = 1, ncol = ncol(gv)))
        names(gv_K) <- names(gv)
        gv_K$name <- "K"
        gv <- rbind(gv, gv_K)

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

findSGSegments <- function(features, cores = 1)
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
    
    list_segments_2 <- mclapply(
        geneIDs_2,
        findSGSegmentsPerGene,
        g = g,
        mc.cores = cores)

    checkApplyResultsForErrors(
        list_segments_2,
        "findSGSegmentsPerGene",
        geneIDs_2)

    segments_2 <- IntegerList(do.call(c, list_segments_2))

    segments <- c(segments_1, segments_2) 
        
    segmentID <- togroup(segments)[match(featureID(features),
        unlist(segments))]

    return(segmentID)

}

findSGSegmentsPerGene <- function(g, geneID)
{

    h <- subgraph(g, geneID)
    
    ## Extract nodes and edges
    
    hv <- nodes(h)
    hd <- edges(h)
    
    ## Find start and end nodes
    
    starts <- which(hv$type != "E" & (hv$type == "S" |
        degree(h, mode = "out") > 1 | degree(h, mode = "in") > 1))
    ends <- which(hv$type != "S" & (hv$type == "E" |
        degree(h, mode = "out") > 1 | degree(h, mode = "in") > 1))

    names(starts) <- NULL
    names(ends) <- NULL
    
    list_neighbors <- neighborhood2(h, 1, starts, "out")
    
    s <- starts[togroup(list_neighbors)]
    n <- unlist(list_neighbors)
    
    r <- is.finite(shortest.paths(h, n, ends, "out"))
    t <- ends[apply(r, 1, function(x) { min(which(x)) })]
    
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

##' @title Find splice variants from splice graph
##' @param features \code{SGFeatures} object
##' @param maxnvariant If more than \code{maxnvariant} variants are
##'   identified in an event, the gene is skipped, resulting in a warning.
##'   Set to \code{NA} to include all genes.
##' @param annotate_events Logical indicating whether identified
##'   splice variants should be annotated in terms of canonical events.
##'   For details see help page for \code{\link{annotateSGVariants}}.
##' @param cores Number of cores available for parallel processing
##' @return An \code{SGVariants} object
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

        warning("features do not include any splice events,
            no splice variants identified")
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
    
    list_variant_info <- mclapply(
        geneIDs,
        findVariantsPerGene,
        g = g,
        maxnvariant = maxnvariant,
        mc.cores = cores)

    checkApplyResultsForErrors(
        list_variant_info,
        "findVariantsPerGene",
        geneIDs)

    variant_info <- rbindListOfDFs(list_variant_info, cores)

    if (!is.na(maxnvariant) && nrow(variant_info) == 0) {

        return(SGVariants())

    }

    ## obtain variants in terms of SGFeatures
    variant_featureID <- variant_info$featureID
    variant_featureID <- gsub("(", "", variant_featureID, fixed = TRUE)
    variant_featureID <- gsub(")", "", variant_featureID, fixed = TRUE)
    variant_featureID <- gsub("|", ",", variant_featureID, fixed = TRUE)
    variant_featureID <- strsplit(variant_featureID, ",", fixed = TRUE)
    variant_featureID <- as(variant_featureID, "CompressedIntegerList")
    variant_featureID <- unique(variant_featureID)
    variants <- split(features[match(unlist(variant_featureID),
        featureID(features))], togroup(variant_featureID))
    variant_gr <- unlist(range(variants))

    ## sort by genomic position
    i <- order(variant_gr)
    variant_info <- variant_info[i, ]
    variants <- variants[i]
    variant_gr <- variant_gr[i]
    
    ## obtain event IDs and variant IDs
    ft <- paste(variant_info$from, variant_info$to)
    variant_info$eventID <- as.integer(factor(ft, levels = unique(ft)))
    variant_info$variantID <- seq_len(nrow(variant_info))
    
    ## replace source nodes with corresponding start nodes
    x <- flank(variant_gr, -1, TRUE)
    i <- which(variant_info$from == "R")
    if (length(i) > 0) { variant_info$from[i] <- paste0("S:", gr2pos(x[i])) }

    ## replace sink nodes with corresponding end nodes
    x <- flank(variant_gr, -1, FALSE)
    i <- which(variant_info$to == "K")
    if (length(i) > 0) { variant_info$to[i] <- paste0("E:", gr2pos(x[i])) }

    ## obtain representative feature IDs
    variant_info$featureID5p <- getRepresentativeFeatureIDs(
        variant_info, features, TRUE)
    variant_info$featureID3p <- getRepresentativeFeatureIDs(
        variant_info, features, FALSE)

    ## create SGVariants object
    mcols(variants) <- variant_info
    variants <- SGVariants(variants)
    variants <- annotatePaths(variants)

    return(variants)

}

getRepresentativeFeatureIDs <- function(variant_info, features, start = TRUE)
{

    variant_rep_id <- vector("list", nrow(variant_info))
    variant_rep_id <- as(variant_rep_id, "CompressedIntegerList")
    
    if (start) {

        index <- grep("^D", variant_info$from)
        tmp_info <- variant_info[index, ]
        tmp_node <- tmp_info$from
        tmp_informative <- tmp_info$closed3p
        
    } else {

        index <- grep("^A", variant_info$to)
        tmp_info <- variant_info[index, ]
        tmp_node <- tmp_info$to
        tmp_informative <- tmp_info$closed5p
        
    }

    if (length(index) == 0) {

        return(variant_rep_id)

    }
    
    tmp_rep_id <- getTerminalFeatureIDs(tmp_info$featureID, start)

    ## replace exons with splice sites

    tmp_rep_id_unlisted <- unlist(tmp_rep_id)
    tmp_rep_id_unlisted_i <- match(tmp_rep_id_unlisted, featureID(features))

    i_E <- which(type(features)[tmp_rep_id_unlisted_i] == "E")

    if (length(i_E) > 0) {
    
        tmp_rep_id_unlisted_i[i_E] <- match(
            tmp_node[togroup(tmp_rep_id)][i_E], feature2name(features))

    }
    
    tmp_rep_id_unlisted <- featureID(features)[tmp_rep_id_unlisted_i]
    tmp_rep_id <- relist(tmp_rep_id_unlisted, tmp_rep_id)
    
    ## exclude variants due to open events or ambiguous features

    event_dup <- tapply(unlist(tmp_rep_id),
        tmp_info$eventID[togroup(tmp_rep_id)],
        function(x) { any(duplicated(x)) })
    i <- which(!tmp_informative |
        tmp_info$eventID %in% names(which(event_dup)))
    tmp_rep_id[i] <- vector("list", length(i))    

    variant_rep_id[index] <- tmp_rep_id
    
    return(variant_rep_id)
    
}

findVariantsPerGene <- function(g, geneID, maxnvariant)
{

    ## Extract subgraph corresponding to geneID
    h <- subgraph(g, geneID)
        
    ## Add unique source and sink nodes    
    h <- addSourceAndSinkNodes(h)

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

        from <- hv$name[b$start[k]]
        to <- hv$name[b$end[k]]

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

        ## nodes within event (including start and end nodes)
        bv <- intersect(subcomponent(h, b$start[k], "out"),
            subcomponent(h, b$end[k], "in"))

        ## nodes outside of event
        ov <- setdiff(seq_len(nrow(hv)), bv)

        ## nodes within event, excluding start and end nodes
        bv <- setdiff(bv, c(b$start[k], b$end[k]))

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

    ## remove references to artifical source and sink nodes
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

    i_start <- setNames(which(degree(g, mode = "out") > 1), NULL)
    
    list_start <- vector()
    list_end <- vector()
    
    for (s in i_start) {

        ## Initialize vector of end nodes completing events from s
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
        
        list_start <- c(list_start, rep(s, length(t)))
        list_end <- c(list_end, t)

    }
    
    ## exclude special case when start and end node are
    ## source and sink node, respectively (cf. PLEKHA5)

    exclude <- intersect(which(gv$name[list_start] == "R"),
        which(gv$name[list_end] == "K"))
    
    if (length(exclude) > 0) {
        
        list_start <- list_start[-exclude]
        list_end <- list_end[-exclude]
        
    }
    
    events <- data.frame(start = list_start, end = list_end)
    
    return(events)
    
}

sortEvents <- function(events, g)
{

    name_split <- strsplit(nodes(g)$name, ":", fixed = TRUE)
    gv_type <- sapply(name_split, "[", 1)
    gv_pos <- as.integer(sapply(name_split, "[", 3))
    st <- setdiff(sapply(name_split, "[", 4), NA)

    start_type <- gv_type[events$start]
    start_pos <- gv_pos[events$start]
    
    end_type <- gv_type[events$end]
    end_pos <- gv_pos[events$end]

    i_I <- which(start_type != "R" & end_type != "K")
    i_I <- i_I[order(end_pos[i_I] - start_pos[i_I],
        decreasing = (st == "-"))]

    i_R <- which(start_type == "R")
    i_R <- i_R[order(end_pos[i_R], decreasing = (st == "-"))]

    i_K <- which(end_type == "K")
    i_K <- i_K[order(start_pos[i_K], decreasing = (st == "+"))]
    
    i <- c(i_I, i_R, i_K)

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

##' Annotate splice variants in terms of canonical events.
##'
##' The following events are considered:
##'
##' \describe{
##'   \item{\dQuote{SE}}{skipped exon}
##'   \item{\dQuote{S2E}}{two consecutive exons skipped}
##'   \item{\dQuote{RI}}{retained intron}
##'   \item{\dQuote{MXE}}{mutually exclusive exons}
##'   \item{\dQuote{A5SS}}{alternative 5' splice site}
##'   \item{\dQuote{A3SS}}{alternative 3' splice site}
##'   \item{\dQuote{AFE}}{alternative first exon}
##'   \item{\dQuote{ALE}}{alternative last exon}
##'   \item{\dQuote{AS}}{alternative start other than \dQuote{AFE}}
##'   \item{\dQuote{AE}}{alternative end other than \dQuote{ALE}}
##' }
##'
##' For events \dQuote{SE} and \dQuote{S2E}, suffixes \dQuote{I} and
##' \dQuote{S} indicate inclusion and skipping, respectively.
##' For event \dQuote{RI} suffixes \dQuote{E} and \dQuote{R} indicate
##' exclusion and retention, respectively.
##' For events \dQuote{A5SS} and \dQuote{A3SS}, suffixes \dQuote{P} and
##' \dQuote{D} indicate use of the proximal (intron-shortening) and
##' distal (intron-lengthening) splice site, respectively. 
##' 
##' All considered events are binary events defined by two alternative
##' variants. A variant is annotated as a canonical event if it coincides
##' with one of the two variants in the canonical event, and there is at
##' least one variant in the same event that coincides with the second
##' variant of the canonical event.
##'
##' @title Annotate splice variants in terms of canonical events
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
    
    path_event_id <- eventID(variants)
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
        event_id <- intersect(path_event_id[i_1], path_event_id[i_2])
        i_1 <- i_1[path_event_id[i_1] %in% event_id]
        i_2 <- i_2[path_event_id[i_2] %in% event_id]

        path_event[i_1] <- pc(path_event[i_1], rep(event_1, length(i_1)))
        path_event[i_2] <- pc(path_event[i_2], rep(event_2, length(i_2)))

    }

    ## symmetric events
    
    for (k in seq_along(list_se_type)) {

        event <- list_se_event[k]
        type <- list_se_type[k]

        i <- unique(togroup(path_type)[grep(type, unlist(path_type))])

        if (event == "AFE") { i <- i[grep("^S", from(variants)[i])] }
        if (event == "ALE") { i <- i[grep("^E", to(variants)[i])] }
        
        if (length(i) > 0) {

            event_id_n <- table(path_event_id[i])
            event_id <- names(which(event_id_n >= 2))
            i <- i[path_event_id[i] %in% event_id]
            path_event[i] <- pc(path_event[i], rep(event, length(i)))
            
        }
        
    }

    ## alternative start/end

    i <- intersect(grep("^S", from(variants)),
        which(!any(path_event == "AFE")))
    path_event[i] <- pc(path_event[i], rep("AS", length(i)))

    i <- intersect(grep("^E", to(variants)),
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

    out <- as(tapply(as.integer(x), x_f, unique, simplify = FALSE), "list")
    
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

##' Create interpretable splice variant names
##' taking format GENE_EVENT_VARIANT/ORDER_TYPE.
##' GENE is based on geneName if available, and geneID otherwise.
##' EVENT and VARIANT enumerate events and variants for the same gene
##' and event, respectively. ORDER indicates the total number of
##' variants in the same event (e.g. 1/2 refers to the first out of two
##' splice variants in the event). TYPE is based on variantType.
##' @title Create interpretable splice variant names
##' @param variants \code{SGVariants} object
##' @return Character vector with splice variant names
##' @examples
##' makeVariantNames(sgv)
##' @keywords internal
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
