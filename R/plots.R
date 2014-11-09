exonGraph <- function(features, tx_view)
{

    if (tx_view) {

        tx_name <- txName(features)
        features <- features[togroup(tx_name)]
        mcols(features)$tx_name <- unlist(tx_name)

    } 
    
    E <- features[type(features) == "E"]
    J <- features[type(features) == "J"]
    
    v <- exonGraphNodes(E, J, tx_view)
    d <- exonGraphEdges(v, J, tx_view)    

    ## Reorder exons 

    tmp <- strsplit(v$coordinates, split = ":", fixed = TRUE)
    tmp <- strsplit(sapply(tmp, "[", 2), "-", fixed = TRUE)
    start <- as.integer(sapply(tmp, "[", 1))
    end <- as.integer(sapply(tmp, "[", 2))
    v <- v[order(start, end), ]

    g <- graph.data.frame(d = d, directed = TRUE, vertices = v)
    
    invisible(g)
    
}

exonGraphNodes <- function(E, J, tx_view)
{

    v <- data.frame(
        name = seq_along(E),
        coordinates = gr2co(E),
        type = type(E),
        stringsAsFactors = FALSE)

    v$color <- mcols(E)$color
    v$label <- mcols(E)$label

    if (tx_view) {

        v$tx_name <- mcols(E)$tx_name

    }
    
    v <- addDummyNodes(v, E, J, tx_view)
    
    return(v)

}

addDummyNodes <- function(v, E, J, tx_view)
{
    
    E_A <- gr2co(flank(E, -1, TRUE))
    E_D <- gr2co(flank(E, -1, FALSE))
    J_D <- gr2co(flank(J, -1, TRUE))
    J_A <- gr2co(flank(J, -1, FALSE))

    if (tx_view) {

        E_A <- paste0(E_A, "_", mcols(E)$tx_name)
        E_D <- paste0(E_D, "_", mcols(E)$tx_name)
        J_D <- paste0(J_D, "_", mcols(J)$tx_name)
        J_A <- paste0(J_A, "_", mcols(J)$tx_name)

    }
    
    dummy_nodes <- union(setdiff(J_D, E_D), setdiff(J_A, E_A))

    if (length(dummy_nodes) > 0) {

        if (tx_view) {

            coordinates <- sapply(strsplit(dummy_nodes, "_"), "[", 1)

        } else {

            coordinates <- dummy_nodes

        }
        
        v_dummy <- data.frame(
            name = nrow(v) + seq_along(dummy_nodes),
            coordinates = coordinates,
            type = NA,
            stringsAsFactors = FALSE)

        if (!is.null(v$color)) {

            v_dummy$color <- NA

        }
        if (!is.null(v$label)) {

            v_dummy$label <- NA

        }        
        if (tx_view) {

            v_dummy$tx_name <- sapply(strsplit(dummy_nodes, "_"), "[", 2)

        }
        
        v <- rbind(v, v_dummy)

    }

    return(v)

}

exonGraphEdges <- function(v, J, tx_view)
{
    
    E_A <- flank(co2gr(v$coordinates), -1, TRUE)
    E_D <- flank(co2gr(v$coordinates), -1, FALSE)

    if (tx_view) {

        mcols(E_A)$tx_name <- v$tx_name
        mcols(E_D)$tx_name <- v$tx_name

    }

    J_D <- flank(J, -1, TRUE)
    J_A <- flank(J, -1, FALSE)

    ol <- findOverlaps(E_D, J_D)

    if (tx_view) {

        qH <- queryHits(ol)
        sH <- subjectHits(ol)
        ol <- ol[mcols(E_D)$tx_name[qH] == mcols(J_D)$tx_name[sH]]

    }
    
    i_from <- queryHits(ol)
    i_junc <- subjectHits(ol)

    ol <- findOverlaps(J_A[i_junc], E_A)

    if (tx_view) {

        qH <- queryHits(ol)
        sH <- subjectHits(ol)
        ol <- ol[mcols(J_A)$tx_name[i_junc][qH] == mcols(E_A)$tx_name[sH]]

    }

    i_from <- i_from[queryHits(ol)]
    i_junc <- i_junc[queryHits(ol)]
    i_to <- subjectHits(ol)

    d <- data.frame(
        from = v$name[i_from],
        to = v$name[i_to],
        coordinates = gr2co(J)[i_junc],
        type = "J",
        stringsAsFactors = FALSE)

    d$color <- mcols(J)$color[i_junc]
    d$label <- mcols(J)$label[i_junc]

    if (tx_view) {

        d$tx_name <- v$tx_name[i_from]

    }
    
    return(d)
    
}

##' Plot splice graph implied by splice junctions and exon bins.
##'
##' By default, splice graph feature color is determined by annotation
##' (see arguments \code{color}, \code{color_novel}) and labels are 
##' generated automatically (see argument \code{label}). 
##'
##' Alternatively, colors and labels can be specified via elementMetadata
##' columns \dQuote{color} and \dQuote{label}, respectively.
##'
##' @title Plot splice graph
##' @param x \code{SGFeatures} or \code{TxVariants} object
##' @param geneID Single gene identifier used to subset \code{x}
##' @param geneName Single gene name used to subset \code{x}
##' @param eventID Single event identifier used to subset \code{x}
##' @param which \code{GRanges} used to subset \code{x}
##' @param toscale Controls which parts of the splice graph are drawn to
##'   scale. Possible values are \dQuote{none} (exonic and intronic regions
##'   have constant length), \dQuote{exon} (exonic regions are drawn to scale)
##'   and \dQuote{gene} (both exonic and intronic regions are drawn to scale).
##' @param label Format of exon/splice junction labels,
##'   possible values are \dQuote{id} (format E1,... J1,...), \dQuote{name}
##'   (format type:chromosome:start-end:strand), \dQuote{label} for labels
##'   specified in elementMetadata column \dQuote{label}, or \dQuote{none}
##'   for no labels.
##' @param color Color used for plotting the splice graph. Ignored if features
##'   elementMetadata column \dQuote{color} is not \code{NULL}.
##' @param color_novel Features with missing annotation are
##'   highlighted in \code{color_novel}. Ignored if features
##'   elementMetadata column \dQuote{color} is not \code{NULL}.
##' @param color_alpha Controls color transparency
##' @param color_labels Logical indicating whether label colors should
##'   be the same as feature colors
##' @param border Determines the color of exon borders, can be \dQuote{fill}
##'   (same as exon color), \dQuote{none} (no border) or a valid color name
##' @param cexLab Scale factor for feature labels
##' @param cexExon Scale factor for exon height
##' @param score \code{RLeList} containing nucleotide-level scores to be
##'   plotted with the splice graph
##' @param score_color Color used for plotting scores
##' @param score_ylim y-axis range used for plotting scores
##' @param score_ypos Numeric vector of length two, indicating the vertical
##'   position and height of the score panel, specificed as fractions of the
##'   height of the plotting region
##' @param score_nbin Number of bins for plotting scores
##' @param main Plot title
##' @param cexMain Scale factor for plot title
##' @param tx_view Plot transcripts instead of splice graph (experimental)
##' @param tx_dist Vertical distance between transcripts as fraction of height
##'   of plotting region
##' @param tx_cex Scale factor for transcript labels
##' @param asp Aspect ratio of graphics region
##' @return \code{data.frame} with plotting information for exons and
##'   splice junctions in the splice graph
##' @examples
##' \dontrun{
##'   sgf_annotated <- annotate(sgf, txf)
##'   plotSpliceGraph(sgf_annotated)
##' }
##' \dontrun{
##'   txv_annotated <- annotate(txv, txf)
##'   plotSpliceGraph(txv_annotated)
##' }
##' @author Leonard Goldstein

plotSpliceGraph <- function(x, geneID = NULL, geneName = NULL,
    eventID = NULL, which = NULL, toscale = c("exon", "none", "gene"),
    label = c("id", "name", "label", "none"), color = "grey",
    color_novel = "red", color_alpha = 0.8, color_labels = FALSE,
    border = "fill", cexLab = 1, cexExon = 1, score = NULL,
    score_color = "darkblue", score_ylim = NULL, score_ypos = c(0.2, 0.1),
    score_nbin = 400, main = NULL, cexMain = 1, tx_view = FALSE,
    tx_dist = 0.2, tx_cex = 1, asp = 1)
{

    toscale <- match.arg(toscale)
    label <- match.arg(label)

    if (!is(x, "SGFeatures") && !is(x, "TxVariants")) {

        stop("x must be an SGFeatures or TxVariants object")

    }
    
    x <- subsetFeatures(x, geneID, eventID, which, geneName)

    if (length(x) == 0) { return() }

    if (is(x, "TxVariants")) {
        
        x <- extractFeaturesFromVariants(x)

    }

    x <- setFeatureColors(x, color, color_novel, color_alpha)
        
    g <- exonGraph(x, tx_view)

    pars <- getPlottingParameters(g, toscale, border, asp, tx_view, tx_dist)

    plot(g,
        xlim = c(-1, 1) / 1.08,
        ylim = c(-1, 1) / 1.08,
        rescale = FALSE,
        layout = cbind(pars$vertex_x, pars$vertex_y),
        vertex.size = pars$vertex_size * 100,
        vertex.size2 = cexExon * 16,
        vertex.color = pars$vertex_color,
        vertex.frame.color = pars$vertex_frame_color,
        vertex.shape = pars$vertex_shape,
        vertex.label = NA,
        edge.color = pars$edge_color,
        edge.width = 2,
        edge.arrow.size = 0.5,
        edge.curved = pars$edge_curved,
        edge.label = NA,
        edge.arrow.mode = pars$edge_arrow_mode,
        asp = asp)

    if (!is.null(score)) {

        plotScore(pars, score, score_color, score_ylim, score_ypos, score_nbin)

    }
    
    text(x = 0, y = 0.95, labels = main, pos = 1, offset = 0, font = 2,
        cex = cexMain)

    df <- getGraphInfo(g, pars, color_labels, tx_view)
        
    if (label != "none") {

        text(x = df$x, y = df$y, labels = df[, label],
            pos = c(E = 1, J = 3)[substr(df$id, 1, 1)],
            offset = 1, cex = cexLab, col = df$color)
        
    }

    if (tx_view) {

        df_collapsed <- unique(df[, c("y", "tx_name")])
        mtext(side = 4, at = df_collapsed$y, text = df_collapsed$tx_name,
            cex = 0.8 * tx_cex, line = 0.5, las = 1)

    }
    
    invisible(df)
    
}

getPlottingParameters <- function(g, toscale, border, asp, tx_view, tx_dist)
{
    
    ## data frames of nodes and edges
    gv <- nodes(g)
    gd <- edges(g)

    i_from <- match(gd$from, gv$name)
    i_to <- match(gd$to, gv$name)
    
    ## exon and gene coordinates
    gr_exon <- co2gr(gv$coordinates)
    gr_gene <- range(gr_exon)

    if (length(gr_gene) > 1) {

        stop("features must be on the same chromosome and strand")

    }
    
    ## exon coordinates relative to gene locus
    ir_exon <- ranges(mapCoords(gr_exon, split(gr_gene, 1)))

    ## exonic regions
    ir_exonic <- reduce(ir_exon)

    ## exon x-coordinates and widths
    exon_coordinates <- xcoordinates(ir_exon, ir_exonic, toscale)
    exon_width <- exon_coordinates$w
    exon_x <- exon_coordinates$x
    
    if (tx_view) {

        exon_y <- as.integer(as.factor(as.character(gv$tx_name)))
        exon_y <- exon_y - 1
        exon_y <- exon_y * tx_dist * 2
        exon_y <- exon_y - mean(unique(exon_y))
        edge_curved <- rep(0, nrow(gd))

    } else {

        exon_y <- rep(0, nrow(gv))
        edge_curved <- as.numeric(rank(exon_x)[i_to] >
            rank(exon_x)[i_from] + 1)
        edge_curved <- edge_curved / asp

    }

    ## vertex plotting parameters
    exon_shape <- ifelse(!is.na(gv$type), "rectangle", "none")
    exon_color <- gv$color
    exon_border_color <- switch(border, fill = gv$color, none = NA, border)
    
    ## edge plotting parameters
    edge_color <- gd$color
    edge_arrow_mode <- rep(0, nrow(gd))

    pars <- list(
        vertex_x = exon_x,
        vertex_y = exon_y,
        vertex_size = exon_width,
        vertex_shape = exon_shape,
        vertex_color = exon_color,
        vertex_frame_color = exon_border_color,
        edge_color = edge_color,
        edge_arrow_mode = edge_arrow_mode,
        edge_curved = edge_curved,
        gene_chrom = as.character(seqnames(gr_gene)),
        gene_strand = as.character(strand(gr_gene)),
        exon_start = start(gr_exon),
        exon_end = end(gr_exon)
    )

    if (tx_view) {

        pars$vertex_tx_name <- gv$tx_name

    }
    
    return(pars)

}

plotScore <- function(pars, score, color, ylim, ypos, nbin)
{

    n_exon <- length(pars$vertex_x)    

    tmp <- lapply(seq_len(n_exon), scorePerExon, pars, score)

    bp_x_start <- do.call(c, lapply(tmp, "[[", "x_start"))
    bp_x_end <- do.call(c, lapply(tmp, "[[", "x_end"))
    bp_y <- do.call(c, lapply(tmp, "[[", "y"))
    
    bin_size <- 2 / nbin
    bin_breaks <- seq(-1, 1, length.out = nbin + 1)

    bp_bin <- IRanges(
        as.integer(cut(bp_x_start, bin_breaks, right = FALSE)),
        as.integer(cut(bp_x_end, bin_breaks, right = TRUE)))

    bin_y <- tapply(bp_y[rep(seq_along(bp_bin), width(bp_bin))],
        as.integer(bp_bin), mean)
    
    if (is.null(ylim)) { ylim <- range(bin_y, na.rm = TRUE, finite = TRUE) }

    bin_y[bin_y < ylim[1]] <- ylim[1]
    bin_y[bin_y > ylim[2]] <- ylim[2]

    ypos_2 <- c(-1 + ypos[1] * 2, ypos[2] * 2)
    
    bin_y <- bin_y - ylim[1]
    bin_y <- bin_y / diff(range(ylim)) * ypos_2[2]
    bin_y <- bin_y + ypos_2[1]
        
    mapply(rect, xleft = bin_breaks[-length(bin_breaks)],
        xright = bin_breaks[-1], ytop = bin_y,
        MoreArgs = list(ybottom = ypos_2[1], col = color, border = NA))

    labels <- names(ylim)
    
    if (is.null(labels)) { labels <- ylim }
        
    axis(side = 2, at = c(ypos_2[1], ypos_2[1] + ypos_2[2]), labels = labels,
        mgp = c(3, 0.5, 0.5), tcl = -0.25, las = 1)
    
}

scorePerExon <- function(i, pars, score)
{

    n_exon <- length(pars$vertex_x)    

    ir_exon <- IRanges(pars$exon_start[i], pars$exon_end[i])
    y <- as.numeric(score[[pars$gene_chrom]][ir_exon])
    if (pars$gene_strand == "-") { y <- rev(y) }
    
    ex_start <- pars$vertex_x[i] - 0.5 * pars$vertex_size[i]
    ex_end <- pars$vertex_x[i] + 0.5 * pars$vertex_size[i]
    
    d <- (ex_end - ex_start) / length(y)
    
    ex_bp_x_start <- seq(ex_start, ex_end - d, length.out = length(y))
    ex_bp_x_end <- seq(ex_start + d, ex_end, length.out = length(y))
    ex_bp_y <- y
    
    if (i < n_exon && pars$exon_start[i + 1] - pars$exon_end[i] > 1) {
        
        ir_intron <- IRanges(pars$exon_end[i], pars$exon_start[i + 1]) - 1
        y <- as.numeric(score[[pars$gene_chrom]][ir_intron])
        if (pars$gene_strand == "-") { y <- rev(y) }
        
        if (pars$gene_strand == "+") {
            
            in_start <- pars$vertex_x[i] + 0.5 * pars$vertex_size[i]
            in_end <- pars$vertex_x[i + 1] - 0.5 * pars$vertex_size[i + 1]
            
        } else {
            
            in_start <- pars$vertex_x[i + 1] + 0.5 * pars$vertex_size[i + 1]
            in_end <- pars$vertex_x[i] - 0.5 * pars$vertex_size[i]
            
        }
        
        d <- (in_end - in_start) / length(y)
        
        in_bp_x_start <- seq(in_start, in_end - d, length.out = length(y))
        in_bp_x_end <- seq(in_start + d, in_end, length.out = length(y))
        in_bp_y <- y
        
        out_bp_x_start <- c(ex_bp_x_start, in_bp_x_start)
        out_bp_x_end <- c(ex_bp_x_end, in_bp_x_end)
        out_bp_y <- c(ex_bp_y, in_bp_y)
        
    } else {
        
        out_bp_x_start <- ex_bp_x_start
        out_bp_x_end <- ex_bp_x_end
        out_bp_y <- ex_bp_y
        
    }
    
    out_bp_x_start[which(out_bp_x_start < -1)] <- -1
    out_bp_x_start[which(out_bp_x_start > 1)] <- 1
    out_bp_x_end[which(out_bp_x_end < -1)] <- -1
    out_bp_x_end[which(out_bp_x_end > 1)] <- 1

    out <- list(x_start = out_bp_x_start, x_end = out_bp_x_end, y = out_bp_y)

    return(out)
    
}

getGraphInfo <- function(g, pars, color_labels, tx_view)
{

    ## data frames of nodes and edges
    gv <- nodes(g)
    gd <- edges(g)

    i_from <- match(gd$from, gv$name)
    i_to <- match(gd$to, gv$name)

    x_exon <- pars$vertex_x
    y_exon <- pars$vertex_y
    w_exon <- pars$vertex_size
    
    ## coordinates for junction labels
    x_junc <- 0.5 * (x_exon[i_from] + 0.5 * w_exon[i_from] +
        x_exon[i_to] - 0.5 * w_exon[i_to])
    w_junc <- x_exon[i_to] - 0.5 * w_exon[i_to] -
        (x_exon[i_from] + 0.5 * w_exon[i_from])
    y_junc <- 0.5 * (y_exon[i_from] + y_exon[i_to]) +
        pars$edge_curved * w_junc * 1/3

    ## output
    i_exon <- which(!is.na(gv$type))
    i_exon <- i_exon[order(x_exon)]
    i_junc <- order(x_exon[i_from], w_junc)

    co_E <- gv$coordinates[i_exon]
    co_J <- gd$coordinates[i_junc]
    
    df <- data.frame(
        id = c(paste0("E", as.integer(factor(co_E, levels = unique(co_E)))),
            paste0("J", as.integer(factor(co_J, levels = unique(co_J))))),
        name = paste0(c(gv$type[i_exon], gd$type[i_junc]), ":",
            c(gv$coordinates[i_exon], gd$coordinates[i_junc])),
        x = c(x_exon[i_exon], x_junc[i_junc]),
        y = c(y_exon[i_exon], y_junc[i_junc]))
    
    if (color_labels) {

        df$color <- addAlpha(c(gv$color[i_exon], gd$color[i_junc]), 1)

    } else {

        df$color <- "black"

    }
    
    if (!is.null(gv$label) && !is.null(gd$label)) {
        
        df$label <- c(gv$label[i_exon], gd$label[i_junc])
        
    }

    if (tx_view) {

        df$tx_name <- c(gv$tx_name[i_exon], gd$tx_name[i_junc])

    }

    return(df)

}

xcoordinates <- function(ir, ir_exonic, toscale)
{

    n_exonic <- length(ir_exonic)

    ## intronic regions
    ir_intronic <- IRanges(end(ir_exonic)[-n_exonic], start(ir_exonic)[-1]) - 1 
    n_intronic <- length(ir_intronic)
    
    ## set widths for exonic/intronic regions
    ## if toscale is 'none' or 'exon', half of the x-axis range is used
    ## for plotting exonic regions and half for intronic regions
    if (toscale == "none") {
        
        w_exonic <- rep(1 / n_exonic, n_exonic)
        w_intronic <- rep(1 / n_intronic, n_intronic)
        
    } else if (toscale == "exon") {
        
        w_exonic <- width(ir_exonic) / sum(width(ir_exonic))
        w_intronic <- rep(1 / n_intronic, n_intronic)
        
    } else if (toscale == "gene") {
        
        w_exonic <- 2 * width(ir_exonic) / width(range(ir_exonic))
        w_intronic <- 2 * width(ir_intronic) / width(range(ir_exonic))
        
    }

    ## if there are no intronic regions, use full x-axis range
    ## for exonic regions
    if (toscale %in% c("none", "exon") && n_intronic == 0) {

        w_exonic <- 2

    }
    
    ## x-coordinates for start of exonic regions
    x_exonic <- c(0, cumsum(w_exonic[-n_exonic] + w_intronic)) - 1
    
    ## map to exonic regions 
    ol <- unlist(as.list(findOverlaps(ir, ir_exonic)))
    
    ## width
    rel_width <- width(ir) / width(ir_exonic)[ol]
    w <- rel_width * w_exonic[ol]
    
    ## x-coordinates (center)
    rel_offset <- (start(ir) - start(ir_exonic)[ol]) / width(ir_exonic)[ol]
    x <- x_exonic[ol] + rel_offset * w_exonic[ol] + 0.5 * w

    out <- list("x" = x, "w" = w)
    
    return(out)
    
}

setFeatureColors <- function(features, color, color_novel, alpha)
{

    if (is.null(mcols(features)$color)) {

        features_color <- rep(color, length(features))

        if (!is.null(color_novel) && !is.null(txName(features))) {

            txName <- txName(features)
            i_novel <- which(elementLengths(txName) == 0)
            features_color[i_novel] <- color_novel

        }

    } else {

        features_color <- mcols(features)$color

    }

    if (!is.null(alpha)) {

        features_color <- addAlpha(features_color, alpha)
            
    }
    
    mcols(features)$color <- features_color

    return(features)

}

addAlpha <- function(col, alpha)
{

    col_rgb <- col2rgb(col)/255
    rgb(col_rgb[1, ], col_rgb[2, ], col_rgb[3, ], alpha)
    
}

##' @title Plot splice graph and heatmap of expression values
##' @inheritParams plotSpliceGraph
##' @param x \code{SGFeatureCounts} object 
##' @param assay Name of assay to be plotted in the heatmap
##' @param include Include \dQuote{exons}, \dQuote{junctions} or
##'   \dQuote{both} in the heatmap
##' @param transform Transformation applied to assay data
##' @param Rowv Determines order of rows. Either a vector of values used to
##'   reorder rows, or \code{NA} to suppress reordering, or \code{NULL} for
##'   hierarchical clustering.
##' @param distfun Distance function used for hierarchical clustering
##'   of rows (samples)
##' @param hclustfun Clustering function used for hierarchical clustering
##'   of rows (samples)
##' @param margin Width of right-hand margin as fraction of width of the
##'   graphics device. Ignored if \code{square} is \code{TRUE}.
##' @param RowSideColors Character vector (or list of character vectors)
##'   with length(s) equal to \code{ncol(x)} containing color names for
##'   horizontal side bars for sample annotation
##' @param square Logical, if \code{TRUE} margins are set such that
##'   cells in the heatmap are square
##' @param cexRow Scale factor for row (sample) labels
##' @param cexCol Scale factor for column (feature) labels
##' @param labRow Character vector of row (sample) labels
##' @param col Heatmap colors
##' @param zlim Range of values for which colors should be plotted,
##'   if \code{NULL} range of finite values
##' @param heightTopPanel Height of top panel as fraction of height of the
##'   graphics device
##' @return Return value of \code{plotSpliceGraph}
##' @examples
##' \dontrun{
##'   sgfc_annotated <- annotate(sgfc, txf)
##'   plotFeatures(sgfc_annotated)
##' }
##' @author Leonard Goldstein

plotFeatures <- function(x, geneID = NULL, geneName = NULL,
    which = NULL, toscale = c("exon", "none", "gene"), color = "grey",
    color_novel = "red", color_alpha = 0.8, color_labels = FALSE,
    border = "fill", cexLab = 1, cexExon = 1, score = NULL,
    score_color = "darkblue", score_ylim = NULL, score_ypos = c(0.2, 0.1),
    score_nbin = 400, main = NULL, cexMain = 1, tx_view = FALSE,
    tx_dist = 0.1, tx_cex = 1, assay = "FPKM",
    include = c("junctions", "exons", "both"),
    transform = function(x) log2(x + 1), Rowv = NULL,
    distfun = dist, hclustfun = hclust, margin = 0.2,
    RowSideColors = NULL, square = FALSE, cexRow = 1, cexCol = 1,
    labRow = colnames(x), col = colorRampPalette(c("black", "gold"))(256),
    zlim = NULL, heightTopPanel = 0.3)
{

    toscale <- match.arg(toscale)
    include <- match.arg(include)
    
    if (!is(x, "SGFeatureCounts")) {

        stop("x must be an SGFeatureCounts object")

    }

    x <- subsetFeatures(x, geneID = geneID, which = which,
        geneName = geneName)

    if (nrow(x) == 0) { return() }

    features <- rowData(x)
    
    if (!is.null(RowSideColors) && !is.list(RowSideColors)) {

        RowSideColors <- list(RowSideColors)

    }    

    n_samples <- ncol(x)
    include_type <- switch(include, junctions = "J", exons = "E",
        both = c("E", "J"))
    n_features <- length(which(type(features) %in% include_type))
    
    pars <- getLayoutParameters(n_samples, n_features, margin, heightTopPanel,
        RowSideColors, square, tx_view)

    layout(pars$mat, heights = pars$hei, widths = pars$wid)

    par(mai = pars$mai)
    
    df <- plotSpliceGraph(x = features, toscale = toscale,
        label = "id", color = color, color_novel = color_novel,
        color_alpha = color_alpha, color_labels = color_labels,
        border = border, cexLab = cexLab, cexExon = cexExon,
        score = score, score_color = score_color,
        score_ylim = score_ylim, score_ypos = score_ypos,
        score_nbin = score_nbin, main = main, cexMain = cexMain,
        tx_view = tx_view, tx_dist = tx_dist, tx_cex = tx_cex,
        asp = pars$asp)

    i_df <- switch(include,
        exons = grep("^E", df$name),
        junctions = grep("^J", df$name),
        both = seq_len(nrow(df)))
    i_df <- i_df[!duplicated(df$name[i_df])]
    i <- match(df$name[i_df], feature2name(features))
    
    X <- transform(assay(x[i, ], assay))
    labCol <- df$id[i_df]
    colLabCol <- df$color[i_df]
        
    plotImage(X, Rowv = Rowv, distfun = distfun,
        hclustfun = hclustfun, RowSideColors = RowSideColors,
        cexRow = cexRow, cexCol = cexCol, labRow = labRow,
        labCol = labCol, colLabCol = colLabCol,
        col = col, zlim = zlim)

    invisible(df)

}

##' @title Plot splice graph and heatmap of variant frequencies
##' @inheritParams plotSpliceGraph
##' @inheritParams plotFeatures
##' @param x \code{TxVariantCounts} object 
##' @param transform Transformation applied to variant frequencies
##' @param expand_variants Experimental option - leave set to \code{FALSE}
##' @return Return value of \code{plotSpliceGraph}
##' @examples
##' \dontrun{
##'   txvc_annotated <- annotate(txvc, txf)
##'   plotVariants(txvc_annotated)
##' }
##' @author Leonard Goldstein

plotVariants <- function(x, eventID = NULL,
    toscale = c("exon", "none", "gene"), color = "grey", color_novel = "red", 
    color_alpha = 0.8, color_labels = FALSE, border = "fill", 
    cexLab = 1, cexExon = 1, score = NULL, score_color = "darkblue",
    score_ylim = NULL, score_ypos = c(0.2, 0.1), score_nbin = 400,
    main = NULL, cexMain = 1, tx_view = FALSE, tx_dist = 0.1, tx_cex = 1,
    transform = function(x) asin(sqrt(x)), Rowv = NULL,
    distfun = dist, hclustfun = hclust, margin = 0.2,
    RowSideColors = NULL, square = FALSE, cexRow = 1, cexCol = 1,
    labRow = colnames(x), col = colorRampPalette(c("black", "gold"))(256),
    zlim = NULL, heightTopPanel = 0.3, expand_variants = FALSE)
{

    toscale <- match.arg(toscale)

    if (!is(x, "TxVariantCounts")) {

        stop("x must be a TxVariantCounts object")

    }
    
    x <- subsetFeatures(x, eventID = eventID,
        expand_variants = expand_variants)

    if (nrow(x) == 0) { return() }    

    if (expand_variants) {
    
        x <- expandSE(x, eventID)

    }

    X <- transform(variantFreq(x))

    not_quantifiable <- elementLengths(featureID5p(x)) == 0 &
        elementLengths(featureID3p(x)) == 0
    
    if (all(not_quantifiable)) {

        stop("none of the selected transcript variants could be quantified")

    }

    if (any(not_quantifiable)) {

        warning("one or more transcript variants could not be quantified")

    }

    if (any(elementLengths(featureID5p(x)) == 0 &
        elementLengths(featureID3p(x)) == 0)) {

        warning("one or more transcript variants could not be quantified")

    }
        
    if (!is.null(RowSideColors) && !is.list(RowSideColors)) {

        RowSideColors <- list(RowSideColors)

    }    

    n_samples <- ncol(x)
    n_features <- nrow(x)
    
    pars <- getLayoutParameters(n_samples, n_features, margin, heightTopPanel,
        RowSideColors, square, tx_view)

    layout(pars$mat, heights = pars$hei, widths = pars$wid)

    par(mai = pars$mai)
    
    df <- plotSpliceGraph(x = rowData(x), toscale = toscale,
        label = "label", color = color, color_novel = color_novel,
        color_alpha = color_alpha, color_labels = color_labels,
        border = border, cexLab = cexLab, cexExon = cexExon,
        score = score, score_color = score_color,
        score_ylim = score_ylim, score_nbin = score_nbin,
        score_ypos = score_ypos, main = main, cexMain = cexMain,
        tx_view = tx_view, tx_dist = tx_dist, tx_cex = tx_cex,
        asp = pars$asp)
    
    labCol <- seq_len(nrow(X))
    colLabCol <- "black"
    
    plotImage(X, Rowv = Rowv, distfun = distfun,
        hclustfun = hclustfun, RowSideColors = RowSideColors,
        cexRow = cexRow, cexCol = cexCol, labRow = labRow,
        labCol = labCol, colLabCol = colLabCol,
        col = col, zlim = zlim)

    invisible(df)

}

extractFeaturesFromVariants <- function(variants)
{
    
    id2variant <- tapply(togroup(variants),
        featureID(unlist(variants)), unique, simplify = FALSE)
    id2variant <- id2variant[elementLengths(id2variant) == 1]
    features <- uniqueFeatures(unlist(variants))
    mcols(features)$label <- id2variant[match(
        featureID(features), names(id2variant))]

    return(features)
    
}

getLayoutParameters <- function(n_samples, n_features, margin, heightTopPanel,
    RowSideColors, square, tx_view)
{

    n_RSC <- length(RowSideColors)
    
    mat <- rbind(rep(1, 3), rep(0, 3), c(0, 2, 0), rep(0, 3))
    
    for (i in seq_len(n_RSC)) {
        
        n <- max(mat) + 1            
        mat <- cbind(c(1, 0, 0, 0), c(1, 0, n, 0), mat)
        
    }

    mar_wid <- c(0.05, margin)
    mar_hei <- c(0.05, 0.05) * (1 - heightTopPanel)
    
    img_wid <- 1 - sum(mar_wid) - 0.025 * n_RSC * 2
    img_hei <- 1 - heightTopPanel - sum(mar_hei)
    
    wid <- c(mar_wid[1], rep(0.025, n_RSC * 2), img_wid, mar_wid[2])
    hei <- c(heightTopPanel, mar_hei[1], img_hei, mar_hei[2])

    if (square) {

        padding <- getPadding(n_samples, n_features, wid, hei)
        
        i_wid_mar <- length(wid) - c(2, 0)
        i_wid_img <- length(wid) - 1

        wid[i_wid_mar] <- wid[i_wid_mar] + 0.5 * padding[1]
        wid[i_wid_img] <- wid[i_wid_img] - padding[1]

        i_hei_mar <- c(2, 4)
        i_hei_img <- 3

        hei[i_hei_mar] <- hei[i_hei_mar] + 0.5 * padding[2]
        hei[i_hei_img] <- hei[i_hei_img] - padding[2]
        
    }

    ds <- dev.size()
    dev_wid <- ds[1]
    dev_hei <- ds[2]

    mai_hei <- c(0.05, 0.05) * dev_hei * heightTopPanel
    
    if (tx_view) {

        mai_wid <- c(0.05, margin) * dev_wid
        
    } else {
        
        mai_wid <- c(0.05, 0.05) * dev_wid
        
    }

    mai <- c(mai_hei[1], mai_wid[1], mai_hei[2], mai_wid[2])
    asp <- (dev_hei * heightTopPanel - sum(mai_hei)) / (dev_wid - sum(mai_wid))
    
    list(mat = mat, wid = wid, hei = hei, mai = mai, asp = asp)
    
}

getPadding <- function(n_samples, n_features, wid, hei)
{

    ds <- dev.size()
    dev_wid <- ds[1]
    dev_hei <- ds[2]

    img_wid <- wid[length(wid) - 1]
    img_hei <- hei[3]

    img_wid_real <- img_wid / sum(wid) * dev_wid
    img_hei_real <- img_hei / sum(hei) * dev_hei
    
    cel_wid_real <- img_wid_real / n_features
    cel_hei_real <- img_hei_real / n_samples
    
    if (cel_wid_real > cel_hei_real) {

        pad_wid <- (1 - cel_hei_real / cel_wid_real) * img_wid
        pad_hei <- 0

    } else {

        pad_wid <- 0
        pad_hei <- (1 - cel_wid_real / cel_hei_real) * img_hei
        
    }

    pad <- c(pad_wid, pad_hei)

    return(pad)

}

plotImage <- function(x, Rowv = NA, distfun, hclustfun, RowSideColors,
    cexRow, cexCol, labRow = NULL, labCol = NULL, colLabCol, col, zlim)
{
    
    if (is.null(Rowv)) {

        if (ncol(x) > 1) {
            
            j <- hclustfun(distfun(t(x)))$order

        } else {

            j <- 1

        }

    } else if (is.na(Rowv)) {

        j <- rev(seq_len(ncol(x)))

    } else {

        j <- Rowv

    }

    x <- x[, j]
    RowSideColors <- lapply(RowSideColors, "[", j)
    
    if (!is.matrix(x)) { x <- matrix(x, ncol = length(j)) }        
    if (is.null(zlim)) { zlim <- range(x[is.finite(x)]) }

    par(mai = c(0, 0, 0, 0))
    image(x, col = col, zlim = zlim, axes = FALSE)

    if (!is.null(labCol)) {
    
        mtext(side = 3,
            at = seq(from = 0, to = 1, length.out = length(labCol)),
            text = labCol, cex = 0.8 * cexCol, col = colLabCol,
            line = 0.25, las = 3)

    }

    if (!is.null(labRow)) {

        mtext(side = 4, at = seq(from = 0, to = 1, length.out = ncol(x)),
            text = labRow[j], cex = 0.8 * cexRow, line = 0.5, las = 1)

    }

    for (r in seq_along(RowSideColors)) {

        par(mai = c(0, 0, 0, 0))
        image(matrix(seq_len(ncol(x)), nrow = 1),
            col = RowSideColors[[r]], axes = FALSE)
        mtext(side = 3, text = names(RowSideColors)[r],
            cex = 0.8 * cexCol, line = 0.25, las = 3)
        
    }

}

subsetFeatures <- function(x, geneID = NULL, eventID = NULL, which = NULL,
    geneName = NULL, expand_variants = FALSE)
{

    if (is(x, "Counts")) {

        y <- rowData(x)

    } else {

        y <- x

    }

    if (!is.null(geneID) || !is.null(eventID) || !is.null(geneName)) {
    
        if (!is.null(geneID)) {
            
            i <- which(geneID(y) %in% geneID)

        } else if (!is.null(eventID)) {

            i <- which(eventID(y) %in% eventID)

            if (length(eventID) > 1) {

                featureIDs <- featureID(y)
                featureIDs <- gsub("(", "", featureIDs, fixed = TRUE)
                featureIDs <- gsub(")", "", featureIDs, fixed = TRUE)
                featureIDs <- gsub("|", ",", featureIDs, fixed = TRUE)
                featureIDs <- unlist(strsplit(featureIDs, ",", fixed = TRUE))

                if (any(duplicated(featureIDs))) {
                    
                    stop("cannot plot overlapping events")

                }

            }
            
            if (expand_variants) {
            
                featureIDs <- unique(unlist(strsplit(
                    unlist(expandPath(featureID(y)[i])), ",", fixed = TRUE)))
                i <- unlist(lapply(featureIDs,
                    function (id) { grep(paste0("(^|,)", id, "(,|$)"),
                        featureID(y)) }))
                eventIDs <- unique(eventID(y)[i])
                i <- which(eventID(y) %in% eventIDs)

            }
            
        } else if (!is.null(geneName)) {

            if (is.null(geneName(y))) {
                
                stop("missing geneName")

            }
        
            if (length(geneName) > 1) {
                
                stop("geneName must be of length 1")

            }
            
            i <- which(any(geneName(y) == geneName))

        }

    } else if (!is.null(which)) {

        i <- which(y %over% which)
                        
    } else {

        return(x)

    }

    if (is(x, "Counts")) {

        x <- x[i, ]
        
    } else {
        
        x <- x[i]
        
    }
    
    return(x)

}
