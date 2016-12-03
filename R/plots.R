exonGraph <- function(features, tx_view)
{

    if (tx_view) {

        tx_name <- txName(features)
        i <- which(elementNROWS(tx_name) == 0)
        tx_name[i] <- as(feature2name(features[i]), "CharacterList")
        features <- features[togroup0(tx_name)]
        mcols(features)$tx_name <- unlist(tx_name)

    }

    E <- features[type(features) == "E"]
    J <- features[type(features) == "J"]

    v <- exonGraphNodes(E, J, tx_view)
    d <- exonGraphEdges(v, J, tx_view)

    g <- graph.data.frame(d = d, directed = TRUE, vertices = v)

    invisible(g)

}

exonGraphNodes <- function(E, J, tx_view)
{

    v <- data.frame(
        name = seq_along(E),
        coordinates = gr2co(E),
        type = type(E),
        featureID = featureID(E),
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
            featureID = NA_integer_,
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

        v <- rbindDfsWithoutRowNames(v, v_dummy)

    }

    v <- v[order(co2gr(v$coordinates)), ]

    return(v)

}

exonGraphEdges <- function(v, J, tx_view)
{

    if (length(J) == 0) {

        d <- data.frame(
            from = character(),
            to = character(),
            stringsAsFactors = FALSE)

        return(d)

    }

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
        featureID = featureID(J)[i_junc],
        stringsAsFactors = FALSE)

    d$color <- mcols(J)$color[i_junc]
    d$label <- mcols(J)$label[i_junc]

    if (tx_view) {

        d$tx_name <- v$tx_name[i_from]

    }

    return(d)

}

##' Plot the splice graph implied by splice junctions and exon bins.
##' Invisibly returns a \code{data.frame} with details of plotted
##' features, including genomic coordinates.
##'
##' By default, the color of features in the splice graph is
##' determined by annotation status (see arguments \code{color},
##' \code{color_novel}) and feature labels are generated automatically
##' (see argument \code{label}). Alternatively, colors and labels can
##' be specified via metadata columns \dQuote{color} and
##' \dQuote{label}, respectively.
##'
##' @title Plot splice graph
##' @param x \code{SGFeatures} or \code{SGVariants} object
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
##'   specified in metadata column \dQuote{label}, or \dQuote{none}
##'   for no labels.
##' @param color Color used for plotting the splice graph. Ignored if features
##'   metadata column \dQuote{color} is not \code{NULL}.
##' @param color_novel Features with missing annotation are
##'   highlighted in \code{color_novel}. Ignored if features
##'   metadata column \dQuote{color} is not \code{NULL}.
##' @param color_alpha Controls color transparency
##' @param color_labels Logical indicating whether label colors should
##'   be the same as feature colors
##' @param border Determines the color of exon borders, can be \dQuote{fill}
##'   (same as exon color), \dQuote{none} (no border), or a valid color name
##' @param curvature Numeric determining curvature of plotted splice junctions.
##' @param ypos Numeric vector of length two, indicating the vertical
##'   position and height of the exon bins in the splice graph,
##'   specificed as fraction of the height of the plotting region
##'   (not supported for \code{tx_view = TRUE})
##' @param score \code{RLeList} containing nucleotide-level scores
##'   to be plotted with the splice graph
##' @param score_color Color used for plotting scores
##' @param score_ylim Numeric vector of length two, determining y-axis range
##'   for plotting scores
##' @param score_ypos Numeric vector of length two, indicating the vertical
##'   position and height of the score panel, specificed as fraction of the
##'   height of the plotting region
##' @param score_nbin Number of bins for plotting scores
##' @param score_summary Function used to calculate per-bin score summaries
##' @param score_label Label used to annotate score panel
##' @param ranges \code{GRangesList} to be plotted with the splice graph
##' @param ranges_color Color used for plotting ranges
##' @param ranges_ypos Numeric vector of length two, indicating the vertical
##'   position and height of the ranges panel, specificed as fraction of the
##'   height of the plotting region
##' @param main Plot title
##' @param tx_view Plot transcripts instead of splice graph (experimental)
##' @param tx_dist Vertical distance between transcripts as fraction of height
##'   of plotting region
##' @param short_output Logical indicating whether the returned data frame
##'   should only include information that is likely useful to the user
##' @return \code{data.frame} with information on exon bins and
##'   splice junctions included in the splice graph
##' @examples
##' \dontrun{
##' sgf_annotated <- annotate(sgf_pred, txf_ann)
##' plotSpliceGraph(sgf_annotated)
##' }
##' \dontrun{
##' sgv_annotated <- annotate(sgv_pred, txf_ann)
##' plotSpliceGraph(sgv_annotated)
##' }
##' NULL
##' @author Leonard Goldstein

plotSpliceGraph <- function(x, geneID = NULL, geneName = NULL,
    eventID = NULL, which = NULL, toscale = c("exon", "none", "gene"),
    label = c("id", "name", "label", "none"), color = "gray",
    color_novel = color, color_alpha = 0.8, color_labels = FALSE,
    border = "fill", curvature = NULL, ypos = c(0.5, 0.1),
    score = NULL, score_color = "darkblue",
    score_ylim = NULL, score_ypos = c(0.3, 0.1), score_nbin = 200,
    score_summary = mean, score_label = NULL, ranges = NULL,
    ranges_color = "darkblue", ranges_ypos = c(0.1, 0.1),
    main = NULL, tx_view = FALSE, tx_dist = 0.2, short_output = TRUE)
{

    toscale <- match.arg(toscale)
    label <- match.arg(label)

    if (!is(x, "SGFeatures") && !is(x, "SGVariants")) {

        stop("x must be an SGFeatures or SGVariants object")

    }

    if (is(x, "SGVariants")) x <- updateObject(x, verbose = TRUE)

    if (!is.null(score) && !is(score, "RleList")) {

        stop("score must be an RleList or NULL")

    }

    if (!is.null(ranges) && !is(ranges, "GRangesList")) {

        stop("ranges must be a GRangesList or NULL")

    }

    x <- restrictFeatures(x, geneID, eventID, which, geneName)

    if (length(x) == 0) { return() }

    if (is(x, "SGVariants")) {

        x <- extractFeaturesFromVariants(x)

    }

    x <- setFeatureColors(x, color, color_novel, color_alpha)

    g <- exonGraph(x, tx_view)

    exon_coordinates <- getExonCoordinates(g, toscale)

    plot(NA, xlim = c(-1, 1), ylim = c(-1, 1), xaxs = "i", yaxs = "i",
        axes = FALSE, xlab = NA, ylab = NA)

    df <- plotExonGraph(g, exon_coordinates, ypos, "both", border,
        curvature, label, color_labels, 1, 3, tx_view, tx_dist)

    if (!is.null(score)) {

        plotTrackScore(exon_coordinates, score, score_color, score_ylim,
            score_ypos, score_nbin, score_summary, score_label)

    }

    if (!is.null(ranges)) {

        plotTrackRanges(exon_coordinates, toscale, ranges, ranges_color,
            ranges_ypos)

    }

    text(x = 0, y = 0.95, labels = main, pos = 1, offset = 0, font = 2)

    if (short_output) { df <- df[c("id", "name", "type", "featureID")] }

    invisible(df)

}

plotExonGraph <- function(g, exon_coordinates, ypos, include, border,
    curvature, label, color_labels, label_pos_exon, label_pos_junction,
    tx_view, tx_dist)
{

    ## data frames of nodes and edges
    gv <- nodes(g)
    gd <- edges(g)

    i_from <- match(gd$from, gv$name)
    i_to <- match(gd$to, gv$name)

    if (tx_view) {

        vertex_y <- as.integer(as.factor(as.character(gv$tx_name)))
        vertex_y <- vertex_y - 1
        vertex_y <- vertex_y * 2 * tx_dist
        vertex_y <- vertex_y - mean(unique(vertex_y))
        edge_short <- rep(TRUE, nrow(gd))

    } else {

        vertex_y <- rep(-1 + 2 * ypos[1], nrow(gv))
        edge_short <- rank(exon_coordinates$vertex_x)[i_to] ==
            rank(exon_coordinates$vertex_x)[i_from] + 1

    }

    vertex_pos <- rep(label_pos_exon, nrow(gv))

    if (is.null(label_pos_junction)) {

        edge_pos <- ifelse(edge_short, 1, 3)

    } else {

        edge_pos <- rep(label_pos_junction, nrow(gd))

    }

    if (is.null(curvature)) {

        edge_curved <- as.numeric(!edge_short)

    } else {

        edge_curved <- rep(curvature, nrow(gd))

    }

    pin <- par()$pin
    asp <- pin[2] / pin[1]
    edge_curved <- edge_curved / asp

    vertex_frame_color <- switch(border, fill = gv$color, none = NA, border)
    vertex_shape <- ifelse(!is.na(gv$type), "rectangle", "none")
    edge_lty <- rep(1, nrow(gd))

    if (include == "exons") {

        edge_lty <- rep(0, nrow(gd))

    } else if (include == "junctions") {

        vertex_shape <- rep("none", nrow(gv))

    }

    plot(g,
        add = TRUE,
        rescale = FALSE,
        layout = cbind(exon_coordinates$vertex_x, vertex_y),
        vertex.size = exon_coordinates$vertex_width * 100,
        vertex.size2 = 2 * ypos[2] * 100,
        vertex.color = gv$color,
        vertex.frame.color = vertex_frame_color,
        vertex.shape = vertex_shape,
        vertex.label = NA,
        edge.color = gd$color,
        edge.width = 2 * par()$cex,
        edge.lty = edge_lty,
        edge.curved = edge_curved,
        edge.label = NA,
        edge.arrow.mode = rep(0, nrow(gd)))

    df <- getGraphInfo(g, exon_coordinates$vertex_x, vertex_y,
        exon_coordinates$vertex_width, vertex_pos, edge_curved,
        edge_pos, color_labels, tx_view)

    if (label != "none") {

        i <- seq_len(nrow(df))

        if (include == "exons") {

            i <- which(df$type == "E")

        } else if (include == "junctions") {

            i <- which(df$type == "J")

        }

        if (length(i) > 0) {

            df_tmp <- df[i, ]

            text(x = df_tmp$x,
                y = df_tmp$y + c(E = -1, J = 0)[df_tmp$type] * ypos[2],
                labels = df_tmp[, label],
                pos = df_tmp$pos,
                offset = c(E = 0.5, J = 0.5)[df_tmp$type],
                col = df_tmp$color)

        }

    }

    if (tx_view) {

        df_collapsed <- unique(df[, c("y", "tx_name")])
        mtext(side = 4, at = df_collapsed$y, text = df_collapsed$tx_name,
            cex = par()$cex, line = 0.5, las = 1)

    }

    cols <- c("id", "name", "type", "featureID", "label", "color")
    cols <- cols[cols %in% names(df)]
    df <- df[, cols]

    return(df)

}

getExonCoordinates <- function(g, toscale)
{

    ## data frames of nodes
    gv <- nodes(g)

    ## exon and gene coordinates
    gr_exon <- co2gr(gv$coordinates)
    gr_gene <- range(gr_exon)

    if (length(gr_gene) > 1) {

        stop("features must be on the same chromosome and strand")

    }

    ## exon coordinates relative to gene locus
    ir_exon <- ranges(mapToTranscripts(gr_exon, split(gr_gene, 1),
        ignore.strand = FALSE))

    ## exonic regions
    ir_exonic <- reduce(ir_exon)

    ## exon x-coordinates and widths
    exon_coordinates <- xcoordinates(ir_exon, ir_exonic, toscale)
    exon_width <- exon_coordinates$w
    exon_x <- exon_coordinates$x

    coordinates <- list(
        gene_chrom = as.character(seqnames(gr_gene)),
        gene_strand = as.character(strand(gr_gene)),
        exon_start = start(gr_exon),
        exon_end = end(gr_exon),
        vertex_x = exon_x,
        vertex_width = exon_width
    )

    return(coordinates)

}

plotTrackScore <- function(exon_coordinates, score, color, ylim, ypos,
    nbin, summary, label)
{

    if (is(score, "RleList")) score <- score[[exon_coordinates$gene_chrom]]

    tmp <- lapply(seq_along(exon_coordinates$vertex_x),
        scorePerExon, exon_coordinates, score)

    bp_x_start <- do.call(c, lapply(tmp, "[[", "x_start"))
    bp_x_end <- do.call(c, lapply(tmp, "[[", "x_end"))
    bp_y <- do.call(c, lapply(tmp, "[[", "y"))

    excl <- which(bp_x_start == bp_x_end)

    if (length(excl) > 0) {

        bp_x_start <- bp_x_start[-excl]
        bp_x_end <- bp_x_end[-excl]
        bp_y <- bp_y[-excl]

    }

    bin_size <- 2 / nbin
    bin_breaks <- seq(-1, 1, length.out = nbin + 1)

    bp_bin <- IRanges(
        as.integer(cut(bp_x_start, bin_breaks, right = FALSE)),
        as.integer(cut(bp_x_end, bin_breaks, right = TRUE)))

    bin_y <- tapply(bp_y[rep(seq_along(bp_bin), width(bp_bin))],
        as.integer(bp_bin), summary)

    if (is.null(ylim)) {

        ymin <- floor(min(bin_y, na.rm = TRUE, finite = TRUE))
        ymax <- ceiling(max(bin_y, na.rm = TRUE, finite = TRUE))

        if (ymin >= 0) {

            ylim <- c(0, ymax)

        } else {

            ylim <- c(ymin, ymax)

        }

    }

    bin_y[bin_y < ylim[1]] <- ylim[1]
    bin_y[bin_y > ylim[2]] <- ylim[2]

    ypos_2 <- c(-1 + 2 * ypos[1] - ypos[2], 2 * ypos[2])

    bin_y <- bin_y - ylim[1]
    bin_y <- bin_y / diff(range(ylim)) * ypos_2[2]
    bin_y <- bin_y + ypos_2[1]

    px <- rep(bin_breaks, rep(2, length(bin_breaks)))
    py <- c(ypos_2[1], rep(bin_y, rep(2, length(bin_y))), ypos_2[1])
    polygon(px, py, border = NA, col = color)

    axis(side = 2, at = c(ypos_2[1], ypos_2[1] + ypos_2[2]), labels = ylim,
        mgp = c(3, 0.5, 0.5), tcl = -0.25, las = 1)

    mtext(side = 2, at = ypos_2[1] + 0.5 * ypos_2[2],
        cex = par()$cex, text = label, line = 3, las = 1)

}

scorePerExon <- function(i, co, score)
{

    ir_exon <- IRanges(co$exon_start[i], co$exon_end[i])
    y <- as.numeric(score[ir_exon])
    if (co$gene_strand == "-") { y <- rev(y) }

    ex_start <- co$vertex_x[i] - 0.5 * co$vertex_width[i]
    ex_end <- co$vertex_x[i] + 0.5 * co$vertex_width[i]

    d <- (ex_end - ex_start) / length(y)

    ex_bp_x_start <- seq(ex_start, ex_end - d, length.out = length(y))
    ex_bp_x_end <- seq(ex_start + d, ex_end, length.out = length(y))
    ex_bp_y <- y

    if (i < length(co$vertex_x) && co$exon_start[i + 1] - co$exon_end[i] > 1) {

        ir_intron <- IRanges(co$exon_end[i], co$exon_start[i + 1]) - 1
        y <- as.numeric(score[ir_intron])

        if (co$gene_strand == "-") {

          y <- rev(y)

        }

        if (co$gene_strand == "+") {

            in_start <- co$vertex_x[i] + 0.5 * co$vertex_width[i]
            in_end <- co$vertex_x[i + 1] - 0.5 * co$vertex_width[i + 1]

        } else {

            in_start <- co$vertex_x[i + 1] + 0.5 * co$vertex_width[i + 1]
            in_end <- co$vertex_x[i] - 0.5 * co$vertex_width[i]

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

plotTrackRanges <- function(exon_coordinates, toscale, score, color, ypos)
{

    gr_exon <- GRanges(exon_coordinates$gene_chrom,
        IRanges(exon_coordinates$exon_start,
        exon_coordinates$exon_end), exon_coordinates$gene_strand)
    gr_gene <- range(gr_exon)

    score <- score[which(unlist(range(score)) %over% gr_gene)]
    score <- endoapply(score, restrict, start(gr_gene), end(gr_gene))

    ypos_2 <- c(-1 + 2 * ypos[1] - ypos[2], 2 * ypos[2])

    n <- length(score)
    h <- ypos_2[2] / (2 * n - 1)
    y <- seq(ypos_2[1] + ypos_2[2] - 0.5 * h, ypos_2[1] + 0.5 * h,
        length.out = n)

    x <- mapply(plotRanges, score, y,
        MoreArgs = list(exon_coordinates = exon_coordinates,
            toscale = toscale, h = h, color = color))

    if (!is.null(names(score))) text(x, y, names(score), pos = 1)

}

plotRanges <- function(exon_coordinates, toscale, score, color, y, h)
{

    gr_exon <- GRanges(exon_coordinates$gene_chrom,
        IRanges(exon_coordinates$exon_start,
        exon_coordinates$exon_end), exon_coordinates$gene_strand)
    gr_gene <- range(gr_exon)

    ir_score <- ranges(mapToTranscripts(score, split(gr_gene, 1),
        ignore.strand = FALSE))
    ir_exon <- ranges(mapToTranscripts(gr_exon, split(gr_gene, 1),
        ignore.strand = FALSE))
    ir_exonic <- reduce(ir_exon)

    coords <- xcoordinates(ir_score, ir_exonic, toscale)

    xleft <- coords$x - 0.5 * coords$w
    xright <- coords$x + 0.5 * coords$w
    ybottom <- y - 0.5 * h
    ytop <- y + 0.5 * h

    mapply(rect, xleft = xleft, xright = xright,
        MoreArgs = list(ybottom = ybottom, ytop = ytop,
            col = color, border = NA))

    x <- 0.5 * (min(xleft) + max(xright))

    return(x)

}

getGraphInfo <- function(g, x_exon, y_exon, w_exon, vertex_pos,
    edge_curved, edge_pos, color_labels, tx_view)
{

    ## data frames of nodes and edges
    gv <- nodes(g)
    gd <- edges(g)

    i_from <- match(gd$from, gv$name)
    i_to <- match(gd$to, gv$name)

    ## coordinates for junction labels
    x_junc <- 0.5 * (x_exon[i_from] + 0.5 * w_exon[i_from] +
        x_exon[i_to] - 0.5 * w_exon[i_to])
    w_junc <- x_exon[i_to] - 0.5 * w_exon[i_to] -
        (x_exon[i_from] + 0.5 * w_exon[i_from])
    y_junc <- 0.5 * (y_exon[i_from] + y_exon[i_to]) +
        edge_curved * w_junc * 1/3

    ## output
    i_exon <- which(!is.na(gv$type))
    i_exon <- i_exon[order(x_exon[i_exon])]
    i_junc <- order(x_exon[i_from], w_junc)

    co_E <- gv$coordinates[i_exon]
    co_J <- gd$coordinates[i_junc]

    df <- data.frame(
        id = c(
            paste0(rep("E", length(co_E)),
                as.integer(factor(co_E, levels = unique(co_E)))),
            paste0(rep("J", length(co_J)),
                as.integer(factor(co_J, levels = unique(co_J))))),
        name = paste0(c(gv$type[i_exon], gd$type[i_junc]), ":",
            c(gv$coordinates[i_exon], gd$coordinates[i_junc])),
        type = c(rep("E", length(co_E)), rep("J", length(co_J))),
        featureID = c(gv$featureID[i_exon], gd$featureID[i_junc]),
        x = c(x_exon[i_exon], x_junc[i_junc]),
        y = c(y_exon[i_exon], y_junc[i_junc]),
        pos = c(vertex_pos[i_exon], edge_pos[i_junc]))

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
    ir_intronic <- IRanges(end(ir_exonic)[-n_exonic],
        start(ir_exonic)[-1]) - 1
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

    o <- order(c(ir_exonic, ir_intronic))
    ir_block <- c(ir_exonic, ir_intronic)[o]
    w_block <- c(w_exonic, w_intronic)[o]
    x_block <- c(0, cumsum(w_block[-length(w_block)])) - 1

    ol <- findOverlaps(ir, ir_block)
    ol_min <- as.integer(tapply(subjectHits(ol), queryHits(ol), min))
    ol_max <- as.integer(tapply(subjectHits(ol), queryHits(ol), max))

    x_1 <- x_block[ol_min] + (start(ir) - start(ir_block)[ol_min]) /
        width(ir_block)[ol_min] * w_block[ol_min]
    x_2 <- x_block[ol_max] + (end(ir) - start(ir_block)[ol_max] + 1) /
        width(ir_block)[ol_max] * w_block[ol_max]

    out <- list("x" = x_1 + 0.5 * (x_2 - x_1), "w" = x_2 - x_1)

    return(out)

}

setFeatureColors <- function(features, color, color_novel, alpha)
{

    if (is.null(mcols(features)$color)) {

        features_color <- rep(color, length(features))

        if (!is.null(color_novel) && !is.null(txName(features))) {

            txName <- txName(features)
            i_novel <- which(elementNROWS(txName) == 0)
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

##' Plot splice graph and heatmap of expression values.
##'
##' @title Plot splice graph and heatmap of expression values
##' @inheritParams plotSpliceGraph
##' @param x \code{SGFeatureCounts} object
##' @param cex Scale parameter for feature labels and annotation
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
##' @param heightPanels Numeric vector of length two indicating height of
##'   the top and bottom panels.
##' @param ... further parameters passed to \code{plotSpliceGraph}
##' @return \code{data.frame} with information on exon bins and
##'   splice junctions included in the splice graph
##' @examples
##' \dontrun{
##' sgfc_annotated <- annotate(sgfc_pred, txf_ann)
##' plotFeatures(sgfc_annotated)
##' }
##' NULL
##' @author Leonard Goldstein

plotFeatures <- function(x, geneID = NULL, geneName = NULL,
    which = NULL, tx_view = FALSE, cex = 1,
    assay = "FPKM", include = c("junctions", "exons", "both"),
    transform = function(x) { log2(x + 1) }, Rowv = NULL,
    distfun = dist, hclustfun = hclust, margin = 0.2,
    RowSideColors = NULL, square = FALSE, cexRow = 1, cexCol = 1,
    labRow = colnames(x), col = colorRampPalette(c("black", "gold"))(256),
    zlim = NULL, heightPanels = c(1, 2), ...)
{

    include <- match.arg(include)

    if (!is(x, "SGFeatureCounts")) {

        stop("x must be an SGFeatureCounts object")

    }

    x <- restrictFeatures(x, geneID = geneID, which = which,
        geneName = geneName)

    if (nrow(x) == 0) { return() }

    features <- rowRanges(x)

    if (!is.null(RowSideColors) && !is.list(RowSideColors)) {

        RowSideColors <- list(RowSideColors)

    }

    n_sample <- ncol(x)
    include_type <- switch(include, junctions = "J", exons = "E",
        both = c("E", "J"))
    n_feature <- length(which(type(features) %in% include_type))

    pars <- getLayoutParameters(n_sample, n_feature, margin, heightPanels,
        RowSideColors, square, tx_view)
    layout(pars$mat, heights = pars$hei, widths = pars$wid)
    par(cex = cex, mai = pars$mai)

    df <- plotSpliceGraph(x = features, label = "id", tx_view = tx_view,
        short_output = FALSE, ...)

    i_df <- switch(include,
        exons = which(df$type == "E"),
        junctions = which(df$type == "J"),
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

    df <- df[c("id", "name", "type", "featureID")]

    invisible(df)

}

##' Plot splice graph and heatmap of splice variant frequencies.
##'
##' @title Plot splice graph and heatmap of splice variant frequencies
##' @inheritParams plotSpliceGraph
##' @inheritParams plotFeatures
##' @param x \code{SGVariantCounts} object
##' @param cex Scale parameter for feature labels and annotation
##' @param transform Transformation applied to splice variant frequencies
##' @param expand_variants Experimental option - leave set to \code{FALSE}
##' @param ... further parameters passed to \code{plotSpliceGraph}
##' @return \code{data.frame} with information on exon bins and
##'   splice junctions included in the splice graph
##' @examples
##' \dontrun{
##' sgvc_annotated <- annotate(sgvc_pred, txf_ann)
##' plotVariants(sgvc_annotated)
##' }
##' NULL
##' @author Leonard Goldstein

plotVariants <- function(x, eventID = NULL, tx_view = FALSE,
    cex = 1, transform = function(x) { x }, Rowv = NULL,
    distfun = dist, hclustfun = hclust, margin = 0.2,
    RowSideColors = NULL, square = FALSE, cexRow = 1, cexCol = 1,
    labRow = colnames(x), col = colorRampPalette(c("black", "gold"))(256),
    zlim = c(0, 1), heightPanels = c(1, 2), expand_variants = FALSE, ...)
{

    if (!is(x, "SGVariantCounts")) {

        stop("x must be an SGVariantCounts object")

    }

    x <- updateObject(x, verbose = TRUE)

    x <- restrictFeatures(x, eventID = eventID,
        expand_variants = expand_variants)

    if (nrow(x) == 0) { return() }

    variant_not_quantifiable <- elementNROWS(featureID5p(x)) == 0 &
        elementNROWS(featureID3p(x)) == 0
    event_not_quantifiable <- tapply(variant_not_quantifiable, eventID(x),
        any)

    if (all(event_not_quantifiable)) {

        warning("none of the selected splice events could be quantified")

    } else if (any(event_not_quantifiable)) {

        warning("one or more splice events could not be quantified")

    }

    if (expand_variants) {

        x <- expandSGVariantCounts(x, eventID)

    }

    X <- transform(variantFreq(x))

    if (!is.null(RowSideColors) && !is.list(RowSideColors)) {

        RowSideColors <- list(RowSideColors)

    }

    exclude <- which(colSums(is.na(X)) == nrow(X))
    n_exclude <- length(exclude)

    if (n_exclude > 0) {

        labRow <- labRow[-exclude]
        x <- x[, -exclude]
        X <- X[, -exclude]
        RowSideColors <- lapply(RowSideColors, "[", -exclude)
        warning(paste("excluded", n_exclude, "samples with all values NA"))

    }

    n_sample <- ncol(x)
    n_feature <- nrow(x)

    pars <- getLayoutParameters(n_sample, n_feature, margin, heightPanels,
        RowSideColors, square, tx_view)
    layout(pars$mat, heights = pars$hei, widths = pars$wid)
    par(cex = cex, mai = pars$mai)

    df <- plotSpliceGraph(x = rowRanges(x), label = "label",
        tx_view = tx_view, ...)

    if (all(event_not_quantifiable) || n_sample == 0) {

        invisible(df)

    }

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

    id2variant <- split(togroup0(variants), featureID(unlist(variants)))
    id2variant <- unstrsplit(sort(unique(CharacterList(id2variant))), ",")
    features <- uniqueFeatures(unlist(variants))
    mcols(features)$label <- id2variant[match(
        featureID(features), names(id2variant))]

    return(features)

}

getLayoutParameters <- function(n_sample, n_feature, margin, heightPanels,
    RowSideColors, square, tx_view)
{

    n_RSC <- length(RowSideColors)

    mat <- rbind(
        c(1, 1, 1),
        c(0, 0, 0),
        c(0, 2, 0),
        c(0, 0, 0))

    for (i in seq_len(n_RSC)) {

        n <- max(mat) + 1
        mat <- cbind(c(1, 0, 0, 0), c(1, 0, n, 0), mat)

    }

    heightTopPanel <- heightPanels[1] / sum(heightPanels)
    heightBottomPanel <- heightPanels[2] / sum(heightPanels)

    mar_wid <- c(0.05, margin)
    mar_hei <- c(0.05, 0.15) * heightBottomPanel

    img_wid <- 1 - sum(mar_wid) - 0.025 * n_RSC * 2
    img_hei <- heightBottomPanel - sum(mar_hei)

    wid <- c(mar_wid[1], rep(0.025, n_RSC * 2), img_wid, mar_wid[2])
    hei <- c(heightTopPanel, mar_hei[1], img_hei, mar_hei[2])

    if (square) {

        padding <- getPadding(n_sample, n_feature, wid, hei)

        i_wid_mar <- c(1, length(wid))
        i_wid_img <- length(wid) - 1

        wid[i_wid_mar] <- wid[i_wid_mar] + 0.5 * padding[1]
        wid[i_wid_img] <- wid[i_wid_img] - padding[1]

        i_hei_mar <- c(2, 4)
        i_hei_img <- 3

        hei[i_hei_mar] <- hei[i_hei_mar] + 0.5 * padding[2]
        hei[i_hei_img] <- hei[i_hei_img] - padding[2]

    }

    ## add color key

    k <- ncol(mat)
    r <- nrow(mat)
    mat <- cbind(mat[, seq_len(k - 2)], mat[, k - 1], mat[, k - 1], mat[, k])
    wid <- c(wid[seq_len(k - 2)], wid[k - 1] * 2/3, wid[k - 1] * 1/3, wid[k])
    mat <- rbind(mat, 0, 0)
    mat[nrow(mat) - 1, ncol(mat) - 1] <- max(mat) + 1
    hei <- c(hei[seq_len(r - 1)], rep(0.05, 2) * heightBottomPanel,
        hei[r] - 0.1 * heightBottomPanel)

    ## splice graph margins and aspect ratio

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

    list(mat = mat, wid = wid, hei = hei, mai = mai)

}

getPadding <- function(n_sample, n_feature, wid, hei)
{

    ds <- dev.size()
    dev_wid <- ds[1]
    dev_hei <- ds[2]

    img_wid <- wid[length(wid) - 1]
    img_hei <- hei[3]

    img_wid_real <- img_wid / sum(wid) * dev_wid
    img_hei_real <- img_hei / sum(hei) * dev_hei

    cel_wid_real <- img_wid_real / n_feature
    cel_hei_real <- img_hei_real / n_sample

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

        mtext(side = 3, at = seq(from = 0, to = 1,
            length.out = length(labCol)), text = labCol,
            cex = par()$cex * cexCol, col = colLabCol,
            line = 0.25, las = 3)

    }

    if (!is.null(labRow)) {

        mtext(side = 4, at = seq(from = 0, to = 1, length.out = ncol(x)),
            text = labRow[j], cex = par()$cex * cexRow,
            line = 0.5, las = 1)

    }

    for (r in seq_along(RowSideColors)) {

        par(mai = c(0, 0, 0, 0))
        image(matrix(seq_len(ncol(x)), nrow = 1),
            col = RowSideColors[[r]], axes = FALSE)
        mtext(side = 3, text = names(RowSideColors)[r],
            cex = par()$cex * cexCol, line = 0.25, las = 3)

    }

    par(mai = c(0, 0, 0, 0))
    image(matrix(seq_along(col), ncol = 1), col = col, axes = FALSE)
    mtext(side = 1, at = c(0, 1), text = format(zlim, digits = 2),
        cex = par()$cex, line = 0.25, las = 1)

}

restrictFeatures <- function(x, geneID = NULL, eventID = NULL, which = NULL,
    geneName = NULL, expand_variants = FALSE)
{

    if (is(x, "Counts")) {

        y <- rowRanges(x)

    } else {

        y <- x

    }

    if (!is.null(geneID) || !is.null(eventID) || !is.null(geneName)) {

        if (!is.null(geneID)) {

            index <- which(geneID(y) %in% geneID)

        } else if (!is.null(eventID)) {

            index <- which(eventID(y) %in% eventID)

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
                    unlist(expandString(featureID(y)[index])), ",",
                        fixed = TRUE)))
                index <- unlist(lapply(featureIDs,
                    function (id) { grep(paste0("(^|,)", id, "(,|$)"),
                        featureID(y)) }))
                eventIDs <- unique(eventID(y)[index])
                index <- which(eventID(y) %in% eventIDs)

            }

        } else if (!is.null(geneName)) {

            if (is.null(geneName(y))) {

                stop("missing geneName")

            }

            if (length(geneName) > 1) {

                stop("geneName must be of length 1")

            }

            index <- which(any(geneName(y) == geneName))

            if (all(type(y)[index] %in% c("D", "A"))) index <- vector()
            
        }

        y <- y[index]

    } else if (!is.null(which)) {

        ol <- findOverlaps(y, which, type = "within")

        index <- which(
            type(y) == "E" & (y %over% which) |
            type(y) == "J" & (seq_along(y) %in% queryHits(ol)))

        y <- restrict(y[index], start = start(which), end = end(which))


    } else {

        return(x)

    }

    if (is(x, "Counts")) {

        x <- x[index, ]
        rowRanges(x) <- y
        return(x)

    } else {

        return(y)

    }

}

##' Plot read coverage and splice junction read counts for an individual
##' sample or averaged across samples.
##'
##' @title Plot read coverage and splice junction read counts
##' @inheritParams plotSpliceGraph
##' @param x \code{SGFeatureCounts} or \code{SGFeatures} object.
##'   If \code{x} is an \code{SGFeatureCounts} object that includes
##'   multiple samples, average coverage and splice junction counts
##'   are obtained.
##' @param sample_info Data frame with sample information.
##'   If \code{x} is an \code{SGFeatureCounts} object, sample information
##'   is obtained from \code{colData(x)}. If \code{sample_info} includes
##'   multiple samples, average coverage and splice junction counts
##'   are obtained.
##' @param sizefactor Numeric vector with length equal to the number of
##'   samples in \code{sample_info}. Used to scale coverages and splice
##'   junction counts before plotting, or before averaging across samples.
##'   Set to \code{NA} to disable scaling. If \code{NULL}, size factors
##'   are calculated as the number of bases sequenced (the product of library
##'   size and average number of bases sequenced per read or fragment),
##'   plotted coverages and splice junction counts are per 1 billion
##'   sequenced bases.
##' @param color Color used for plotting coverages
##' @param ylim Numeric vector of length two, determining y-axis range used
##'   for plotting coverages.
##' @param nbin Number of bins for plotting coverages
##' @param summary Function used to calculate per-bin coverage summaries
##' @param label Optional y-axis label
##' @param min_anchor Integer specifiying minimum anchor length
##' @param cores Number of cores available for parallel processing.
##' @return \code{data.frame} with information on splice junctions included
##'   in the splice graph
##' @examples
##' \dontrun{
##' par(mfrow = c(4, 1))
##' for (j in seq_len(4)) plotCoverage(sgfc_pred[, j])
##' }
##' NULL
##' @author Leonard Goldstein

plotCoverage <- function(x, geneID = NULL, geneName = NULL,
    eventID = NULL, which = NULL, sample_info = NULL, sizefactor = NA,
    toscale = c("exon", "none", "gene"), color = "darkblue",
    ylim = NULL, label = NULL, nbin = 200, summary = mean,
    curvature = 1, main = NULL, min_anchor = 1, cores = 1)
{

    toscale <- match.arg(toscale)

    x <- restrictFeatures(x, geneID, eventID, which, geneName)

    if (is(x, "SGFeatures") && !is.null(sample_info)) {

        sgf <- x
        sgfc <- getSGFeatureCounts(sample_info, sgf,
            min_anchor = min_anchor, cores = cores)

    } else if (is(x, "SGFeatureCounts")) {

        sample_info <- colData(x)
        sgf <- rowRanges(x)
        sgfc <- x

    } else {

        stop("either x must be an SGFeatureCounts object,
            or x must be an SGFeatures object and sample_info not NULL")

    }

    if (is.null(sizefactor)) {

        sizefactor <- calculateSizeFactor(sample_info)

    } else if (length(sizefactor) == 1 && is.na(sizefactor)) {

        sizefactor <- rep(1, nrow(sample_info))

    } else if (length(sizefactor) != nrow(sample_info)) {

        stop("sizefactor must have length equal to the number of samples")

    }

    scaled_cov <- getCoverage(sample_info, range(sgf), sizefactor, cores)
    average_cov <- Reduce("+", scaled_cov) / length(scaled_cov)

    scaled_counts <- sweep(counts(sgfc), 2, sizefactor, FUN = "/")
    average_counts <- rowSums(scaled_counts) / ncol(scaled_counts)

    sgf <- setFeatureColors(sgf, color, color, 1)
    sgf$label <- round(average_counts)
    g <- exonGraph(sgf, FALSE)
    exon_coordinates <- getExonCoordinates(g, toscale)

    plot(NA, xlim = c(-1, 1), ylim = c(-1, 1), xaxs = "i", yaxs = "i",
        axes = FALSE, xlab = NA, ylab = NA)
    plotTrackScore(exon_coordinates, average_cov, color, ylim, c(0.5, 1),
        nbin, summary, label)
    df <- plotExonGraph(g, exon_coordinates, c(0, 1), "junctions", "none",
        curvature, "label", FALSE, 1, NULL, FALSE, NA)
    text(x = 0, y = 0.95, labels = main, pos = 1, offset = 0, font = 2)

    invisible(df)

}
