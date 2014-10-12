feature2name <- function(features, collapse_terminal = FALSE)
{

    if (is(features, "Features")) {
    
        features_type <- type(features)

    } else {

        features_type <- mcols(features)$type

    }
    
    name <- rep(NA, length(features))
    
    if (collapse_terminal) {

        i <- which(features_type %in% c("J", "I", "U", "E"))
        name[i] <- paste0(features_type[i], ":", gr2co(features[i]))
        i <- which(features_type %in% c("F", "L"))
        start <- ifelse(features_type[i] == "F", FALSE, TRUE)
        name[i] <- paste0(features_type[i], ":",
            gr2co(flank(features[i], -1, start)))
        
    } else {
    
        i <- which(features_type %in% c("J", "I", "F", "L", "U", "E"))
        name[i] <- paste0(features_type[i], ":", gr2co(features[i]))
        
    }

    i <- which(features_type %in% c("D", "A"))
    name[i] <- paste0(features_type[i], ":", gr2pos(features[i]))

    return(name)

}

gr2co <- function(x)
{

    paste0(seqnames(x), ":", start(x), "-", end(x), ":", strand(x))
    
}

co2gr <- function(co)
{

    x <- strsplit(co, split = ":", fixed = TRUE)
    r <- strsplit(sapply(x, "[", 2), split = "-", fixed = TRUE)
    sn <- sapply(x, "[", 1)
    start <- as.integer(sapply(r, "[", 1))
    end <- as.integer(sapply(r, "[", 2))
    st <- sapply(x, "[", 3)
    GRanges(sn, IRanges(start, end), st)
    
}

gr2pos <- function(x)
{

    paste0(seqnames(x), ":", start(x), ":", strand(x))

}

pos2gr <- function(x)
{

    x <- strsplit(x, split = ":", fixed = TRUE)    
    sn <- sapply(x, "[", 1)
    pos <- as.integer(sapply(x, "[", 2))
    st <- sapply(x, "[", 3)    
    GRanges(sn, IRanges(pos, pos), st)

}

readGap <- function(file, paired_end, which = NULL)
{

    if (is.null(which)) {

        param <- ScanBamParam(tag = "XS")
        
    } else {

        param <- ScanBamParam(tag = "XS", which = reduce(which))

    }

    if (paired_end) {

        gap <- suppressWarnings(readGAlignmentPairs(file = file,
            param = param))
        gap <- propagateXS(gap)

    } else {

        gap <- suppressWarnings(readGAlignments(file = file, param = param))

    }

    gap <- filterGap(gap)
    
    strand(gap) <- XS2strand(mcols(gap)$XS)
    
    return(gap)

}

propagateXS <- function(gap)
{

    first_xs <- mcols(first(gap))$XS
    last_xs <- mcols(last(gap))$XS
    xs <- first_xs
    xs[is.na(xs)] <- last_xs[is.na(xs)]
    mcols(gap)$XS <- xs
    return(gap)

}

XS2strand <- function(xs)
{

    s <- xs
    s[is.na(s)|s == "?"] <- "*"
    return(s)
    
}

filterGap <- function(gap)
{

    if (is(gap, "GAlignments")) {

        exclude <- filterGa(gap)

    }    
    if (is(gap, "GAlignmentPairs")) {
        
        exclude <- filterGa(left(gap)) | filterGa(right(gap))

    }
    
    gap <- gap[!exclude]

    return(gap)

}

filterGa <- function(ga)
{

    grepl("(\\d+D\\d+N)|(\\d+N\\d+D)", cigar(ga))
    
}

dropMcols <- function(x)
{

    mcols(x) <- NULL
    return(x)
    
}

completeMcols <- function(x, include_counts, retain_coverage)
{

    mcol <- "type"
    
    if (include_counts) {

        mcol <- c(mcol, "N")

    }
    if (retain_coverage) {

        mcol <- c(mcol, "N_splicesite", "coverage")

    }
        
    for (m in setdiff(mcol, names(mcols(x)))) {

        if (m == "N") {

            mcols(x)[, m] <- NA_integer_

        } else if (m == "N_splicesite") {

            mcols(x)[, m] <- IntegerList(as.integer())

        } else if (m == "coverage") {

            mcols(x)[, m] <- RleList(as.integer(), compress = TRUE)

        }
        
    }
    
    mcols(x) <- mcols(x)[, mcol, drop = FALSE]
    names(mcols(x)) <- mcol        
    
    return(x)

}

getBamInfoPerSample <- function(file_bam, yieldSize = NULL)
{

    if (is.null(yieldSize)) {

        file <- BamFile(file_bam)
        
    } else {

        file <- BamFile(file_bam, yieldSize = yieldSize)

    }

    param <- ScanBamParam(what = c("qname", "qwidth", "mrnm", "isize"))
    
    bam <- scanBam(file = file, param = param)[[1]]
    paired_end <- !all(is.na(bam$mrnm))
    read_length <- median(bam$qwidth, na.rm = TRUE)

    if (paired_end) {
        
        isize <- bam$isize
        frag_length <- median(isize[which(isize > 0)], na.rm = TRUE)

    } else {

        frag_length <- NA

    }

    x <- data.frame(
       paired_end = paired_end,
       read_length = read_length,
       frag_length = frag_length)

    if (is.null(yieldSize)) {
        
        x$lib_size <- length(unique(bam$qname))
        
    } 

    return(x)
    
}

expandUnstrandedRanges <- function(x)
{

    i <- which(strand(x) == "*")

    if (length(i) > 0) {
    
        additional <- x[i]
        strand(additional) <- "-"
        strand(x)[i] <- "+"
        x <- c(x, additional)

    }

    return(x)

}

uniqueFeatures <- function(features)
{

    i_duplicated <- vector()
    
    for (type in levels(type(features))) {

        i_type <- which(type(features) == type)
        i <- i_type[which(duplicated(features[i_type]))]
        i_duplicated <- c(i_duplicated, i)

    }

    if (length(i_duplicated) > 0) {

        features <- features[-i_duplicated]

    }

    return(features)
    
}

##' Export features to BED format. Splice sites are not included.
##'
##' @title Export to BED format
##' @param features \code{TxFeatures} or \code{SGFeatures} object
##' @param file Character string specifying output file
##' @return \code{NULL}
##' @examples
##' \dontrun{
##'   exportFeatures(txf, "txf.bed")
##'   exportFeatures(sgf, "sgf.bed")
##' }
##' @author Leonard Goldstein

exportFeatures <- function(features, file)
{

    if (!is(features, "Features")) {

        stop("features must be a TxFeatures or SGFeatures object")

    }

    features <- asGRanges(features)
    
    i_splicesite <- which(mcols(features)$type %in% c("D", "A"))

    if (length(i_splicesite) > 0) {

        features <- features[-i_splicesite]

    }

    i_junction <- which(mcols(features)$type == "J")
    color <- mcols(features)$color
    mcols(features) <- NULL
    
    bed <- split(features, seq_along(features))

    if (length(i_junction) > 0) {
    
        bed[i_junction] <- psetdiff(
           split(features[i_junction], seq_along(i_junction)),
           split(features[i_junction] - 1, seq_along(i_junction)))

    }

    if (!is.null(color)) {

        itemRgb <- rgb(t(col2rgb(color)), maxColorValue = 255)
        mcols(bed)$itemRgb <- itemRgb
        
    }

    names(bed) <- feature2name(features)
    
    export(object = bed, con = file, format = "BED")

    return()
    
}

nextFrame <- function(f, w, prev = FALSE)
{

    if (is(f, "list") || is(f, "List")) { 
    
        f_unlisted <- unlist(f)
        w_unlisted <- w[togroup(f)]
        n_unlisted <- nextFrame(f_unlisted, w_unlisted, prev)
        n <- relist(n_unlisted, f)

    } else {

        if (prev) {

            n <- ifelse(f != -1, (f - w) %% 3, -1)

        } else {
            
            n <- ifelse(f != -1, (f + w) %% 3, -1)

        }

    }

    return(n)
    
}

collapseCharacterList <- function(x, f)
{

    if (!is(f, "factor")) {

        stop("f must be a factor")

    }

    x_unlisted <- setNames(unlist(x), NULL)
    f_unlisted <- f[togroup(x)]
    y <- CharacterList(split(x_unlisted, f_unlisted))
    y <- unique(y)

    return(y)

}

asGRanges <- function(from)
{

    granges(from, TRUE)
    
}

asGRangesList <- function(from)
{

    as(from, "GRangesList")
    
}

reorderFeatures <- function(x)
{

    x_names <- names(x)
    x_mc <- mcols(x)
    features <- unlist(x, use.names = FALSE)
    features_x <- togroup(x)
    i_pos <- which(strand(features) == "+" | strand(features) == "*")
    i_neg <- which(strand(features) == "-")
    i_pos <- i_pos[order(features[i_pos])]
    i_neg <- i_neg[order(features[i_neg], decreasing = TRUE)]
    i_all <- c(i_pos, i_neg)
    x <- split(features[i_all], features_x[i_all])
    names(x) <- x_names
    mcols(x) <- x_mc
    
    return(x)

}

pfirst <- function(x, use_names = FALSE)
{

    unlist(phead(x, 1), use.names = use_names)

}

plast <- function(x, use_names = FALSE)
{

    unlist(ptail(x, 1), use.names = use_names)

}
