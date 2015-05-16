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

co2str <- function(seqlevel, start, end, strand)
{

    paste0(seqlevel, ":", start, "-", end, ":", strand)
    
}

gr2co <- function(x)
{

    if (length(x) == 0) {

        return()

    } else {
      
        co2str(seqnames(x), start(x), end(x), strand(x))

    }
    
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

pos2str <- function(seqlevel, position, strand)
{

    paste0(seqlevel, ":", position, ":", strand)

}

gr2pos <- function(x)
{

    if (length(x) == 0) {

        return()

    } else {
      
        pos2str(seqnames(x), start(x), strand(x))

    }
    
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

    ## the following flags are set by functions
    ## readGAlignments and readGAlignmentPairs
    ## - isUnmappedQuery
    ## - isPaired
    ## - hasUnmappedMate

    flag <- scanBamFlag(isSecondaryAlignment = FALSE)
    param <- ScanBamParam(flag = flag, tag = "XS")
    
    if (!is.null(which)) {

        bamWhich(param) <- reduce(which)

    }

    if (paired_end) {
      
        gap <- suppressWarnings(readGAlignmentPairs(file = file,
            param = param))

        ## scanBam workaround start
        ## bamWhat(param) <- c("flag", "mrnm", "mpos")
        ## ga <- readGAlignments(file = file, use.names = TRUE, param = param)
        ## gap <- makeGAlignmentPairs(ga, use.names = TRUE, use.mcols = TRUE)
        ## names(gap) <- NULL
        ## scanBam workaround end
      
        gap <- propagateXS(gap)

    } else {

        gap <- suppressWarnings(readGAlignments(file = file, param = param))
        
    }

    gap <- filterGap(gap)
    
    mcols(gap)$strand <- XS2strand(mcols(gap)$XS)

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
        
        exclude <- filterGa(first(gap)) | filterGa(last(gap))

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

getBamInfoPerSample <- function(file_bam, yieldSize, sample_name)
{

    if (is.null(yieldSize)) file <- BamFile(file_bam)
    else file <- BamFile(file_bam, yieldSize = yieldSize)

    flag <- scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE)
    what <- c("qname", "flag", "qwidth", "isize")
    param <- ScanBamParam(flag = flag, what = what, tag = "XS")
    bam <- scanBam(file = file, param = param)[[1]]

    XS <- !is.null(bam$tag$XS)
    paired_end <- any(bamFlagTest(bam$flag, "isPaired"))
    read_length <- median(bam$qwidth, na.rm = TRUE)

    if (paired_end) {
        
        isize <- bam$isize
        frag_length <- median(isize[which(isize > 0)], na.rm = TRUE)

    } else {

        frag_length <- NA_real_

    }

    x <- data.frame(
       sample_name = sample_name,
       file_bam = file_bam,
       XS = XS,
       paired_end = paired_end,
       read_length = read_length,
       frag_length = frag_length)

    if (is.null(yieldSize)) {
        
        x$lib_size <- length(unique(bam$qname))
        
    } 

    generateCompleteMessage(sample_name)
    
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
##' exportFeatures(txf, "txf.bed")
##' exportFeatures(sgf, "sgf.bed")
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

rbindListOfDFs <- function(x, cores)
{

    k <- names(x[[1]])

    df <- vector("list", length(k))

    for (j in seq_along(k)) {

        df[[j]] <- do.call(c, mclapply(x, "[[", j, mc.cores = cores))

    }

    names(df) <- k
    
    df <- DataFrame(df, check.names = FALSE)

    return(df)
  
}

checkApplyResultsForErrors <- function(out, fun_name, items)
{
  
    failed <- sapply(out, is, "try-error")

    if (any(failed)) {
  
        i <- which(failed)
        err <- makeErrorMessage(fun_name, items[i], out[i])
        err <- paste0("\n", err)
        stop(err, call. = FALSE)
        
    }

}

makeErrorMessage <- function(fun_name, items, msgs)
{

    msg <- paste0("Error in ", fun_name, " for ", items, ":", "\n", msgs)
    msg <- paste(msg, collapse = "\n")

    return(msg)
 
}

makeWarningMessage <- function(fun_name, items, msgs)
{

    msg <- paste0("Warning in ", fun_name, " for ", items, ":", "\n", msgs)
    msg <- paste(msg, collapse = "\n")

    return(msg)
 
}

makeCompleteMessage <- function(item) {

    paste(item, "complete.")

}

generateWarningMessage <- function(fun_name, item, msg)
{

    message(makeWarningMessage(fun_name, item, msg))

}

generateCompleteMessage <- function(item)
{

    message(makeCompleteMessage(item))

}

getCoverage <- function(sample_info, which, sizefactor, cores)
{

  if (!is(which, "GRanges") || length(which) > 1) {

      stop("which must be a GRanges object of length 1")

  }
  
  list_cov <- mcmapply(
      getCoveragePerSample,
      file_bam = sample_info$file_bam,
      paired_end = sample_info$paired_end,
      sizefactor = sizefactor,
      MoreArgs = list(which = which),
      SIMPLIFY = FALSE,
      USE.NAMES = FALSE,
      mc.preschedule = FALSE,
      mc.cores = cores) 

  return(list_cov)

}

getCoveragePerSample <- function(file_bam, paired_end, sizefactor, which)
{

  st <- as.character(strand(which))  
  gap <- readGap(file_bam, paired_end, which)
  gap <- gap[mcols(gap)$strand %in% c(st, "*")]
  irl <- ranges(grglist(gap, drop.D.ranges = TRUE))
  ir <- unlist(reduce(irl))
  cov <- coverage(ir, width = end(which)) / sizefactor
  
  return(cov)
  
}

calculateSizeFactor <- function(sample_info)
{
        
    E <- rep(NA, nrow(sample_info))

    i_PE <- which(sample_info$paired_end)
    
    if (length(i_PE) > 0) {

        R_PE <- sample_info$read_length[i_PE]
        F_PE <- sample_info$frag_length[i_PE]
        I_PE <- F_PE - 2 * R_PE
        E[i_PE] <- F_PE - pmax(I_PE, 0)
    
    }

    i_SE <- which(!sample_info$paired_end)
        
    if (length(i_SE) > 0) {

        E[i_SE] <- sample_info$read_length[i_SE]
    
    }

    sizefactor <- sample_info$lib_size * E * 1e-9
    
    return(sizefactor)

}
