##' The effect of each splice variant is assessed with respect to individual
##' protein-coding transcripts. 
##' 
##' @title Predict the effect of splice variants on protein-coding transcripts
##' @param sgv \code{SGVariants} object
##' @param tx A \code{TxDb} object, or a \code{GRangesList} of
##'   exons grouped by transcripts with metadata columns \code{cdsStart} and
##'   \code{cdsEnd} indicating genomic start and end of the coding sequence,
##'   respectively. By convention, cdsStart < cdsEnd, regardless of strand. 
##' @param genome \code{BSgenome} object
##' @param annotate Logical indicating whether \code{sgv} should be annotated
##'   with respect to transcripts provided in argument \code{tx}.
##'   Can be set to \code{FALSE} if variants are already annotated.
##' @param summarize Logical indicating whether results should be
##'   summarized per variant
##' @param cores Number of cores available for parallel processing
##' @return For \code{summarize = FALSE} a \code{data.frame} with rows
##'   corresponding to a variant-transcript pair. The \code{data.frame}
##'   includes columns for variant identifier, transcript name, type of
##'   alteration, protein sequences for the reference transcript and the
##'   transcript variant, as well as coordinates of the variant in the 
##'   protein sequences. Start and end coordinates are 0- and 1-based,
##'   respectively, to allow for the specification of deletions.
##'   For \code{summarize = TRUE} a character vector matching argument
##'   \code{sgv} with comma-separated predicted alterations for individual
##'   transcripts.
##' @examples
##' require(BSgenome.Hsapiens.UCSC.hg19)
##' seqlevelsStyle(Hsapiens) <- "NCBI"
##' predictVariantConsequences(sgv_pred, tx, Hsapiens)
##' @author Leonard Goldstein

predictVariantConsequences <- function(sgv, tx, genome, annotate = TRUE,
    summarize = TRUE, cores = 1)
{

    if (is(tx, "TxDb")) {

        tx <- convertToTranscripts(tx)
      
    } else if (is(tx, "GRangesList")) {

        checkTranscripts(tx)

        if (is.null(names(tx))) {

            names(tx) <- seq_along(tx)

        }

    } else {

        stop("tx must be a TxDb object or GRangesList of exons
            grouped by transcripts")

    }

    common_seqlevels <- intersect(seqlevels(sgv), seqlevels(tx))
    sgv <- keepSeqlevels(sgv, common_seqlevels)
    tx <- keepSeqlevels(tx, common_seqlevels)

    if (annotate) {
    
        sgv <- annotate(sgv, convertToTxFeatures(tx))

    }

    sgv_vid <- variantID(sgv)

    eid_tx <- split(unlist(txName(sgv)), eventID(sgv)[togroup(txName(sgv))])
    
    excl <- grep("(", type(sgv), fixed = TRUE)

    if (length(excl) > 0) {

        sgv <- sgv[-excl]

    }
    
    sgv_tx <- eid_tx[match(eventID(sgv), names(eid_tx))]
    sgv_tx <- mcmapply(
        setdiff,
        sgv_tx,
        txName(sgv),
        SIMPLIFY = FALSE,
        mc.cores = cores)
    
    expanded_sgv <- sgv[rep(seq_along(sgv), elementLengths(sgv_tx))]
    expanded_tx <- tx[match(unlist(sgv_tx), names(tx))]

    list_expanded_sgv <- split(expanded_sgv, seq_along(expanded_sgv))
    list_expanded_tx <- split(expanded_tx, seq_along(expanded_tx))
    
    res <- do.call(rbind, mcmapply(
        predictVariantConsequencesPerVariantAndTranscript,
        list_expanded_sgv,
        list_expanded_tx,
        MoreArgs = list(genome = genome),
        SIMPLIFY = FALSE,
        mc.cores = cores))
    rownames(res) <- NULL
    
    if (summarize) {

        res$consequence <- paste0(res$ref_name, ":", res$alt)
        vid_consequence <- unstrsplit(split(res$consequence, res$var_id), ",")
        res <- vid_consequence[match(sgv_vid, names(vid_consequence))]
        names(res) <- NULL
        
    }
    
    return(res)
    
}

predictVariantConsequencesPerVariantAndTranscript <- function(sgv, ref, genome)
{

    var <- getTranscriptVariant(sgv, ref, genome)

    n <- length(var)
    
    res <- data.frame(
        var_id = rep(variantID(sgv), n),
        ref_name = rep(names(ref), n),
        alt = alt(var),
        stringsAsFactors = FALSE)

    list_var <- split(var, seq_along(var))
    
    regions <- do.call(rbind, lapply(
        list_var,
        getVariantRegions,
        ref = ref,
        sgv = sgv,
        genome = genome))

    res <- cbind(res, regions)

    return(res)
    
}

getTranscriptVariant <- function(sgv, ref, genome) 
{

    ref_loc <- range(ref[[1]])
    ref_cds <- range(cds(ref))
    
    utr_5p <- range(c(flank(ref_loc, -1, TRUE), flank(ref_cds, 1, TRUE)))
    utr_3p <- range(c(flank(ref_cds, 1, FALSE), flank(ref_loc, -1, FALSE)))

    event <- getEventLocation(ref, sgv)
    event_start <- flank(event, -1, TRUE)
    event_end <- flank(event, -1, FALSE)
    
    del <- intersect(ref[[1]], event)
    ins <- granges(sgv[[1]][type(sgv[[1]]) == "E"])
    var <- GRangesList(reduce(c(setdiff(ref[[1]], del), ins)))

    if (event_start %over% utr_5p && event_end %over% utr_5p) {

        var <- findCDS(var, NA_integer_, cdsEnd(ref), genome)

        if (identical(cds(var), cds(ref))) {

            alt(var) <- "5p_UTR_variant"

        } else {

            alt(var) <- "CDS_upstream_start"

        }
        
    } else if (event_start %over% utr_5p && event_end %over% ref_cds) {

        var <- findCDS(var, NA_integer_, cdsEnd(ref), genome)

        if (identical(cds(var), cds(ref))) {

            alt(var) <- "CDS_unaffected"

        } else if (is.na(cdsStart(var))) {

            alt(var) <- "CDS_start_lost"

        } else {

            alt(var) <- "CDS_5p_alteration"

        }
        
    } else if (event_start %over% utr_5p && event_end %over% utr_3p) {

        cdsStart(var) <- NA_integer_
        cdsEnd(var) <- NA_integer_
        alt(var) <- "CDS_lost"

    } else if (event_start %over% ref_cds && event_end %over% ref_cds) {
      
        var_3p <- findCDS(var, cdsStart(ref), NA_integer_, genome)
        var_5p <- findCDS(var, NA_integer_, cdsEnd(ref), genome)

        del_w <- sum(width(del))
        ins_w <- sum(width(ins))

        if (del_w > 0 && ins_w > 0) {
          
            alt_tx <- "alteration"

        } else if (del_w > 0) {

            alt_tx <- "deletion"

        } else if (ins_w > 0) {

            alt_tx <- "insertion"

        }

        if (del_w %% 3 == ins_w %% 3) {

            alt_cds <- "in-frame"

        } else {

            alt_cds <- "frame-shift"

        }
        
        prefix <- paste0("CDS_", alt_tx, "_", alt_cds)

        if (alt_cds == "in-frame") {

            if (cdsEnd(var_3p) == cdsEnd(ref)) {

                alt(var_3p) <- prefix
                var <- var_3p
                
            } else {

                alt(var_3p) <- paste0(prefix, "_premature_stop")
                var <- var_3p

                if (!is.na(cdsStart(var_5p))) {
                 
                    alt(var_5p) <- paste0(prefix, "_5p_alteration")
                    var <- c(var, var_5p)

                }
                
            }

        } else if (alt_cds == "frame-shift") {

            if (is.na(cdsEnd(var_3p)) && is.na(cdsStart(var_5p))) {

                alt(var_3p) <- paste0(prefix, "_CDS_lost")
                var <- var_3p

            } else if (!is.na(cdsEnd(var_3p)) && is.na(cdsStart(var_3p))) {

                alt(var_3p) <- paste0(prefix, "_3p_alteration")
                var <- var_3p

            } else if (is.na(cdsEnd(var_3p)) && !is.na(cdsStart(var_3p))) {

                alt(var_5p) <- paste0(prefix, "_5p_alteration")
                var <- var_5p

            } else if (!is.na(cdsEnd(var_3p)) && !is.na(cdsStart(var_3p))) {

                alt(var_3p) <- paste0(prefix, "_3p_alteration")
                alt(var_5p) <- paste0(prefix, "_5p_alteration")
                var <- c(var_3p, var_5p)
            
            }

        }
        
    } else if (event_start %over% ref_cds && event_end %over% utr_3p) {

        var <- findCDS(var, cdsStart(ref), NA_integer_, genome)

        if (identical(cds(var), cds(ref))) {

            alt(var) <- "CDS_unaffected"

        } else if (is.na(cdsEnd(var))) {

            alt(var) <- "CDS_stop_lost"

        } else {

            alt(var) <- "CDS_3p_alteration"

        }
      
    } else if (event_start %over% utr_3p && event_end %over% utr_3p) {

        cdsStart(var) <- cdsStart(ref)
        cdsEnd(var) <- cdsEnd(ref)
        alt(var) <- "3p_UTR_variant"
        
    }
    
    return(var)
    
}

findCDS <- function(tx, cdsStart, cdsEnd, genome)
{

    start_codons <- c("ATG")
    stop_codons <- c("TAG", "TAA", "TGA")        

    chrom <- as.character(seqnames(tx[[1]][1]))
    strand <- as.character(strand(tx[[1]][1]))

    names(tx) <- "1"
    tx[[1]] <- sort(tx[[1]], decreasing = (strand == "-"))
    tx_seq <- as.character(do.call(c, getSeq(genome, tx[[1]])))

    if (is.na(cdsStart)) {
      
        cdsEnd_gr <- GRanges(chrom, IRanges(cdsEnd, cdsEnd), strand)
        cdsEnd_tx <- start(mapToTranscripts(cdsEnd_gr, tx))
        p <- seq(from = cdsEnd_tx %% 3 + 1, to = cdsEnd_tx - 2, by = 3)
        codons <- mapply(substr, p, p + 2, MoreArgs = list(x = tx_seq))
        i_start <- which(codons %in% start_codons)
        i_stop <- which(codons[-length(codons)] %in% stop_codons)

        if (length(i_start) > 0 && length(i_stop) > 0) {
        
            i_start <- i_start[i_start > max(i_stop)]

        }
        
        if (length(i_start) > 0) {
              
            cdsStart_tx <- p[min(i_start)]
            x <- GRanges(1, IRanges(cdsStart_tx, cdsStart_tx), "*")
            cdsStart_gr <- mapFromTranscripts(x, tx)
            cdsStart <- start(cdsStart_gr)
            
        } else {

            cdsStart <- NA_integer_
            cdsEnd <- NA_integer_
                
        }
                                
    } else {

        cdsStart_gr <- GRanges(chrom, IRanges(cdsStart, cdsStart), strand)
        cdsStart_tx <- start(mapToTranscripts(cdsStart_gr, tx))
        p <- seq(from = cdsStart_tx, to = nchar(tx_seq), by = 3)
        codons <- mapply(substr, p, p + 2, MoreArgs = list(x = tx_seq))
        i_stop <- which(codons %in% stop_codons)
        
        if (length(i_stop) > 0) {
          
            cdsEnd_tx <- (p + 2)[min(i_stop)]
            x <- GRanges(1, IRanges(cdsEnd_tx, cdsEnd_tx), "*")
            cdsEnd_gr <- mapFromTranscripts(x, tx)
            cdsEnd <- start(cdsEnd_gr)
            
        } else {

            cdsStart <- NA_integer_
            cdsEnd <- NA_integer_
            
        }
        
    }

    cdsStart(tx) <- cdsStart
    cdsEnd(tx) <- cdsEnd
    
    return(tx)

}

getEventLocation <- function(x, sgv)
{

    start_type <- substr(from(sgv), 1, 1)
    end_type <- substr(to(sgv), 1, 1)
    
    if (start_type == "D") {

        start_gr <- pos2gr(sub("^D:", "", from(sgv)))
        start_gr <- flank(start_gr, 1, FALSE)
      
    } else if (start_type == "S") {

        start_gr <- flank(range(x[[1]]), -1, TRUE)
      
    }

    if (end_type == "A") {
      
        end_gr <- pos2gr(sub("^A:", "", to(sgv)))
        end_gr <- flank(end_gr, 1, TRUE)
        
    } else if (end_type == "E") {

        end_gr <- flank(range(x[[1]]), -1, FALSE)
      
    }

    event <- range(c(start_gr, end_gr))
        
    return(event)
    
}

getAASeq <- function(tx, genome)
{

    if (is.na(cdsStart(tx)) || is.na(cdsEnd(tx))) {

        return(NA_character_)

    }
    
    strand <- as.character(strand(tx[[1]][1]))
    cds <- sort(cds(tx), decreasing = (strand == "-"))
    cds_seq <- do.call(c, getSeq(genome, cds))
    aa_seq <- as.character(translate(cds_seq))

    return(aa_seq)
    
}

getVariantRegions <- function(ref, var, sgv, genome)
{

    ref_aa <- getAASeq(ref, genome)

    ## note the following always holds for alt == "CDS_lost"
    
    if (is.na(cdsLeft(var)) || is.na(cdsRight(var))) {
      
        data.frame(
            ref_aa_start = 0,
            ref_aa_end = nchar(ref_aa),
            var_aa_start = NA_integer_,
            var_aa_end = NA_integer_,
            ref_aa = ref_aa,
            var_aa = NA_character_,
            stringsAsFactors = FALSE)

    }

    var_aa <- getAASeq(var, genome)

    ref_loc <- range(ref[[1]])
    ref_cds <- range(cds(ref))
    
    var_loc <- range(var[[1]])
    var_cds <- range(cds(var))

    ref_cds_5p <- granges(flank(ref_cds, 1, TRUE))
    ref_cds_3p <- granges(flank(ref_cds, 1, FALSE))
    var_cds_5p <- granges(flank(var_cds, 1, TRUE))
    var_cds_3p <- granges(flank(var_cds, 1, FALSE))

    event <- reduce(c(
        getEventLocation(ref, sgv),
        getEventLocation(var, sgv)))
    
    if (!identical(ref_cds_5p, var_cds_5p)) {

        event <- range(c(event, ref_cds_5p, var_cds_5p))
      
    } else if (!identical(ref_cds_3p, var_cds_3p)) {
      
        event <- range(c(event, ref_cds_3p, var_cds_3p))
      
    }
    
    ref_event_tx <- mapEventToProtein(ref, event)
    var_event_tx <- mapEventToProtein(var, event)

    data.frame(
        ref_aa_start = ref_event_tx$start,
        ref_aa_end = ref_event_tx$end,
        var_aa_start = var_event_tx$start,
        var_aa_end = var_event_tx$end,
        ref_aa = ref_aa,
        var_aa = var_aa,
        stringsAsFactors = FALSE)
    
}

mapEventToProtein <- function(tx, event)
{

    cds <- restrict(tx, cdsLeft(tx), cdsRight(tx))
    cds <- setNames(cds, "1")
    event <- intersect(event + 1, cds[[1]])
    event <- ranges(reduce(mapToTranscripts(event, cds)))
    end(event) <- end(event) - 1
    
    if (length(event) == 0) {

        out <- data.frame(start = NA_integer_, end = NA_integer_)

    } else {
    
        out <- data.frame(
            start = floor(start(event) / 3),
            end = ceiling(end(event) / 3))

    }

    return(out)
    
}

convertToTranscripts <- function(txdb)
{

    tx <- exonsBy(txdb, "tx", use.names = TRUE)
    cds <- unlist(range(cdsBy(txdb, "tx", use.names = TRUE)))
    tx <- tx[match(names(cds), names(tx))]
    cdsLeft(tx) <- start(cds)
    cdsRight(tx) <- end(cds)

    return(tx)
    
}

checkTranscripts <- function(tx)
{

    if (!exonsOnSameChromAndStrand(tx)) {

        msg <- "All ranges in the same element of tx\n
            must be on the same chromosome and strand"
        stop(msg, call. = FALSE)

    }

    if (is.null(mcols(tx)$cdsStart) || is.null(mcols(tx)$cdsEnd)) {

        msg <- "tx must have metadata columns cdsStart and cdsEnd"
        stop(msg, call. = FALSE)

    }

    if (any(mcols(tx)$cdsStart > mcols(tx)$cdsEnd, na.rm = TRUE)) {

        msg <- "All transcripts must have cdsStart < cdsEnd"
        stop(msg, call. = FALSE)

    }

}

cdsLeft <- function(tx)
{

    mcols(tx)$cdsStart

}

'cdsLeft<-' <- function(tx, value)
{

    mcols(tx)$cdsStart <- value
    
    return(tx)
    
}

cdsRight <- function(tx)
{

    mcols(tx)$cdsEnd

}

'cdsRight<-' <- function(tx, value)
{

    mcols(tx)$cdsEnd <- value
    
    return(tx)
    
}

cdsStart <- function(tx)
{

    st <- as.character(strand(range(tx[[1]])))
    ifelse(st == "+", cdsLeft(tx), cdsRight(tx))
  
}

'cdsStart<-' <- function(tx, value)
{

    st <- as.character(strand(range(tx[[1]])))

    if (st == "+") {

        cdsLeft(tx) <- value

    } else {

        cdsRight(tx) <- value

    }

    return(tx)
  
}

cdsEnd <- function(tx) 
{

    st <- as.character(strand(range(tx[[1]])))
    ifelse(st == "+", cdsRight(tx), cdsLeft(tx))
  
}

'cdsEnd<-' <- function(tx, value)
{

    st <- as.character(strand(range(tx[[1]])))

    if (st == "+") {

        cdsRight(tx) <- value

    } else {

        cdsLeft(tx) <- value

    }

    return(tx)
  
}

alt <- function(tx)
{

    mcols(tx)$alt

}

'alt<-' <- function(tx, value)
{

    mcols(tx)$alt <- value

    return(tx)

}

cds <- function(tx)
{

    if (length(tx) > 1) {

        stop("tx must have length 1")

    }

    tx <- restrict(tx, cdsLeft(tx), cdsRight(tx))
    tx <- granges(tx[[1]])
    
    return(tx)
    
}
