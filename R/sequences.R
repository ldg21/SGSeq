##' The effect of each splice variant is assessed with respect to individual
##' protein-coding transcripts.
##'
##' @title Predict the effect of splice variants on protein-coding transcripts
##' @param sgv \code{SGVariants} object
##' @param tx A \code{TxDb} object or \code{GRangesList} of exons
##'   grouped by transcript with metadata columns \code{cdsStart} and
##'   \code{cdsEnd} (by convention, cdsStart < cdsEnd for both strands).
##'   For import from GFF format, use function \code{importTranscripts}.
##' @param genome \code{BSgenome} object
##' @param summarize Logical indicating whether results should be
##'   summarized per variant
##' @param cores Number of cores available for parallel processing
##' @return For \code{summarize = FALSE} a \code{data.frame} with rows
##'   corresponding to a variant-transcript pair. The \code{data.frame}
##'   includes columns for variant identifier, transcript name, type of
##'   alteration, protein sequences for the reference transcript and the
##'   transcript variant, protein lengths and coordinates of the variant
##'   in the protein sequences. Start and end coordinates are 0- and
##'   1-based, respectively, to allow for specification of deletions.
##'   For \code{summarize = TRUE} a character vector matching argument
##'   \code{sgv} with comma-separated predicted alterations for
##'   individual transcripts.
##' @examples
##' require(BSgenome.Hsapiens.UCSC.hg19)
##' seqlevelsStyle(Hsapiens) <- "NCBI"
##' predictVariantEffects(sgv_pred, tx, Hsapiens)
##' @author Leonard Goldstein

predictVariantEffects <- function(sgv, tx, genome, summarize = TRUE, cores = 1)
{

    if (is(tx, "TxDb")) {

        message("Obtaining transcripts from TxDb...")
        tx <- convertToTranscripts(tx)

    } else if (is(tx, "GRangesList")) {

        checkTranscriptFormat(tx)

        if (is.null(names(tx))) {

            names(tx) <- seq_along(tx)

        }

    } else {

        stop("tx must be a TxDb object or GRangesList of exons
            grouped by transcripts")

    }

    sgv_vid <- variantID(sgv)

    excl <- grep("(", type(sgv), fixed = TRUE)

    if (length(excl) > 0) {

        sgv <- sgv[-excl]
        msg <- paste("Excluded", length(excl), "nested variants")
        warning(msg, call. = FALSE, immediate. = TRUE)

    }

    sgv_tx <- mapVariantsToTranscripts(sgv, tx, cores)
    tx <- tx[names(tx) %in% unlist(sgv_tx)]
    tx <- processTranscripts(tx, genome, cores)
    sgv_tx <- lapply(sgv_tx, intersect, names(tx))

    msg <- paste("Predicting effect of", length(sgv),
        "variants on", length(tx), "coding transcripts...")
    message(msg)

    expanded_sgv <- sgv[rep(seq_along(sgv), elementLengths(sgv_tx))]
    expanded_tx <- tx[match(unlist(sgv_tx), names(tx))]

    list_expanded_sgv <- split(expanded_sgv, seq_along(expanded_sgv))
    list_expanded_tx <- split(expanded_tx, seq_along(expanded_tx))

    list_res <- mcmapply(
        predictVariantEffectPerVariantAndTranscript,
        list_expanded_sgv,
        list_expanded_tx,
        MoreArgs = list(genome = genome),
        SIMPLIFY = FALSE,
        USE.NAMES = FALSE,
        mc.preschedule = setPreschedule(cores),
        mc.cores = cores)

    items <- paste("variant", variantID(expanded_sgv),
        "and transcript", names(expanded_tx))

    checkApplyResultsForErrors(
        list_res,
        "predictVariantEffectPerVariantAndTranscript\n",
        items,
        "character")

    res <- do.call(rbindDfsWithoutRowNames, list_res)

    if (summarize) {

        res$effect <- paste0(res$txName, ":", res$effect)
        vid_effect <- unstrsplit(split(res$effect, res$variantID), ",")
        res <- vid_effect[match(sgv_vid, names(vid_effect))]
        names(res) <- NULL

    }

    return(res)

}

mapVariantsToTranscripts <- function(sgv, tx, cores)
{

    gr_from <- pos2gr(sub("^(D|S):", "", from(sgv)))
    gr_to <- pos2gr(sub("^(A|E):", "", to(sgv)))

    list_D <- vector("list", length(sgv))
    i_D <- which(substr(from(sgv), 1, 1) == "D")
    list_D[i_D] <- as.list(findOverlaps(gr_from[i_D], tx))

    list_A <- vector("list", length(sgv))
    i_A <- which(substr(to(sgv), 1, 1) == "A")
    list_A[i_A] <- as.list(findOverlaps(gr_to[i_A], tx))

    list_i <- vector("list", length(sgv))
    i <- intersect(i_D, i_A)
    list_i[i] <- mcmapply(
        intersect,
        list_D[i],
        list_A[i],
        SIMPLIFY = FALSE,
        mc.cores = cores)
    i <- setdiff(i_D, i_A)
    list_i[i] <- list_D[i]
    i <- setdiff(i_A, i_D)
    list_i[i] <- list_A[i]

    sgv_tx <- relist(names(tx)[unlist(list_i)], list_i)

    inc <- tx[match(unique(unlist(sgv_tx)), names(tx))]
    sgv <- annotate(sgv, convertToTxFeatures(inc))

    sgv_tx <- mcmapply(
        setdiff,
        sgv_tx,
        txName(sgv),
        SIMPLIFY = FALSE,
        mc.cores = cores)

    return(sgv_tx)

}

predictVariantEffectPerVariantAndTranscript <- function(sgv, ref, genome)
{

    var <- getTranscriptVariant(sgv, ref, genome)

    n <- length(var)

    res <- data.frame(
        variantID = rep(variantID(sgv), n),
        txName = rep(names(ref), n),
        effect = effect(var),
        stringsAsFactors = FALSE)

    list_var <- split(var, seq_along(var))

    regions <- do.call(rbindDfsWithoutRowNames, lapply(
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

        if (identical(gr2co(cds(var)), gr2co(cds(ref)))) {

            effect(var) <- "5p_UTR_variant"

        } else {

            effect(var) <- "upstream_start_codon"

        }

    } else if (event_start %over% utr_5p && event_end %over% ref_cds) {

        var <- findCDS(var, NA_integer_, cdsEnd(ref), genome)

        if (identical(gr2co(cds(var)), gr2co(cds(ref)))) {

            effect(var) <- "protein_unaffected"

        } else if (is.na(cdsStart(var))) {

            effect(var) <- "start_codon_lost"

        } else {

            effect(var) <- "N-terminal_variant"

        }

    } else if (event_start %over% utr_5p && event_end %over% utr_3p) {

        cdsStart(var) <- NA_integer_
        cdsEnd(var) <- NA_integer_
        effect(var) <- "CDS_lost"

    } else if (event_start %over% ref_cds && event_end %over% ref_cds) {

        var_3p <- findCDS(var, cdsStart(ref), NA_integer_, genome)
        var_5p <- findCDS(var, NA_integer_, cdsEnd(ref), genome)

        del_w <- sum(width(del))
        ins_w <- sum(width(ins))

        if (del_w %% 3 == ins_w %% 3) {

            effect_cds <- "in-frame"

        } else {

            effect_cds <- "frame-shift"

        }

        if (del_w > 0 && ins_w > 0) {

            effect_tx <- "alteration"

        } else if (del_w > 0) {

            effect_tx <- "deletion"

        } else if (ins_w > 0) {

            effect_tx <- "insertion"

        }

        prefix <- paste0(effect_cds, "_", effect_tx)

        if (identical(gr2co(cds(var_3p)), gr2co(cds(ref)))) {

            stop("internal CDS variant does not affect CDS")

        }

        if (effect_cds == "in-frame") {

            if (cdsEnd(var_3p) == cdsEnd(ref)) {

                effect(var_3p) <- prefix
                var <- var_3p

            } else {

                effect(var_3p) <- paste0(prefix, "_premature_stop")
                var <- var_3p

                if (!is.na(cdsStart(var_5p))) {

                    effect(var_5p) <- paste0(prefix, "_N-terminal_variant")
                    var <- c(var, var_5p)

                }

            }

        } else if (effect_cds == "frame-shift") {

            if (is.na(cdsEnd(var_3p)) && is.na(cdsStart(var_5p))) {

                effect(var_3p) <- paste0(prefix, "_CDS_lost")
                var <- var_3p

            } else if (!is.na(cdsEnd(var_3p)) && is.na(cdsStart(var_5p))) {

                effect(var_3p) <- paste0(prefix, "_C-terminal_variant")
                var <- var_3p

            } else if (is.na(cdsEnd(var_3p)) && !is.na(cdsStart(var_5p))) {

                effect(var_5p) <- paste0(prefix, "_N-terminal_variant")
                var <- var_5p

            } else if (!is.na(cdsEnd(var_3p)) && !is.na(cdsStart(var_5p))) {

                effect(var_3p) <- paste0(prefix, "_C-terminal_variant")
                effect(var_5p) <- paste0(prefix, "_N-terminal_variant")
                var <- c(var_3p, var_5p)

            }

        }

    } else if (event_start %over% ref_cds && event_end %over% utr_3p) {

        var <- findCDS(var, cdsStart(ref), NA_integer_, genome)

        if (identical(gr2co(cds(var)), gr2co(cds(ref)))) {

            effect(var) <- "protein_unaffected"

        } else if (is.na(cdsEnd(var))) {

            effect(var) <- "stop_codon_lost"

        } else {

            effect(var) <- "C-terminal_variant"

        }

    } else if (event_start %over% utr_3p && event_end %over% utr_3p) {

        cdsStart(var) <- cdsStart(ref)
        cdsEnd(var) <- cdsEnd(ref)
        effect(var) <- "3p_UTR_variant"

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
    tx_seq <- as.character(do.call(c,
        suppressWarnings(getSeq(genome, tx[[1]]))))

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
    cds_seq <- do.call(c, suppressWarnings(getSeq(genome, cds)))
    aa <- as.character(translate(cds_seq))
    aa <- sub("\\*$", "", aa)

    return(aa)

}

getVariantRegions <- function(ref, var, sgv, genome)
{

    ref_aa <- getAASeq(ref, genome)

    ## note the following always holds for effect == "CDS_lost"

    if (is.na(cdsLeft(var)) || is.na(cdsRight(var))) {

        data.frame(
            ref_aa_start = 0,
            ref_aa_end = nchar(ref_aa),
            ref_aa_length = nchar(ref_aa),
            var_aa_start = NA_integer_,
            var_aa_end = NA_integer_,
            var_aa_length = NA_integer_,
            ref_aa = ref_aa,
            var_aa = NA_character_,
            stringsAsFactors = FALSE)

    }

    var_aa <- getAASeq(var, genome)

    ref_loc <- range(ref[[1]])
    ref_cds <- range(cds(ref))

    var_loc <- range(var[[1]])
    var_cds <- range(cds(var))

    ref_cds_5p_fl <- granges(flank(ref_cds, 1, TRUE))
    ref_cds_3p_fl <- granges(flank(ref_cds, 1, FALSE))
    var_cds_5p_fl <- granges(flank(var_cds, 1, TRUE))
    var_cds_3p_fl <- granges(flank(var_cds, 1, FALSE))

    event <- reduce(c(
        getEventLocation(ref, sgv),
        getEventLocation(var, sgv)))

    if (start(ref_cds_5p_fl) != start(var_cds_5p_fl)) {

        event <- range(c(event, ref_cds_5p_fl, var_cds_5p_fl))

    }

    if (start(ref_cds_3p_fl) != start(var_cds_3p_fl)) {

        event <- range(c(event, ref_cds_3p_fl, var_cds_3p_fl))

    }

    ref_event_tx <- mapEventToProtein(ref, event)
    var_event_tx <- mapEventToProtein(var, event)

    data.frame(
        ref_start = ref_event_tx$start,
        ref_end = min(ref_event_tx$end, nchar(ref_aa)),
        ref_length = nchar(ref_aa),
        var_start = var_event_tx$start,
        var_end = min(var_event_tx$end, nchar(var_aa)),
        var_length = nchar(var_aa),
        ref_aa = ref_aa,
        var_aa = var_aa,
        stringsAsFactors = FALSE)

}

mapEventToProtein <- function(tx, event)
{

    cds <- setNames(restrict(tx, cdsLeft(tx), cdsRight(tx)), "1")
    event <- intersect(event + 1, cds[[1]])
    event <- ranges(reduce(mapToTranscripts(event, cds)))
    event <- flank(event, -width(event) + 1, TRUE)

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
    cdsLeft(tx) <- start(cds)[match(names(tx), names(cds))]
    cdsRight(tx) <- end(cds)[match(names(tx), names(cds))]

    return(tx)

}

checkTranscriptFormat <- function(tx)
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

        msg <- "All coding transcripts must have cdsStart < cdsEnd"
        stop(msg, call. = FALSE)

    }

}

processTranscripts <- function(tx, genome, cores)
{

    tx <- tx[!is.na(cdsLeft(tx)) & !is.na(cdsRight(tx))]

    list_tx <- split(tx, seq_along(tx))

    tx_aa <- unlist(mclapply(list_tx, getAASeq, genome, mc.cores = cores))

    tx_aa_first <- substr(tx_aa, 1, 1)

    list_pr <- mcmapply(
        findCDS,
        list_tx,
        cdsStart(tx),
        rep(NA_integer_, length(tx)),
        MoreArgs = list(genome = genome),
        SIMPLIFY = FALSE,
        mc.cores = cores)

    pr <- setNames(Reduce(c, list_pr), names(tx))

    pr_aa <- unlist(mclapply(list_pr, getAASeq, genome, mc.cores = cores))

    invalid <- which(is.na(tx_aa) | is.na(pr_aa) |
        tx_aa_first != "M" | tx_aa != pr_aa)

    if (length(invalid) > 0) {

        pr <- pr[-invalid]
        msg <- paste0("Excluded ", length(invalid),
            " transcripts with invalid CDS\n",
            paste(names(tx)[invalid], collapse = "\n"))
        warning(msg, call. = FALSE, immediate. = TRUE)

    }

    return(pr)

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

    st <- as.character(strand(unlist(range(tx))))
    ifelse(st == "+", cdsLeft(tx), cdsRight(tx))

}

'cdsStart<-' <- function(tx, value)
{

    st <- as.character(strand(unlist(range(tx))))
    i_pos <- which(st == "+")
    cdsLeft(tx)[i_pos] <- value[i_pos]
    i_neg <- which(st == "-")
    cdsRight(tx)[i_neg] <- value[i_neg]

    return(tx)

}

cdsEnd <- function(tx)
{

    st <- as.character(strand(unlist(range(tx))))
    ifelse(st == "+", cdsRight(tx), cdsLeft(tx))

}

'cdsEnd<-' <- function(tx, value)
{

    st <- as.character(strand(unlist(range(tx))))
    i_pos <- which(st == "+")
    cdsRight(tx)[i_pos] <- value[i_pos]
    i_neg <- which(st == "-")
    cdsLeft(tx)[i_neg] <- value[i_neg]

    return(tx)

}

effect <- function(tx)
{

    mcols(tx)$effect

}

'effect<-' <- function(tx, value)
{

    mcols(tx)$effect <- value

    return(tx)

}

cds <- function(tx)
{

    if (length(tx) > 1) {

        stop("tx must have length 1")

    }

    tx <- restrict(tx, cdsLeft(tx), cdsRight(tx))
    tx <- sort(granges(tx[[1]]))

    return(tx)

}
