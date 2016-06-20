##' The effect of a splice variant is predicted for individual
##' protein-coding transcripts.
##'
##' @title Predict the effect of splice variants on protein-coding transcripts
##' @param sgv \code{SGVariants} object
##' @param tx \code{TxDb} object, or \code{GRangesList} of exons
##'   grouped by transcript with metadata columns \code{txName},
##'   \code{geneName}, \code{cdsStart} and \code{cdsEnd}
##'   (by convention, cdsStart < cdsEnd for both strands).
##'   For import from GFF format, use function \code{importTranscripts}.
##' @param genome \code{BSgenome} object
##' @param fix_start_codon Logical indicating whether the annotated start
##'   codon should be considered fixed and the variant transcript should
##'   not be scanned for alternative start codons
##' @param output Character string indicating whether short results or
##'   full results (with additional columns) should be returned
##' @param cores Number of cores available for parallel processing
##' @return \code{data.frame} with rows corresponding to a
##'   variant-transcript pair. The output includes columns for variant
##'   identifier, transcript name, gene name, type of alteration at the
##'   RNA and protein level, and variant description at the RNA and
##'   protein level in HGVS notation. For \code{output = "full"}
##'   additional columns are returned. These include the full-length RNA
##'   and protein sequence for the reference and variant transcript.
##'   Event start and end coordinates in the full output are 0- and
##'   1-based, respectively (to allow for description of deletions).
##'   Coordinates for the last junction in a transcript refer to the
##'   last base of the second-to-last exon.
##' @examples
##' require(BSgenome.Hsapiens.UCSC.hg19)
##' seqlevelsStyle(Hsapiens) <- "NCBI"
##' predictVariantEffects(sgv_pred, tx, Hsapiens)
##' @author Leonard Goldstein

predictVariantEffects <- function(sgv, tx, genome, fix_start_codon = TRUE,
    output = c("short", "full"), cores = 1)
{

    if (missing(genome)) {

        stop("missing argument genome")

    }

    output <- match.arg(output)

    sgv <- updateObject(sgv, verbose = TRUE)

    if (is(tx, "TxDb")) {

        tx <- convertToTranscripts(tx)

    } else if (is(tx, "GRangesList")) {

        checkTranscriptFormat(tx)

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

    expanded_sgv <- sgv[togroup0(sgv_tx)]
    expanded_tx <- tx[match(unlist(sgv_tx), names(tx))]

    list_expanded_sgv <- split(expanded_sgv, seq_along(expanded_sgv))
    list_expanded_tx <- split(expanded_tx, seq_along(expanded_tx))

    list_res <- mcmapply(
        predictVariantEffectPerVariantAndTranscript,
        list_expanded_sgv,
        list_expanded_tx,
        MoreArgs = list(
            genome = genome,
            fix_start_codon = fix_start_codon,
            output = output),
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

    inc <- tx[unique(unlist(list_i))]
    sgv <- annotate(sgv, convertToTxFeatures(inc))

    sgv_tx <- relist(names(tx)[unlist(list_i)], list_i)
    sgv_tx <- mcmapply(
        setdiff,
        sgv_tx,
        txName(sgv),
        SIMPLIFY = FALSE,
        mc.cores = cores)

    return(sgv_tx)

}

predictVariantEffectPerVariantAndTranscript <- function(sgv, ref, genome,
    fix_start_codon, output)
{

    ref_loc <- range(ref[[1]])
    ref_cds <- range(cds(ref))

    ref_utr_5p <- range(c(flank(ref_loc, -1, TRUE), flank(ref_cds, 1, TRUE)))
    ref_utr_3p <- range(c(flank(ref_cds, 1, FALSE), flank(ref_loc, -1, FALSE)))

    event <- getEventLocation(sgv, ref)
    event_start <- flank(event, -1, TRUE)
    event_end <- flank(event, -1, FALSE)

    del <- intersect(ref[[1]], event)
    ins <- reduce(granges(sgv[[1]][type(sgv[[1]]) == "E"]))

    var <- GRangesList(reduce(c(setdiff(ref[[1]], del), ins)))
    var_3p <- getVariant(var, ref, sgv, genome, fix = "start")
    var_5p <- getVariant(var, ref, sgv, genome, fix = "stop")

    alt_RNA <- getHGVSVariantRNA(var_3p, ref)

    if (event_start %over% ref_utr_5p && event_end %over% ref_utr_5p) {

        if (fix_start_codon) {

            var <- var_3p
            alt_protein <- getHGVSVariantProteinUnaffected()

        } else {

            var <- var_5p

            if (cdsStart(var) == cdsStart(ref)) {

                alt_protein <- getHGVSVariantProteinUnaffected()

            } else {

                alt_protein <- getHGVSVariantNTerminalExtension(var)

            }

        }

    } else if (event_start %over% ref_utr_5p && event_end %over% ref_cds) {

        var <- var_5p

        if (identical(gr2co(cds(var)), gr2co(cds(ref)))) {

            alt_protein <- getHGVSVariantProteinUnaffected()

        } else if (fix_start_codon || is.na(cdsStart(var))) {

            alt_protein <- getHGVSVariantNoProtein()

        } else {

            alt_protein <- getHGVSVariantNTerminalVariant(var)

        }

    } else if (event_start %over% ref_utr_5p && event_end %over% ref_utr_3p) {

        var <- var_3p
        cdsStart(var) <- NA_integer_
        alt_protein <- getHGVSVariantNoProtein()

    } else if (event_start %over% ref_cds && event_end %over% ref_cds) {

        if (sum(width(del)) %% 3 == sum(width(ins)) %% 3) {

            var <- var_3p
            alt_protein <- getHGVSVariantDeletionInsertion(var_3p)

            if (cdsEnd(var_3p) != cdsEnd(ref)
                && !fix_start_codon && !is.na(cdsStart(var_5p))) {

                var <- c(var, var_5p)
                alt_protein <- pc(alt_protein,
                    getHGVSVariantNTerminalVariant(var_5p))

            }

        } else {

            var <- var_3p
            alt_protein <- getHGVSVariantFrameshift(var_3p)

            if (!fix_start_codon && !is.na(cdsStart(var_5p))) {

                var <- c(var, var_5p)
                alt_protein <- pc(alt_protein,
                    getHGVSVariantNTerminalVariant(var_5p))

            }

        }

    } else if (event_start %over% ref_cds && event_end %over% ref_utr_3p) {

        var <- var_3p

        if (identical(gr2co(cds(var)), gr2co(cds(ref)))) {

            alt_protein <- getHGVSVariantProteinUnaffected()

        } else {

            alt_protein <- getHGVSVariantFrameshift(var)

        }

    } else if (event_start %over% ref_utr_3p && event_end %over% ref_utr_3p) {

        var <- var_3p
        alt_protein <- getHGVSVariantProteinUnaffected()

    }

    n <- length(var)

    res <- data.frame(
        variantID = rep(variantID(sgv), n),
        txName = rep(mcols(ref)$txName, n),
        geneName = rep(mcols(ref)$geneName, n),
        RNA_change = rep(alt_RNA$HGVS, n),
        RNA_variant_type = rep(alt_RNA$type, n),
        protein_change = alt_protein$HGVS,
        protein_variant_type = alt_protein$type,
        stringsAsFactors = FALSE)

    res <- cbind(res, as.data.frame(mcols(var)))

    if (output == "short") {

        res <- res[c(
            "variantID",
            "txName",
            "geneName",
            "RNA_change",
            "RNA_variant_type",
            "protein_change",
            "protein_variant_type")]

    } else if (output == "full") {

        res <- res[c(
            "variantID",
            "txName",
            "geneName",
            "RNA_change",
            "RNA_variant_type",
            "RNA_ref_length",
            "RNA_ref_cds_start",
            "RNA_ref_cds_end",
            "RNA_ref_last_junction",
            "RNA_ref_event_start",
            "RNA_ref_event_end",
            "RNA_var_length",
            "RNA_var_cds_start",
            "RNA_var_cds_end",
            "RNA_var_last_junction",
            "RNA_var_event_start",
            "RNA_var_event_end",
            "protein_change",
            "protein_variant_type",
            "protein_ref_length",
            "protein_ref_event_start",
            "protein_ref_event_end",
            "protein_var_length",
            "protein_var_event_start",
            "protein_var_event_end",
            "RNA_ref_seq",
            "RNA_var_seq",
            "protein_ref_seq",
            "protein_var_seq")]

    }

    return(res)

}

getHGVSVariantRNA <- function(var, ref)
{

    ref_start <- mcols(var)$RNA_ref_event_start
    ref_end <- mcols(var)$RNA_ref_event_end
    var_start <- mcols(var)$RNA_var_event_start
    var_end <- mcols(var)$RNA_var_event_end

    type <- c()

    if (ref_end > ref_start) {

        hgvs <- getHGVSDeletion(var, "RNA")
        type <- c(type, "deletion")

    }

    if (var_end > var_start) {

        if (ref_end == ref_start) {

            hgvs <- getHGVSFlanking(var, "RNA")

        }

        hgvs <- paste0(hgvs, getHGVSInsertion(var, "RNA", ref))
        type <- c(type, "insertion")

    }

    list(type = paste(type, collapse = "/"), HGVS = hgvs)

}

getHGVSVariantNTerminalExtension <- function(var)
{

    ext <- mcols(var)$protein_var_length - mcols(var)$protein_ref_length
    hgvs <- paste0("p.M1ext-", ext)
    list(type = "N-terminal_extension", HGVS = hgvs)

}

getHGVSVariantNTerminalVariant <- function(var) {

    hgvs <- getHGVSDeletion(var, "protein")
    type <- "N-terminal_deletion"

    var_l <- mcols(var)$protein_var_length
    ref_l <- mcols(var)$protein_ref_length
    del_l <- mcols(var)$protein_ref_event_end -
        mcols(var)$protein_ref_event_start
    ext <- var_l - (ref_l - del_l)

    if (ext > 0) {

        hgvs <- paste0(hgvs, "ext-", ext)
        type <- "N-terminal_variant"

    }

    list(type = type, HGVS = hgvs)

}

getHGVSVariantDeletionInsertion <- function(var)
{

    ref_start <- mcols(var)$protein_ref_event_start
    ref_end <- mcols(var)$protein_ref_event_end
    var_start <- mcols(var)$protein_var_event_start
    var_end <- mcols(var)$protein_var_event_end

    type <- c()

    if (ref_end > ref_start) {

        hgvs <- getHGVSDeletion(var, "protein")
        type <- c(type, "deletion")

    }

    if (var_end > var_start) {

        if (ref_end == ref_start) {

            hgvs <- getHGVSFlanking(var, "protein")

        }

        hgvs <- paste0(hgvs, getHGVSInsertion(var, "protein"))
        type <- c(type, "insertion")

    }

    out <- list(
        type = paste0("in-frame_", paste(type, collapse = "/")),
        HGVS = hgvs)

    return(out)

}

getHGVSVariantFrameshift <- function(var)
{

    ref_start <- mcols(var)$protein_ref_event_start
    var_start <- mcols(var)$protein_var_event_start

    res <- getHGVSRefPos(var, ref_start + 1, "protein")
    alt <- substr(mcols(var)$protein_var_seq, var_start + 1, var_start + 1)
    hgvs <- paste0("p.", res, alt, "fs*")

    if (!is.na(cdsEnd(var))) {

        ## need to add 1 since sequence does not include stop codon
        ext <- mcols(var)$protein_var_length - var_start + 1
        hgvs <- paste0(hgvs, ext)

    } else {

        hgvs <- paste0(hgvs, "?")

    }

    list(type = "frame_shift", HGVS = hgvs)

}

getHGVSVariantProteinUnaffected <- function()
{

    list(type = "no_change", HGVS = "p.=")

}

getHGVSVariantNoProtein <- function()
{

    list(type = "no_protein", HGVS = "p.0")

}

getHGVSDeletion <- function(var, type)
{

    suffix <- switch(type, RNA = "r", protein = "p")
    ref_start <- mcols(var)[[paste0(type, "_ref_event_start")]]
    ref_end <- mcols(var)[[paste0(type, "_ref_event_end")]]
    del_start <- getHGVSRefPos(var, ref_start + 1, type)
    del_end <- getHGVSRefPos(var, ref_end, type)
    paste0(suffix, ".", del_start, "_", del_end, "del")

}

getHGVSFlanking <- function(var, type)
{

    suffix <- switch(type, RNA = "r", protein = "p")
    ref_start <- mcols(var)[[paste0(type, "_ref_event_start")]]
    ref_end <- mcols(var)[[paste0(type, "_ref_event_end")]]
    flanking_start <- getHGVSRefPos(var, ref_start, type)
    flanking_end <- getHGVSRefPos(var, ref_end + 1, type)
    paste0(suffix, ".", flanking_start, "_", flanking_end)

}

getHGVSInsertion <- function(var, type, ref)
{

    if (type == "RNA") {

        strand <- as.character(strand(ref[[1]][1]))

        I <- setdiff(var[[1]], ref[[1]])
        I <- sort(I, decreasing = (strand == "-"))

        ins <- vector()

        for(i in seq_along(I)) {

            I5p <- start(flank(I[i], -1, TRUE))
            I3p <- start(flank(I[i], -1, FALSE))

            ins <- c(ins, paste0(
                getHGVSRefPosIntronic(I5p, ref, var), "_",
                getHGVSRefPosIntronic(I3p, ref, var)))

        }

        if (length(ins) > 1) {

            ins <- paste0("[", paste(ins, collapse = ";"), "]")

        }

    } else if (type == "protein") {

        ins <- substr(mcols(var)$protein_var_seq,
            mcols(var)$protein_var_event_start + 1,
            mcols(var)$protein_var_event_end)

    }

    ins <- paste0("ins", ins)

    return(ins)

}

getHGVSRefPos <- function(var, pos, type)
{

    if (type == "RNA") {

        cdsStart <- mcols(var)$RNA_ref_cds_start
        cdsEnd <- mcols(var)$RNA_ref_cds_end

        if (pos < cdsStart) {

            pos <- pos - cdsStart

        } else if (pos >= cdsStart && pos <= cdsEnd) {

            pos <- pos - cdsStart + 1

        } else {

            pos <- paste0("*", pos - cdsEnd)

        }

        pos <- as.character(pos)

    } else if (type == "protein") {

        pos <- paste0(substr(mcols(var)$protein_ref_seq, pos, pos), pos)

    }

    return(pos)

}

getHGVSRefPosIntronic <- function(pos, ref, var)
{

    names(ref) <- "1"
    strand <- as.character(strand(ref[[1]][1]))

    E <- ref[[1]]
    E <- sort(E, decreasing = (strand == "-"))

    E5p <- flank(E, -1, TRUE)
    E3p <- flank(E, -1, FALSE)

    R5p <- c(start(mapToTranscripts(E5p, ref)), NA)
    R3p <- c(NA, start(mapToTranscripts(E3p, ref)))

    E5p <- c(start(E5p), NA)
    E3p <- c(NA, start(E3p))

    d5p <- (pos - E3p) * ifelse(strand == "-", -1, 1)
    d3p <- (E5p - pos) * ifelse(strand == "-", -1, 1)

    i <- which((is.na(d5p) | d5p > 0) & (is.na(d3p) | d3p > 0))

    if (is.na(d3p[i]) ||
        (!is.na(d5p[i]) && !is.na(d3p[i]) && d5p[i] <= d3p[i])) {

        pos <- paste0(getHGVSRefPos(var, R3p[i], "RNA"), "+", d5p[i])

    } else {

        pos <- paste0(getHGVSRefPos(var, R5p[i], "RNA"), "-", d3p[i])

    }

    return(pos)

}

getVariant <- function(var, ref, sgv, genome, fix = c("start", "stop"))
{

    fix <- match.arg(fix)

    if (fix == "start") {

        var <- findCDS(var, cdsStart(ref), NA_integer_, genome)

    } else {

        var <- findCDS(var, NA_integer_, cdsEnd(ref), genome)

    }

    var <- getVariantRNA(var, ref, sgv, genome)
    var <- getVariantProtein(var, ref, sgv, genome)

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

    if (is.na(cdsEnd)) {

        cdsStart_gr <- GRanges(chrom, IRanges(cdsStart, cdsStart), strand)

        if (!cdsStart_gr %over% tx) {

            cdsStart <- NA_integer_
            cdsEnd <- NA_integer_

        } else {

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

    } else {

        cdsEnd_gr <- GRanges(chrom, IRanges(cdsEnd, cdsEnd), strand)

        if (!cdsEnd_gr %over% tx) {

            cdsStart <- NA_integer_
            cdsEnd <- NA_integer_

        } else {

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

        }

    }

    cdsStart(tx) <- cdsStart
    cdsEnd(tx) <- cdsEnd

    return(tx)

}

getEventLocation <- function(sgv, ref)
{

    start_type <- substr(from(sgv), 1, 1)
    end_type <- substr(to(sgv), 1, 1)

    if (start_type == "D") {

        start_gr <- pos2gr(sub("^D:", "", from(sgv)))
        start_gr <- flank(start_gr, 1, FALSE)

    } else if (start_type == "S") {

        start_gr <- flank(range(ref[[1]]), -1, TRUE)

    }

    if (end_type == "A") {

        end_gr <- pos2gr(sub("^A:", "", to(sgv)))
        end_gr <- flank(end_gr, 1, TRUE)

    } else if (end_type == "E") {

        end_gr <- flank(range(ref[[1]]), -1, FALSE)

    }

    event <- range(c(start_gr, end_gr))

    return(event)

}

getRNASeq <- function(tx, genome)
{

    tx <- tx[[1]]
    strand <- as.character(strand(tx[1]))
    tx <- sort(tx, decreasing = (strand == "-"))
    nt <- do.call(c, suppressWarnings(getSeq(genome, tx)))
    nt <- as.character(nt)
    nt <- gsub("A", "a", nt)
    nt <- gsub("C", "c", nt)
    nt <- gsub("G", "g", nt)
    nt <- gsub("T", "u", nt)
    nt <- gsub("N", "n", nt)

    return(nt)

}

getVariantRNA <- function(var, ref, sgv, genome)
{

    strand <- as.character(strand(ref[[1]][1]))

    names(ref) <- "1"
    names(var) <- "1"

    ref[[1]] <- sort(ref[[1]], decreasing = (strand == "-"))
    var[[1]] <- sort(var[[1]], decreasing = (strand == "-"))

    ref_seq <- getRNASeq(ref, genome)
    var_seq <- getRNASeq(var, genome)

    ref_cds <- reduce(mapToTranscripts(cds(ref), ref))
    var_cds <- reduce(mapToTranscripts(cds(var), var))

    ref_event <- getEventLocation(sgv, ref) + 1
    ref_event <- intersect(ref_event, ref[[1]])
    ref_event <- reduce(mapToTranscripts(ref_event, ref))

    var_event <- getEventLocation(sgv, var) + 1
    var_event <- intersect(var_event, var[[1]])
    var_event <- reduce(mapToTranscripts(var_event, var))

    if (length(ref_event) > 1 || length(var_event) > 1) stop()

    if (grepl("^S", from(sgv))) {

        ref_event_start <- 0L
        var_event_start <- 0L

    } else {

        ref_event_start <- start(ref_event)
        var_event_start <- start(var_event)

    }

    if (grepl("^E", to(sgv))) {

        ref_event_end <- nchar(ref_seq)
        var_event_end <- nchar(var_seq)

    } else {

        ref_event_end <- end(ref_event) - 1L
        var_event_end <- end(var_event) - 1L

    }

    ref_last_exon <- ref[[1]][length(ref[[1]])]
    var_last_exon <- var[[1]][length(var[[1]])]

    ref_last_junction <- start(mapToTranscripts(ref_last_exon, ref)) - 1L
    var_last_junction <- start(mapToTranscripts(var_last_exon, var)) - 1L

    if (length(ref[[1]]) == 1) ref_last_junction <- NA_integer_
    if (length(var[[1]]) == 1) var_last_junction <- NA_integer_

    df <- data.frame(
        RNA_ref_cds_start = start(ref_cds),
        RNA_ref_cds_end = end(ref_cds),
        RNA_ref_last_junction = ref_last_junction,
        RNA_ref_event_start = ref_event_start,
        RNA_ref_event_end = ref_event_end,
        RNA_ref_length = nchar(ref_seq),
        RNA_ref_seq = ref_seq,
        RNA_var_cds_start = start(var_cds),
        RNA_var_cds_end = end(var_cds),
        RNA_var_last_junction = var_last_junction,
        RNA_var_event_start = var_event_start,
        RNA_var_event_end = var_event_end,
        RNA_var_length = nchar(var_seq),
        RNA_var_seq = var_seq,
        stringsAsFactors = FALSE)

    mcols(var) <- cbind(mcols(var), DataFrame(df))

    return(var)

}

getProteinSeq <- function(tx, genome)
{

    if (is.na(cdsStart(tx)) || is.na(cdsEnd(tx))) {

        return(NA_character_)

    }

    cds_seq <- do.call(c, suppressWarnings(getSeq(genome, cds(tx))))
    protein <- as.character(translate(cds_seq))
    protein <- sub("\\*$", "", protein)

    return(protein)

}

getVariantProtein <- function(var, ref, sgv, genome)
{

    ref_seq <- getProteinSeq(ref, genome)

    if (is.na(cdsStart(var)) || is.na(cdsEnd(var))) {

        df <- data.frame(
            protein_ref_event_start = 0L,
            protein_ref_event_end = nchar(ref_seq),
            protein_ref_length = nchar(ref_seq),
            protein_ref_seq = ref_seq,
            protein_var_event_start = NA_integer_,
            protein_var_event_end = NA_integer_,
            protein_var_length = NA_integer_,
            protein_var_seq = NA_character_,
            stringsAsFactors = FALSE)

        mcols(var) <- cbind(mcols(var), DataFrame(df))

        return(var)

    }

    var_seq <- getProteinSeq(var, genome)

    ref_loc <- range(ref[[1]])
    ref_cds <- range(cds(ref))

    var_loc <- range(var[[1]])
    var_cds <- range(cds(var))

    ref_cds_start <- granges(flank(ref_cds, -1, TRUE))
    ref_cds_end <- granges(flank(ref_cds, -1, FALSE))
    var_cds_start <- granges(flank(var_cds, -1, TRUE))
    var_cds_end <- granges(flank(var_cds, -1, FALSE))

    event <- reduce(c(
        getEventLocation(sgv, ref),
        getEventLocation(sgv, var)))

    if (start(ref_cds_start) != start(var_cds_start)) {

        event <- range(c(event, ref_cds_start, var_cds_start))

    }

    if (start(ref_cds_end) != start(var_cds_end)) {

        event <- range(c(event, ref_cds_end, var_cds_end))

    }

    event_ref <- mapEventToProtein(ref, event)
    if (!is.na(event_ref$end) && event_ref$end == nchar(ref_seq) + 1)
        event_ref$end <- nchar(ref_seq)
    event_var <- mapEventToProtein(var, event)
    if (!is.na(event_var$end) && event_var$end == nchar(var_seq) + 1)
        event_var$end <- nchar(var_seq)

    if (!is.na(event_ref$start) && !is.na(event_ref$end) &&
        !is.na(event_var$start) && !is.na(event_var$end)) {

        if (event_ref$end > event_ref$start &&
            event_var$end > event_var$start &&
            start(ref_cds_start) == start(var_cds_start)) {

            del <- substr(ref_seq, event_ref$start + 1, event_ref$end)
            ins <- substr(var_seq, event_var$start + 1, event_var$end)

            D <- nchar(del)
            I <- nchar(ins)

            P <- min(D, I)
            p <- 0

            while (p < P &&
                substr(del, p + 1, p + 1) == substr(ins, p + 1, p + 1)) {

                p <- p + 1

            }

            event_ref$start <- event_ref$start + p
            event_var$start <- event_var$start + p

        }

        if (event_ref$end > event_ref$start &&
            event_var$end > event_var$start &&
            start(ref_cds_end) == start(var_cds_end)) {

            del <- substr(ref_seq, event_ref$start + 1, event_ref$end)
            ins <- substr(var_seq, event_var$start + 1, event_var$end)

            D <- nchar(del)
            I <- nchar(ins)

            P <- min(D, I)
            p <- 0

            while (p < P &&
                substr(del, D - p, D - p) == substr(ins, I - p, I - p)) {

                p <- p + 1

            }

            event_ref$end <- event_ref$end - p
            event_var$end <- event_var$end - p

        }

    }

    df <- data.frame(
        protein_ref_event_start = event_ref$start,
        protein_ref_event_end = event_ref$end,
        protein_ref_length = nchar(ref_seq),
        protein_ref_seq = ref_seq,
        protein_var_event_start = event_var$start,
        protein_var_event_end = event_var$end,
        protein_var_length = nchar(var_seq),
        protein_var_seq = var_seq,
        stringsAsFactors = FALSE)

    mcols(var) <- cbind(mcols(var), DataFrame(df))

    return(var)

}

mapEventToProtein <- function(tx, event)
{

    event_start_0 <- start(flank(event, -1, TRUE))
    event_end_0 <- start(flank(event, -1, FALSE))

    cds <- setNames(restrict(tx, cdsLeft(tx), cdsRight(tx)), "1")
    event <- intersect(event + 1, cds[[1]])
    event <- ranges(reduce(mapToTranscripts(event, cds)))

    if (length(event) == 0) {

        out <- data.frame(start = NA_integer_, end = NA_integer_)

    } else {

        if (event_start_0 == cdsStart(tx)) {

            event_start <- 0L

        } else {

            event_start <- start(event)

        }

        if (event_end_0 == cdsEnd(tx)) {

            event_end <- end(event)

        } else {

            event_end <- end(event) - 1L

        }

        out <- data.frame(
            start = as.integer(floor(event_start / 3)),
            end = as.integer(ceiling(event_end / 3)))

    }

    return(out)

}

## processTranscripts() removes transcripts with invalid CDS
## and ensures that CDS coordinates include the stop codon

processTranscripts <- function(tx, genome, cores)
{

    tx <- tx[!is.na(cdsStart(tx)) & !is.na(cdsEnd(tx))]

    list_tx <- split(tx, seq_along(tx))

    tx_protein <- unlist(mclapply(list_tx,
        getProteinSeq, genome, mc.cores = cores))

    list_processed <- mcmapply(
        findCDS,
        list_tx,
        cdsStart(tx),
        rep(NA_integer_, length(tx)),
        MoreArgs = list(genome = genome),
        SIMPLIFY = FALSE,
        mc.cores = cores)

    processed <- Reduce(c, list_processed)
    names(processed) <- names(tx)

    processed_protein <- unlist(mclapply(list_processed,
        getProteinSeq, genome, mc.cores = cores))

    invalid <- which(is.na(tx_protein) | is.na(processed_protein) |
        substr(tx_protein, 1, 1) != "M" | tx_protein != processed_protein)

    if (length(invalid) > 0) {

        processed <- processed[-invalid]
        msg <- paste0("Excluded ", length(invalid),
            " transcripts with invalid CDS\n",
            paste(names(tx)[invalid], collapse = "\n"))
        warning(msg, call. = FALSE, immediate. = TRUE)

    }

    return(processed)

}

cdsLeft <- function(tx)
{

    mcols(tx)$cdsStart

}

'cdsLeft<-' <- function(tx, value)
{

    mcols(tx)$cdsStart <- value; tx

}

cdsRight <- function(tx)
{

    mcols(tx)$cdsEnd

}

'cdsRight<-' <- function(tx, value)
{

    mcols(tx)$cdsEnd <- value; tx

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

cds <- function(tx)
{

    if (length(tx) > 1) stop("tx must have length 1")
    tx <- restrict(tx, cdsLeft(tx), cdsRight(tx))
    tx <- granges(tx[[1]])
    st <- as.character(strand(tx[1]))
    tx <- sort(tx, decreasing = (st == "-"))

    return(tx)

}
