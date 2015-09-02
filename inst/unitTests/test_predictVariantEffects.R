test_predictVariantEffects <- function()
{

    require(BSgenome.Hsapiens.UCSC.hg19)
    seqlevelsStyle(Hsapiens) <- "NCBI"

    current <- predictVariantEffects(sgv_pred, tx,
        Hsapiens, summarize = FALSE)
    target <- data.frame(
        variantID = c(2L, 2L, 2L),
        txName = c(
            "uc002fjv.3",
            "uc002fjw.3",
            "uc010vot.2"),
        effect = c(
            "in-frame_insertion",
            "in-frame_insertion",
            "upstream_start_codon"),
        ref_start = c(29, 137, 0),
        ref_end = c(30, 138, 0),
        ref_length = c(431L, 539L, 367L),
        var_start = c(29, 137, 0),
        var_end = c(59, 167, 44),
        var_length = c(460L, 568L, 411L),
        stringsAsFactors = FALSE)
    checkIdentical(target, current[names(target)])
    
    current <- predictVariantEffects(sgv_ann, tx,
        Hsapiens, summarize = FALSE)
    target <- data.frame(
        variantID = c(1L, 1L, 2L, 2L, 3L, 3L),
        txName = c(
            "uc002fjv.3",
            "uc010vot.2",
            "uc002fjv.3",
            "uc002fjw.3",
            "uc002fjw.3",
            "uc010vot.2"),
        effect = c(
            "N-terminal_variant",
            "upstream_start_codon",
            "N-terminal_variant",
            "N-terminal_variant",
            "N-terminal_variant",
            "upstream_start_codon"),
        ref_start = c(0, 0, 0, 0, 0, 0),
        ref_end = c(6, 0, 64, 172, 114, 0),
        ref_length = c(431L, 367L, 431L, 539L, 539L, 367L),
        var_start = c(0, 0, 0, 0, 0, 0),
        var_end = c(114, 172, 0, 0, 6, 64),
        var_length = c(539L, 539L, 367L, 367L, 431L, 431L),
        stringsAsFactors = FALSE)
    checkIdentical(target, current[names(target)])

}
