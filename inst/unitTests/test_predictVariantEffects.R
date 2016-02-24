test_predictVariantEffects <- function()
{

    require(BSgenome.Hsapiens.UCSC.hg19)
    seqlevelsStyle(Hsapiens) <- "NCBI"

    current <- predictVariantEffects(sgv_pred, tx, Hsapiens,
        summarize = FALSE, fix_start_codon = FALSE)
    target <- data.frame(
        variantID = c(2L, 2L, 2L),
        txName = c(
            "uc002fjv.3",
            "uc002fjw.3",
            "uc010vot.2"),
        effect = c(
            "CDS:insertion:in-frame",
            "CDS:insertion:in-frame",
            "5p_UTR:insertion:upstream_start"),
        ref_nt_cdsStart = c(573L, 44L, 308L),
        ref_nt_cdsEnd = c(1869L, 1664L, 1412L),
        ref_nt_lastJunction = c(1646L, 1441L, 1189L),
        ref_nt_eventStart = c(661L, 456L, 204L),
        ref_nt_eventEnd = c(661L, 456L, 204L),
        ref_nt_length = c(3821L, 3616L, 3364L),
        var_nt_cdsStart = c(573L, 44L, 263L),
        var_nt_cdsEnd = c(1956L, 1751L, 1499L),
        var_nt_lastJunction = c(1733L, 1528L, 1276L),
        var_nt_eventStart = c(661L, 456L, 204L),
        var_nt_eventEnd = c(748L, 543L, 291L),
        var_nt_length = c(3908L, 3703L, 3451L),
        ref_aa_eventStart = c(29L, 137L, 0L),
        ref_aa_eventEnd = c(30L, 138L, 0L),
        ref_aa_length = c(431L, 539L, 367L),
        var_aa_eventStart = c(29L, 137L, 0L),
        var_aa_eventEnd = c(59L, 167L, 44L),
        var_aa_length = c(460L, 568L, 411L),
        stringsAsFactors = FALSE)
    checkIdentical(target, current)
    
    current <- predictVariantEffects(sgv_ann, tx, Hsapiens,
        summarize = FALSE, fix_start_codon = FALSE)
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
            "5p_UTR/CDS:deletion/insertion:N-terminal_variant",
            "5p_UTR:deletion/insertion:upstream_start",
            "5p_UTR/CDS:deletion/insertion:N-terminal_variant",
            "5p_UTR/CDS:deletion/insertion:N-terminal_variant",
            "5p_UTR/CDS:deletion/insertion:N-terminal_variant",
            "5p_UTR:deletion/insertion:upstream_start"),
        ref_nt_cdsStart = c(573L, 308L, 573L, 44L, 44L, 308L),
        ref_nt_cdsEnd = c(1869L, 1412L, 1869L, 1664L, 1664L, 1412L),
        ref_nt_lastJunction = c(1646L, 1189L, 1646L, 1441L, 1441L, 1189L),
        ref_nt_eventStart = c(0L, 0L, 0L, 0L, 0L, 0L),
        ref_nt_eventEnd = c(589L, 132L, 589L, 384L, 384L, 132L),
        ref_nt_length = c(3821L, 3364L, 3821L, 3616L, 3616L, 3364L),
        var_nt_cdsStart = c(44L, 44L, 308L, 308L, 573L, 573L),
        var_nt_cdsEnd = c(1664L, 1664L, 1412L, 1412L, 1869L, 1869L),
        var_nt_lastJunction = c(1441L, 1441L, 1189L, 1189L, 1646L, 1646L),
        var_nt_eventStart = c(0L, 0L, 0L, 0L, 0L, 0L),
        var_nt_eventEnd = c(384L, 384L, 132L, 132L, 589L, 589L),
        var_nt_length = c(3616L, 3616L, 3364L, 3364L, 3821L, 3821L),
        ref_aa_eventStart = c(0L, 0L, 0L, 0L, 0L, 0L),
        ref_aa_eventEnd = c(6L, 0L, 64L, 172L, 114L, 0L),
        ref_aa_length = c(431L, 367L, 431L, 539L, 539L, 367L),
        var_aa_eventStart = c(0L, 0L, 0L, 0L, 0L, 0L),
        var_aa_eventEnd = c(114L, 172L, 0L, 0L, 6L, 64L),
        var_aa_length = c(539L, 539L, 367L, 367L, 431L, 431L),
        stringsAsFactors = FALSE)
    checkIdentical(target, current)

}
