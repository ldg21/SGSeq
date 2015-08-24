test_predictVariantConsequences <- function()
{

    require(BSgenome.Hsapiens.UCSC.hg19)
    seqlevelsStyle(Hsapiens) <- "NCBI"

    current <- predictVariantConsequences(sgv_pred, tx,
        Hsapiens, summarize = FALSE)
    target <- data.frame(
        var_id = c(2L, 2L, 2L),
        ref_name = c(
            "uc002fjv.3",
            "uc002fjw.3",
            "uc010vot.2"),
        alt = c(
            "CDS_insertion_in-frame",
            "CDS_insertion_in-frame",
            "CDS_upstream_start"),
        ref_aa_start = c(29, 137, 0),
        ref_aa_end = c(30, 138, 0),
        var_aa_start = c(29, 137, 0),
        var_aa_end = c(59, 167, 44),
        stringsAsFactors = FALSE)
    checkIdentical(target, current[names(target)])
    
    current <- predictVariantConsequences(sgv_ann, tx,
        Hsapiens, summarize = FALSE)
    target <- data.frame(
        var_id = c(1L, 1L, 2L, 2L, 3L, 3L),
        ref_name = c(
            "uc010vot.2",
            "uc002fjv.3",
            "uc002fjw.3",
            "uc002fjv.3",
            "uc002fjw.3",
            "uc010vot.2"),
        alt = c(
            "CDS_upstream_start",
            "CDS_5p_alteration",
            "CDS_5p_alteration",
            "CDS_5p_alteration",
            "CDS_5p_alteration",
            "CDS_upstream_start"),
        ref_aa_start = c(0, 0, 0, 0, 0, 0),
        ref_aa_end = c(0, 6, 172, 64, 114, 0),
        var_aa_start = c(0, 0, 0, 0, 0, 0),
        var_aa_end = c(172, 114, 0, 0, 6, 64),
        stringsAsFactors = FALSE)
    checkIdentical(target, current[names(target)])

}
