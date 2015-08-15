test_getSGVariantCounts <- function()
{

    sgvc_from_sgfc <- getSGVariantCounts(sgv_pred, sgfc_pred)
    SGSeq:::checkIdenticalSummarizedExperiment(sgvc_pred, sgvc_from_sgfc)
    
    path <- system.file("extdata", package = "SGSeq")
    si$file_bam <- file.path(path, "bams", si$file_bam)
    sgvc_from_bam <- getSGVariantCounts(sgv_pred,
      features = sgf_pred, sample_info = si)
    SGSeq:::checkIdenticalSummarizedExperiment(
        sgvc_pred_from_bam, sgvc_from_bam)

    assays(sgvc_from_bam)$countsVariant <- NULL
    assays(sgvc_from_bam)$countsTotal <- NULL
    SGSeq:::checkIdenticalSummarizedExperiment(sgvc_pred, sgvc_from_bam)
    
}
