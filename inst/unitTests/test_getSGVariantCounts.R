test_getSGVariantCounts <- function()
{

    sgvc_from_sgfc <- getSGVariantCounts(sgv_pred, feature_counts = sgfc_pred)
    SGSeq:::checkIdenticalSummarizedExperiment(sgvc_pred, sgvc_from_sgfc)

    path <- system.file("extdata", package = "SGSeq")
    si$file_bam <- file.path(path, "bams", si$file_bam)
    sgvc_from_bam <- getSGVariantCounts(sgv_pred, sample_info = si)
    SGSeq:::checkIdenticalSummarizedExperiment(
        sgvc_pred_from_bam, sgvc_from_bam)

    assays(sgvc_from_bam)$counts <- NULL
    SGSeq:::checkIdenticalSummarizedExperiment(sgvc_from_sgfc, sgvc_from_bam)

}
