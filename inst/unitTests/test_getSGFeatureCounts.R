test_getSGFeatureCounts <- function()
{

    path <- system.file("extdata", package = "SGSeq")
    si$file_bam <- file.path(path, "bams", si$file_bam)
    sgfc <- getSGFeatureCounts(si, sgf_pred)
    SGSeq:::checkIdenticalSummarizedExperiment(sgfc_pred, sgfc)
    
}
