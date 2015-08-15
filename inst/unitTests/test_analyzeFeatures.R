test_analyzeFeatures <- function()
{

    path <- system.file("extdata", package = "SGSeq")
    si$file_bam <- file.path(path, "bams", si$file_bam)
    sgfc <- analyzeFeatures(si, gr)
    SGSeq:::checkIdenticalSummarizedExperiment(sgfc_pred, sgfc)
    
}
