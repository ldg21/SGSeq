test_predictTxFeatures <- function()
{

    path <- system.file("extdata", package = "SGSeq")
    si$file_bam <- file.path(path, "bams", si$file_bam)
    txf <- predictTxFeatures(si, gr)
    target <- txf_pred
    current <- txf
    ## target <- as.data.frame(target)
    ## current <- as.data.frame(current)
    checkIdentical(target, current)

}
