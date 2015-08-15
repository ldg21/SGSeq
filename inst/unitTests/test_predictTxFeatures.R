test_predictTxFeatures <- function()
{

    path <- system.file("extdata", package = "SGSeq")
    si$file_bam <- file.path(path, "bams", si$file_bam)
    txf <- predictTxFeatures(si, gr)
    checkIdentical(txf_pred, txf)

}
