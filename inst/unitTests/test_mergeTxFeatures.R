test_mergeTxFeatures <- function()
{

    txf_merged <- mergeTxFeatures(txf_pred, txf_pred)
    target <- txf_pred
    current <- txf_merged
    ## target <- as.data.frame(target)
    ## current <- as.data.frame(current)
    checkIdentical(target, current)
    
}
