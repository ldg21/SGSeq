test_mergeTxFeatures <- function()
{

    txf_merged <- mergeTxFeatures(txf_pred, txf_pred)
    checkIdentical(txf_pred, txf_merged)
    
}
