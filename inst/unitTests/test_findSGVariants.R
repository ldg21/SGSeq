test_findSGVariants <- function()
{

    sgv <- findSGVariants(sgf_pred)
    checkIdentical(sgv_pred, sgv)
    
}
