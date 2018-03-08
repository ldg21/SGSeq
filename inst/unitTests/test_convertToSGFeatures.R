test_convertToSGFeatures <- function()
{

    sgf <- convertToSGFeatures(txf_ann)
    target <- sgf_ann
    current <- sgf
    ## target <- as.data.frame(target)
    ## current <- as.data.frame(current)
    checkIdentical(target, current)
    
}
