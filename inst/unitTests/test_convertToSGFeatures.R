test_convertToSGFeatures <- function()
{

    sgf <- convertToSGFeatures(txf_ann)
    checkIdentical(sgf_ann, sgf)
    
}
