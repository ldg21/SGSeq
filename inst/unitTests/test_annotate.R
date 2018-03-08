test_annotate <- function()
{

    sgf <- annotate(sgf_ann, txf_ann)
    target <- sgf_ann
    current <- sgf
    ## target <- as.data.frame(target)
    ## current <- as.data.frame(current)
    checkIdentical(target, current)

    sgv <- annotate(sgv_ann, txf_ann)
    target <- sgv_ann
    current <- sgv
    ## target <- as.data.frame(target)
    ## current <- as.data.frame(current)
    checkIdentical(target, current)

}
