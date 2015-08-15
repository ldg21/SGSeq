test_annotate <- function()
{

    sgf <- annotate(sgf_ann, txf_ann)
    checkIdentical(sgf_ann, sgf)

    sgv <- annotate(sgv_ann, txf_ann)
    checkIdentical(sgv_ann, sgv)

}
