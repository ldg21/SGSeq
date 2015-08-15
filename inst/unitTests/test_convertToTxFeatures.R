test_convertToTxFeatures <- function()
{
  
    gr <- GRanges(1, IRanges(c(1, 201), c(100, 300)), "+")
    grl <- split(gr, 1)
    res <- convertToTxFeatures(grl)
    exp <- TxFeatures(GRanges(1, IRanges(c(1, 100, 201), c(100, 201, 300)),
      "+", type = c("F", "J", "L")), txName = CharacterList(1, 1, 1))
    checkIdentical(res, exp)
    
}
