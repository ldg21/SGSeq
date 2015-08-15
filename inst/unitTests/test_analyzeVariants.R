test_analyzeVariants <- function()
{
  
    sgvc <- analyzeVariants(sgfc_pred)
    SGSeq:::checkIdenticalSummarizedExperiment(sgvc_pred, sgvc)
    
}
