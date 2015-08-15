test_TxFeatures <- function()
{

    txf <- TxFeatures()
    checkTrue(is(txf, "TxFeatures"))
    
}

test_SGFeatures <- function()
{

    sgf <- SGFeatures()
    checkTrue(is(sgf, "SGFeatures"))
    
}

test_SGVariants <- function()
{

    sgv <- SGVariants()
    checkTrue(is(sgv, "SGVariants"))
    
}

test_SGFeatureCounts <- function()
{

    sgfc <- SGFeatureCounts()
    checkTrue(is(sgfc, "SGFeatureCounts"))
    
}

test_SGVariantCounts <- function()
{

    sgvc <- SGVariantCounts()
    checkTrue(is(sgvc, "SGVariantCounts"))
    
}
