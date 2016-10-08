test_findSGVariants <- function()
{

    ## example data
    sgv <- findSGVariants(sgf_pred)
    checkIdentical(sgv_pred, sgv)

    ## complex event 1
    txf <- TxFeatures(c(
        GRanges("1", IRanges(1, 100), "+", type = "F"),
        GRanges("1", IRanges(201, 300), "+", type = "I"),
        GRanges("1", IRanges(401, 500), "+", type = "I"),
        GRanges("1", IRanges(601, 700), "+", type = "L"),
        GRanges("1", IRanges(100, 201), "+", type = "J"),
        GRanges("1", IRanges(100, 401), "+", type = "J"),
        GRanges("1", IRanges(100, 601), "+", type = "J"),
        GRanges("1", IRanges(300, 401), "+", type = "J"),
        GRanges("1", IRanges(300, 601), "+", type = "J"),
        GRanges("1", IRanges(500, 601), "+", type = "J")))
    sgf <- convertToSGFeatures(txf)
    sgv <- findSGVariants(sgf, include = "all")
    target <- DataFrame(
        from = c(
            "D:1:100:+", "D:1:100:+", "D:1:100:+", "D:1:100:+",
            "D:1:100:+", "D:1:100:+", "D:1:300:+", "D:1:300:+"),
        to = c(
            "A:1:401:+", "A:1:401:+", "A:1:601:+", "A:1:601:+",
            "A:1:601:+", "A:1:601:+", "A:1:601:+", "A:1:601:+"),
        type = c(
            "J", "JEJ", "J", "JEJ", "JEJ", "JEJEJ", "J", "JEJ"),
        featureID = c(
            "4", "3,7,9", "5", "3,7,10",
            "4,12,14", "3,7,9,12,14", "10", "9,12,14"),
        segmentID = c(
            "6", "2,3", "7", "2,8", "6,4", "2,3,4", "8", "3,4"),
        closed5p = c(
            TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE),
        closed3p = c(
            TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE),
        closed5pEvent = c(
            TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
        closed3pEvent = c(
            FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
        geneID = rep(1L, 8),
        eventID = c(
            1L, 1L, 2L, 2L, 2L, 2L, 3L, 3L),
        variantID = 1:8,
        featureID5p = IntegerList(
            4, NULL, 5, NULL, 4, NULL, 10, 9),
        featureID3p = IntegerList(
            4, 9, 5, 10, NULL, NULL, 10, NULL),
        featureID5pEvent = IntegerList(
            NULL, NULL, c(5, 3, 4), c(5, 3, 4),
            c(5, 3, 4), c(5, 3, 4), c(10, 9), c(10, 9)),
        featureID3pEvent = IntegerList(
            c(4, 9), c(4, 9), c(5, 10, 14), c(5, 10, 14),
            c(5, 10, 14), c(5, 10, 14), NULL, NULL))
    current <- mcols(sgv)[names(target)]
    checkIdentical(target, current)

    ## complex event 2
    txf <- TxFeatures(c(
        GRanges("1", IRanges(1, 100), "+", type = "F"),
        GRanges("1", IRanges(201, 300), "+", type = "I"),
        GRanges("1", IRanges(401, 500), "+", type = "I"),
        GRanges("1", IRanges(601, 700), "+", type = "L"),
        GRanges("1", IRanges(801, 900), "+", type = "L"),
        GRanges("1", IRanges(1001, 1100), "+", type = "L"),
        GRanges("1", IRanges(100, 201), "+", type = "J"),
        GRanges("1", IRanges(100, 401), "+", type = "J"),
        GRanges("1", IRanges(100, 601), "+", type = "J"),
        GRanges("1", IRanges(300, 601), "+", type = "J"),
        GRanges("1", IRanges(300, 801), "+", type = "J"),
        GRanges("1", IRanges(500, 601), "+", type = "J"),
        GRanges("1", IRanges(500, 1001), "+", type = "J")))
    sgf <- convertToSGFeatures(txf)
    sgv <- findSGVariants(sgf, include = "all")
    target <- DataFrame(
        from =  c(
            "D:1:100:+", "D:1:100:+", "D:1:100:+", "D:1:100:+",
            "D:1:100:+", "D:1:100:+", "D:1:100:+", "D:1:100:+",
            "D:1:300:+", "D:1:300:+", "D:1:500:+", "D:1:500:+"),
        to = c(
            "A:1:601:+", "A:1:601:+", "A:1:601:+", "E:1:700:+",
            "E:1:700:+", "E:1:700:+", "E:1:900:+", "E:1:1100:+",
            "E:1:700:+", "E:1:900:+", "E:1:700:+", "E:1:1100:+"),
        type = c(
            "J", "JEJ", "JEJ", "JE", "JEJE", "JEJE",
            "JEJE", "JEJE", "JE", "JE", "JE", "JE"),
        featureID = c(
            "5", "3,7,9", "4,12,14", "5,17", "3,7,9,17", "4,12,14,17",
            "3,7,10,19", "4,12,15,21", "9,17", "10,19", "14,17", "15,21"),
        segmentID = c(
            "9", "2,4", "3,6", "9,8", "2,4,8", "3,6,8",
            "2,5", "3,7", "4,8", "5", "6,8", "7"),
        closed5p = c(
            TRUE, TRUE, TRUE, FALSE, FALSE, FALSE,
            TRUE, TRUE, FALSE, TRUE, FALSE, TRUE),
        closed3p = c(
            TRUE, FALSE, FALSE, TRUE, FALSE, FALSE,
            FALSE, FALSE, TRUE, TRUE, TRUE, TRUE),
        closed5pEvent = c(
            TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
            TRUE, TRUE, FALSE, FALSE, FALSE, FALSE),
        closed3pEvent = c(
            FALSE, FALSE, FALSE, TRUE, TRUE, TRUE,
            TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
        geneID = rep(1L, 12),
        eventID = c(
            1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 4L, 4L),
        variantID = 1:12,
        featureID5p = IntegerList(
            5, NULL, NULL, 5, NULL ,NULL,
            NULL, NULL, 9, 10, 14, 15),
        featureID3p = IntegerList(
            5, 9, 14, NULL, NULL, NULL,
            NULL, NULL, NULL, NULL, NULL, NULL),
        featureID5pEvent = IntegerList(
            NULL, NULL, NULL, c(5, 3, 4), c(5, 3, 4), c(5, 3, 4),
            c(5, 3, 4), c(5, 3, 4), c(9, 10), c(9, 10), c(14, 15), c(14, 15)),
        featureID3pEvent = IntegerList(
            c(5, 9, 14), c(5, 9, 14), c(5, 9, 14), NULL, NULL, NULL,
            NULL, NULL, NULL, NULL, NULL, NULL))
    current <- mcols(sgv)[names(target)]
    checkIdentical(target, current)

    ## complex event 3
    txf <- TxFeatures(c(
        GRanges("1", IRanges(1, 100), "+", type = "F"),
        GRanges("1", IRanges(201, 300), "+", type = "F"),
        GRanges("1", IRanges(401, 500), "+", type = "F"),
        GRanges("1", IRanges(601, 700), "+", type = "I"),
        GRanges("1", IRanges(801, 900), "+", type = "I"),
        GRanges("1", IRanges(1001, 1100), "+", type = "L"),
        GRanges("1", IRanges(100, 601), "+", type = "J"),
        GRanges("1", IRanges(300, 801), "+", type = "J"),
        GRanges("1", IRanges(500, 601), "+", type = "J"),
        GRanges("1", IRanges(500, 801), "+", type = "J"),
        GRanges("1", IRanges(500, 1001), "+", type = "J"),
        GRanges("1", IRanges(700, 1001), "+", type = "J"),
        GRanges("1", IRanges(900, 1001), "+", type = "J")))
    sgf <- convertToSGFeatures(txf)
    sgv <- findSGVariants(sgf, include = "all")
    target <- DataFrame(
        from =  c(
            "S:1:1:+", "S:1:1:+", "S:1:201:+", "S:1:201:+",
            "S:1:401:+", "S:1:401:+", "S:1:401:+", "S:1:401:+",
            "S:1:401:+", "D:1:500:+", "D:1:500:+", "D:1:500:+"),
       to = c(
           "A:1:601:+", "A:1:1001:+", "A:1:801:+", "A:1:1001:+",
           "A:1:601:+", "A:1:801:+", "A:1:1001:+", "A:1:1001:+",
           "A:1:1001:+", "A:1:1001:+", "A:1:1001:+", "A:1:1001:+"),
       type = c(
           "EJ", "EJEJ", "EJ", "EJEJ", "EJ", "EJ",
           "EJ", "EJEJ", "EJEJ", "J", "JEJ", "JEJ"),
       featureID = c(
           "1,3", "1,3,13,15", "4,6", "4,6,17,19", "7,9", "7,10",
           "7,11", "7,9,13,15", "7,10,17,19", "11", "9,13,15", "10,17,19"),
        segmentID = c(
            "1", "1,6", "2", "2,7", "3,4", "3,5",
            "3,9", "3,4,6", "3,5,7", "9", "4,6", "5,7"),
        closed5p = c(
            TRUE, FALSE, TRUE, FALSE, TRUE, TRUE,
            TRUE, FALSE, FALSE, TRUE, FALSE, FALSE),
        closed3p = c(
            TRUE, TRUE, TRUE, TRUE, FALSE, FALSE,
            FALSE, FALSE, FALSE, TRUE, TRUE, TRUE),
        closed5pEvent = c(
            TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
            TRUE, TRUE, TRUE, FALSE, FALSE, FALSE),
        closed3pEvent = c(
            FALSE, TRUE, FALSE, TRUE, FALSE, FALSE,
            TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
        geneID = rep(1L, 12),
        eventID = c(
            1L, 2L, 3L, 2L, 1L, 3L, 2L, 2L, 2L, 4L, 4L, 4L),
        variantID = 1:12,
        featureID5p = IntegerList(
            NULL, NULL, NULL, NULL, NULL ,NULL,
            NULL, NULL, NULL, 11, 9, 10),
        featureID3p = IntegerList(
            3, NULL, 6, NULL, 9, 10,
            11, NULL, NULL, 11, NULL, NULL),          
        featureID5pEvent = IntegerList(
            NULL, NULL, NULL, NULL, NULL ,NULL,
            NULL, NULL, NULL, c(11, 9, 10), c(11, 9, 10), c(11, 9, 10)),
        featureID3pEvent = IntegerList(
            c(3, 9), c(15, 19, 11), c(6, 10), c(15, 19, 11), c(3, 9), c(6, 10),
            c(15, 19, 11), c(15, 19, 11), c(15, 19, 11), NULL, NULL, NULL))
    current <- mcols(sgv)[names(target)]
    checkIdentical(target, current)

}
