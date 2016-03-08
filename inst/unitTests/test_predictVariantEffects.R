test_predictVariantEffects <- function()
{

    require(BSgenome.Hsapiens.UCSC.hg19)
    seqlevelsStyle(Hsapiens) <- "NCBI"

    current <- predictVariantEffects(sgv_pred, tx, Hsapiens, FALSE)
    target <- data.frame(
        variantID = c(2L, 2L, 2L),
        txName = c(
            "uc002fjv.3",
            "uc002fjw.3",
            "uc010vot.2"),
        RNA_change = c(
            "r.88_89ins88+1798_89-11161",
            "r.412_413ins412+1798_413-11161",
            "r.-105_-104ins-105+1798_-104-11161"),
        RNA_variant_type = c("insertion", "insertion", "insertion"),
        protein_change = c(
            "p.K29_L30insRINPRVKSGRFVKILPDYEHMAYRDVYTC",
            "p.K137_L138insRINPRVKSGRFVKILPDYEHMAYRDVYTC",
            "p.M1ext-44"),
        protein_variant_type = c(
            "in-frame_insertion",
            "in-frame_insertion",
            "N-terminal_extension"),
        stringsAsFactors = FALSE)
    checkIdentical(target, current)

    current <- predictVariantEffects(sgv_ann, tx, Hsapiens, FALSE)
    target <- data.frame(
        variantID = c(1L, 1L, 2L, 2L, 3L, 3L),
        txName = c(
            "uc002fjv.3",
            "uc010vot.2",
            "uc002fjv.3",
            "uc002fjw.3",
            "uc002fjw.3",
            "uc010vot.2"),
        RNA_change = c(
            "r.-573_16delins89-36538_89-36155",
            "r.-308_-177delins-177+5949_-176-23039",
            "r.-573_16delins[89-44852_89-44833;89-42598_89-42487]",
            "r.-44_340delins[341-31736_341-31717;341-29482_341-29371]",
            "r.-44_340delins341-589_341-1",
            "r.-308_-177delins-176-589_-176-1"),
        RNA_variant_type = c(
            "deletion/insertion",
            "deletion/insertion",
            "deletion/insertion",
            "deletion/insertion",
            "deletion/insertion",
            "deletion/insertion"),
        protein_change = c(
            "p.M1_T5delext-113",
            "p.M1ext-172",
            "p.M1_W64del",
            "p.M1_W172del",
            "p.M1_E113delext-5",
            "p.M1ext-64"),
        protein_variant_type = c(
            "N-terminal_variant",
            "N-terminal_extension",
            "N-terminal_deletion",
            "N-terminal_deletion",
            "N-terminal_variant",
            "N-terminal_extension"),
        stringsAsFactors = FALSE)
    checkIdentical(target, current)

}
