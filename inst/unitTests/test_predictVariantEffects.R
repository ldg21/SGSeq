test_predictVariantEffects <- function()
{

    require(BSgenome.Hsapiens.UCSC.hg19)
    seqlevelsStyle(Hsapiens) <- "NCBI"

    current <- predictVariantEffects(sgv_pred, tx, Hsapiens, FALSE)
    target <- data.frame(
        variantID = rep(2L, 3),
        txName = c(
            "uc002fjv.3",
            "uc002fjw.3",
            "uc010vot.2"),
        geneName = rep("79791", 3),
        RNA_change = c(
            "r.88_89ins88+1798_88+1884",
            "r.412_413ins412+1798_412+1884",
            "r.-105_-104ins-105+1798_-105+1884"),
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
        geneName = rep("79791", 6),
        RNA_change = c(
            "r.-573_16delins-573-22833_-573-22450",
            "r.-308_-177delins-177+5949_-177+6332",
            "r.-573_16delins[-573-31147_-573-31128;-573-28893_-573-28782]",
            "r.-44_340delins[-44-8314_-44-8295;-44-6060_-44-5949]",
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
