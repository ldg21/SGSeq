## NOTE
## Sample information for test BAM files must be available.
## All other example R data objects are created. 

library(SGSeq)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

si_file_bam <- si$file_bam

path <- system.file("extdata", package = "SGSeq")
si$file_bam <- file.path(path, "bams", si$file_bam)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
txdb <- keepSeqlevels(txdb, "chr16")
seqlevelsStyle(txdb) <- "NCBI"

gr <- GRanges("16", IRanges(87362942, 87425708), "-")
save(gr, file = "gr.rda")

## annotated

txf_ann <- convertToTxFeatures(txdb)
txf_ann <- txf_ann[txf_ann %over% gr]
save(txf_ann, file = "txf_ann.rda")

sgf_ann <- convertToSGFeatures(txf_ann)
save(sgf_ann, file = "sgf_ann.rda")

sgfc_ann <- getSGFeatureCounts(si, sgf_ann)
sgfc_ann$file_bam <- si_file_bam
save(sgfc_ann, file = "sgfc_ann.rda")

sgv_ann <- findSGVariants(sgf_ann)
save(sgv_ann, file = "sgv_ann.rda")

sgvc_ann <- getSGVariantCounts(sgv_ann, sgfc_ann)
save(sgvc_ann, file = "sgvc_ann.rda")

sgvc_ann_from_bam <- getSGVariantCounts(sgv_ann,
    features = sgf_ann, sample_info = si)
sgvc_ann_from_bam$file_bam <- si_file_bam
save(sgvc_ann_from_bam, file = "sgvc_ann_from_bam.rda")

## predicted

txf_pred <- predictTxFeatures(si, gr)
save(txf_pred, file = "txf_pred.rda")

sgf_pred <- convertToSGFeatures(txf_pred)
save(sgf_pred, file = "sgf_pred.rda")

sgfc_pred <- getSGFeatureCounts(si, sgf_pred)
sgfc_pred$file_bam <- si_file_bam
save(sgfc_pred, file = "sgfc_pred.rda")

sgv_pred <- findSGVariants(sgf_pred)
save(sgv_pred, file = "sgv_pred.rda")

sgvc_pred <- getSGVariantCounts(sgv_pred, sgfc_pred)
save(sgvc_pred, file = "sgvc_pred.rda")

sgvc_pred_from_bam <- getSGVariantCounts(sgv_pred,
    features = sgf_pred, sample_info = si)
sgvc_pred_from_bam$file_bam <- si_file_bam
save(sgvc_pred_from_bam, file = "sgvc_pred_from_bam.rda")

## transcripts

tx <- SGSeq:::convertToTranscripts(txdb)
tx <- tx[tx %over% gr]
save(tx, file = "tx.rda")
