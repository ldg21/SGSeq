test_getBamInfo <- function()
{
  
    path <- system.file("extdata", package = "SGSeq")
    si$file_bam <- file.path(path, "bams", si$file_bam)
    target <- si
    target$frag_length <- c(270, 197, 212, 197, 281.5, 236, 229, 241)
    target$lib_size <- c(384L, 613L, 564L, 675L, 533L, 687L, 475L, 576L)
    current <- getBamInfo(si)
    checkIdentical(target, current)

}
