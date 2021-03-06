% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{predictTxFeatures}
\alias{predictTxFeatures}
\title{Splice junction and exon prediction from BAM files}
\usage{
predictTxFeatures(sample_info, which = NULL, alpha = 2, psi = 0,
  beta = 0.2, gamma = 0.2, min_junction_count = NULL, min_anchor = 1,
  max_complexity = 20, min_n_sample = 1, min_overhang = NA,
  verbose = FALSE, cores = 1)
}
\arguments{
\item{sample_info}{Data frame with sample information.
Required columns are \dQuote{sample_name}, \dQuote{file_bam},
\dQuote{paired_end}, \dQuote{read_length}, \dQuote{frag_length}
and \dQuote{lib_size}. Library information can be obtained with
function \code{getBamInfo}.}

\item{which}{\code{GRanges} of genomic regions to be considered for
feature prediction, passed to \code{ScanBamParam}}

\item{alpha}{Minimum FPKM required for a splice junction to be
included. Internally, FPKMs are converted to counts, requiring arguments
\code{read_length}, \code{frag_length} and \code{lib_size}.
\code{alpha} is ignored if argument \code{min_junction_count}
is specified.}

\item{psi}{Minimum splice frequency required for a splice junction
to be included}

\item{beta}{Minimum relative coverage required for an internal exon
to be included}

\item{gamma}{Minimum relative coverage required for a terminal exon
to be included}

\item{min_junction_count}{Minimum fragment count required for a splice
junction to be included. If specified, argument \code{alpha} is ignored.}

\item{min_anchor}{Integer specifiying minimum anchor length}

\item{max_complexity}{Maximum allowed complexity. If a locus exceeds
this threshold, it is skipped, resulting in a warning.
Complexity is defined as the maximum number of unique predicted
splice junctions overlapping a given position.
High complexity regions are often due to spurious read alignments
and can slow down processing. To disable this filter, set to \code{NA}.}

\item{min_n_sample}{Minimum number of samples a feature must be
observed in to be included}

\item{min_overhang}{Minimum overhang required to suppress filtering or
trimming of predicted terminal exons (see the manual page for
\code{processTerminalExons}). Use \code{NULL} to disable processing
(disabling processing is useful if results are subsequently merged
with other predictions and processing is postponed until after the
merging step).}

\item{verbose}{If \code{TRUE}, generate messages indicating progress}

\item{cores}{Number of cores available for parallel processing}
}
\value{
\code{TxFeatures} object
}
\description{
Splice junctions and exons are predicted for each sample and merged
across samples. Terminal exons are filtered and trimmed, if applicable.
For details, see the help pages for
\code{\link{predictTxFeaturesPerSample}}, \code{\link{mergeTxFeatures}},
and \code{\link{processTerminalExons}}.
}
\examples{
path <- system.file("extdata", package = "SGSeq")
si$file_bam <- file.path(path, "bams", si$file_bam)
txf <- predictTxFeatures(si, gr)
}
\author{
Leonard Goldstein
}
