% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{analyzeFeatures}
\alias{analyzeFeatures}
\title{Analysis of splice graph features from BAM files}
\usage{
analyzeFeatures(sample_info, which = NULL, features = NULL,
  predict = is.null(features), alpha = 2, psi = 0, beta = 0.2,
  gamma = 0.2, min_junction_count = NULL, min_anchor = 1,
  min_n_sample = 1, min_overhang = NA, annotation = NULL,
  max_complexity = 20, verbose = FALSE, cores = 1)
}
\arguments{
\item{sample_info}{Data frame with sample information.
Required columns are \dQuote{sample_name}, \dQuote{file_bam},
\dQuote{paired_end}, \dQuote{read_length}, \dQuote{frag_length}
and \dQuote{lib_size}. Library information can be obtained with
function \code{getBamInfo}.}

\item{which}{\code{GRanges} of genomic regions to be considered for
feature prediction, passed to \code{ScanBamParam}}

\item{features}{\code{TxFeatures} or \code{SGFeatures} object}

\item{predict}{Logical indicating whether transcript
features should be predicted from BAM files}

\item{alpha}{Minimum FPKM required for a splice junction to be included}

\item{psi}{Minimum splice frequency required for a splice junction
to be included}

\item{beta}{Minimum relative coverage required for an internal exon
to be included}

\item{gamma}{Minimum relative coverage required for a terminal exon
to be included}

\item{min_junction_count}{Minimum fragment count required for a splice
junction to be included. If specified, argument \code{alpha} is ignored.}

\item{min_anchor}{Integer specifiying minimum anchor length}

\item{min_n_sample}{Minimum number of samples a feature must be
observed in to be included}

\item{min_overhang}{Minimum overhang required to suppress filtering or
trimming of predicted terminal exons (see the manual page for
\code{processTerminalExons}). Use \code{NULL} to disable processing
(disabling processing is useful if results are subsequently merged
with other predictions and processing is postponed until after the
merging step).}

\item{annotation}{\code{TxFeatures} object used for annotation}

\item{max_complexity}{Maximum allowed complexity. If a locus exceeds
this threshold, it is skipped, resulting in a warning.
Complexity is defined as the maximum number of unique predicted
splice junctions overlapping a given position.
High complexity regions are often due to spurious read alignments
and can slow down processing. To disable this filter, set to \code{NA}.}

\item{verbose}{If \code{TRUE}, generate messages indicating progress}

\item{cores}{Number of cores available for parallel processing}
}
\value{
\code{SGFeatureCounts} object
}
\description{
High-level function for the prediction and quantification of
splice junctions, exon bins and splice sites from BAM files.
}
\details{
Splice junctions and exons are predicted from BAM files with
\code{\link{predictTxFeatures}}.

Known features can be provided as \code{TxFeatures} or
\code{SGFeatures} via argument \code{features}.

If \code{features} is not \code{NULL} and \code{predict} is
\code{TRUE}, known features are augmented with predictions.

Known and/or predicted transcript features are converted to splice
graph features. For details, see \code{\link{convertToSGFeatures}}.

Optionally, splice graph features can be annotated with respect to
a \code{TxFeatures} object provided via argument \code{annotation}.
For details, see the help page for function \code{\link{annotate}}.

Finally, compatible fragment counts for splice graph features are
obtained from BAM files with \code{\link{getSGFeatureCounts}}.
}
\examples{
path <- system.file("extdata", package = "SGSeq")
si$file_bam <- file.path(path, "bams", si$file_bam)
sgfc <- analyzeFeatures(si, gr)
}
\author{
Leonard Goldstein
}
