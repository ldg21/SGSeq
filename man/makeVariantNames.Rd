% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graphs.R
\name{makeVariantNames}
\alias{makeVariantNames}
\title{Create interpretable splice variant names}
\usage{
makeVariantNames(variants)
}
\arguments{
\item{variants}{\code{SGVariants} object}
}
\value{
Character vector with splice variant names
}
\description{
Create interpretable splice variant names
taking format GENE_EVENT_VARIANT/ORDER_TYPE.
GENE is based on geneName if available, and geneID otherwise.
EVENT and VARIANT enumerate events and variants for the same gene
and event, respectively. ORDER indicates the total number of
variants in the same event (e.g. 1/2 refers to the first out of two
splice variants in the event). TYPE is based on variantType.
}
\examples{
makeVariantNames(sgv_pred)
}
\author{
Leonard Goldstein
}
\keyword{internal}
