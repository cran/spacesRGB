\name{Full Range to SMPTE Range}
\alias{FullRangeToSMPTE.TF}

\docType{data}

\title{
Full Range to SMPTE Range
}

\description{
\tabular{ll}{
\code{FullRangeToSMPTE.TF}   \tab   the standardized SMPTE range transfer function \cr
}
}


\details{
This is a \code{TransferFunction} object that maps from non-linear display signal RGB to itself.
It maps from the full range [0,1] to the 
smaller range [ 64/1023, 940/1023 ] \eqn{\approx}  [ 0.06256109, 0.9188661 ].
The latter is the 10-bit "legal-SMPTE"" range.
It does this in an affine way, and in fact simply uses \code{\link{affine.TF}()}.
}

\seealso{
\code{\link{affine.TF}()},
\code{\link{TransferFunction}}
}

\keyword{datasets}
