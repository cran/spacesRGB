\name{Digital Cinema Distribution Master}
\alias{DCDM.EOTF}


\docType{data}

\title{
Digital Cinema Distribution Master, the EOTF of DCDM
}

\description{
\tabular{ll}{
\code{DCDM.EOTF}        \tab the standardized DCDM transfer function \cr
}
}


\details{
This is a \code{TransferFunction} designed to be applied to XYZ, instead of the usual RGB.
The electrical encoding of XYZ is denoted X'Y'Z'.
The EOTF is:
\deqn{ X = (52.37/48) * (X')^{2.6} }
and similarly for Y and Z.
}



\references{
SMPTE Standard RP 431-2.
D-Cinema Quality - Reference Projector and Environment for the Display of DCDM in Review Rooms and Theaters.
2011.
}


\seealso{
\code{\link{TransferFunction}} 
}



\keyword{datasets}
