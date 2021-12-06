\name{Standard Primaries}
\alias{Standard Primaries}
\alias{AP0_PRI}
\alias{AP1_PRI}
\alias{REC709_PRI}
\alias{REC2020_PRI}
\alias{P3D65_PRI}
\alias{P3D60_PRI}
\alias{P3DCI_PRI}

\docType{data}

\title{
Standard Primaries
}

\description{
xy Chromaticities for some standard primary sets.
These include Red, Green, Blue, and White.

\tabular{ll}{
\code{AP0_PRI}          \tab ACES Scene-Referred Primaries, from SMPTE ST2065-1\cr
\code{AP1_PRI}          \tab working space and rendering primaries for ACES 1.0\cr
\code{REC709_PRI}       \tab Rec.709  (aka BT.709) primaries\cr
\code{REC2020_PRI}      \tab Rec.2020  (aka BT.2020) primaries\cr
\code{P3D65_PRI}        \tab RGB primaries from DCI-P3, with D65 for the whitepoint \cr
\code{P3D60_PRI}        \tab RGB primaries from DCI-P3, with ACES whitepoint (approximately D60)\cr
\code{P3DCI_PRI}        \tab RGB primaries from DCI-P3, with DCI whitepoint\cr
}
}



\details{
All of these are built-in 4x2 matrices, with xy coordinates in the rows,
and in RGBW order.
}



\references{
ST 2065-1:2012.
SMPTE Standard - Academy Color Encoding Specification (ACES).
2013.

SMPTE Standard RP 431-2.
D-Cinema Quality - Reference Projector and Environment for the Display of DCDM in Review Rooms and Theaters.
2011.

Wikipedia. \bold{DCI-P3}.
\url{https://en.wikipedia.org/wiki/DCI-P3}.



BT.709.
Parameter values for the HDTV standards for production and international programme exchange.
June 2015.

BT.2020.
Parameter values for ultra-high definition television systems for production and international programme exchange.
October 2015.

}




\examples{
AP0_PRI
##          x        y
##  R 0.73470  0.26530
##  G 0.00000  1.00000
##  B 0.00010 -0.07700
##  W 0.32168  0.33767
}

\keyword{datasets}
