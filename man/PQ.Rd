\name{Perceptual Quality Transform}
\alias{PQ.EOTF}


\title{
Perceptual Quality Transform
}

\description{
The Perceptual Quantizer is a transfer function that allows for the display of high dynamic range (HDR) video.
}


\usage{
PQ.EOTF( Lmax=10000 )
}

\arguments{
\item{Lmax}{the maximum luminance, in \eqn{cd/m^2}, or nit.
}
}


\value{
\code{PQ.EOTF()} returns a \code{TransferFunction} that maps signal-display RGB to linear-display RGB.
The interval [0,1] maps to [0,\code{Lmax}].
}


\references{
ST-2084.
SMPTE Standard - High Dynamic Range Electro-Optical Transfer Function of Mastering Reference Displays.
2014.
}

\seealso{
\code{\link{TransferFunction}}
}
