\name{Hybrid Log-Gamma Transform}
\alias{HLG.OETF}
\alias{HLG.OOTF}

\title{
Hybrid Log-Gamma Transform
}

\description{
The Hybrid Log-Gamma OETF is a transfer function that allows for the display of high dynamic range (HDR) video.
The version here is that supported by the 
ACES (Academy Color Encoding System) and 
HEVC (High Efficiency Video Coding) standards.
\cr
For use with ACES, a specialized HLG-based OOTF is provided that references the Perceptual Quality (PQ) EOTF.
It converts the ST.2084 (PQ) output to HLG using the method specified in Section 7 of BT.2390-0. 
}


\usage{
HLG.OETF()
HLG.OOTF( gamma=1.2, Lb=0, Lw=1000 )
}

\arguments{
\item{gamma}{the applied exponent, from scene linear to display linear}
\item{Lb}{the luminance of black, in \eqn{cd/m^2}, or nit.}
\item{Lw}{the luminance of white, in \eqn{cd/m^2}, or nit.}
}


\value{
\code{HLG.OETF()} returns a univariate \code{TransferFunction} that maps
linear scene RGB to signal display RGB.
The interval [0,1] maps to [0,1] (as in the HVEC standard).
\cr
\code{HLG.OOTF()} returns a multivariate \code{TransferFunction} of dimension 3
that maps linear scene RGB to linear display RGB.
It maps the cube [0,1]\eqn{^3} to the cube [\code{Lb},\code{Lw}]\eqn{^3},
but the image is only a proper subset of the cube.
}

\details{
\code{HLG.OOTF()} is 3D and does not operate on each channel independently.
It uses a scaling factor based on these RGB weights - (0.2627,0.6780,0.0593) - from Section 7 of BT.2390.
}




\references{
ST-2084.
SMPTE Standard - High Dynamic Range Electro-Optical Transfer Function of Mastering Reference Displays.
2014.

BT.2390.
High dynamic range television for production and international programme exchange.
April 2018.

H.265 : High Efficiency Video Coding.
https://www.itu.int/rec/T-REC-H.265-201802-I/en.
2018-02-13.
}

\seealso{
\code{\link{TransferFunction}},
\code{\link{PQ.EOTF}}
}
