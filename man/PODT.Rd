\name{Partial Output Device Transform, parameterized}
\alias{general.PODT}


\title{
Partial Output Device Transform, general 
}

\description{
A partial Output Device Transform (PODT) maps from OCES to linear display RGB
(both of these are optical in nature).
The adjective "partial" is used because this
is an ODT that omits the final OETF (which maps from linear display RGB to signal display RGB).
\cr
This PODT is parameterized.
}


\usage{
general.PODT( display_pri, Ymax=1, observerWP=NULL, surround='dark', limiting_pri=NULL )
}

\arguments{
\item{display_pri}{a 4x2 matrix containing the display primaries, or a numeric vector of length 8
that can be converted to such a matrix, by row.
Some built-in matrices are \code{\link{REC709_PRI}}, etc.
\cr
\code{display_pri} can also be \code{NULL}, which means that the PODT maps to XYZ, instead of RGB.
This is used in the case of DCDM (Digital Cinema Distribution Master).
See the User Guide Appendix for examples of this.
\cr
\code{display_pri}  is stored in the metadata of the returned object and later used
in \code{\link{installRGB}()} (if the PODT is passed in an argument).
}

\item{Ymax}{the maximum luminance of the output device, in \eqn{cd/m^2} (or nits).
This has no effect on the PODT itself.
It is stored in the metadata and later used
in \code{\link{installRGB}()} (if the PODT is passed in an argument) when computing the
3x3 matrix that transforms from display RGB to display XYZ.
}

\item{observerWP}{the xy chromaticity of the assumed observer whitepoint.
This is used to make a Chromatic Adaptation Transform (CAT) from the ACES whitepoint (approximately D60)
to the assumed observer whitepoint.
If \code{observerWP} is \code{NULL}, it is taken from \code{display_pri}.
If \code{display_pri} is \code{NULL}, then it is taken from \code{limiting_pri}.
If \code{limiting_pri} is \code{NULL}, or if two whitepoints are the same, then there is no CAT.
}

\item{surround}{The level of the surround luminance.
Valid values are \code{'dark'} and \code{'dim'}.
If the level is \code{'dark'} there is no special color compensation.
Partial matching is enabled and matching is case-insensitive.
}

\item{limiting_pri}{a 4x2 matrix containing the limiting primaries, or a numeric vector of length 8
that can be converted to such a matrix, by row.
If \code{limiting_pri} is not \code{NULL}, and not equal to \code{display_pri},
then the output RGB is clamped to the RGB cube that corresponds to \code{limiting_pri}. 
}
}

\value{
\code{general.PODT()} returns a \code{TransferFunction} of dimension 3
that maps OCES RGB to linear display RGB.
The domain is [0,10000]\eqn{^3} and the range is [0,1]\eqn{^3}.
\cr
The \code{\link{metadata}} contains the display primaries and whitepoint,
which is useful in \code{\link{installRGB}()}.
}

\details{
The transfer is complicated; here is a summary of the steps:
\enumerate{
\item matrix conversion from AP0 RGB \rarrow AP1 RGB
\item clamp to non-negative RGB
\item segmented spline, assuming \code{CINEMA_WHITE}=48 nit
\item absolute luminance to linear code-value, in cube [0,1]\eqn{^3}
\item scale and roll-white to avoid clipping (optional, only when \code{observerWP} is ACES whitepoint and display whitepoint is D65 or DCI whitepoint)
\item dim surround compensation with conversion to XYZ and back again (optional, see argument \code{surround})
\item matrix conversion from AP1 RGB to XYZ
\item adapt XYZ from ACES whitepoint to observer whitepoint (optional, see argument \code{observerWP})
\item clamp XYZ to limiting primaries (optional, see argument \code{limiting_pri})
\item matrix conversion from XYZ to linear display RGB (but not for DCDM, see argument \code{display_pri})
\item clamp linear display RGB (or XYZ for DCDM) to the cube [0,1]\eqn{^3}
}    
}


\source{
This function was based on source code at:
\url{https://github.com/ampas/aces-dev};
especially the files under the folder \code{aces-dev-master/transforms/ctl/odt/}.
}


\references{
ST 2065-1:2012.
SMPTE Standard - Academy Color Encoding Specification (ACES).
2013.
}

\seealso{
\code{\link{TransferFunction}},
\code{\link{installRGB}()},
\code{\link{metadata}()},
\code{\link{RRT.TF}},
\link{Standard Primaries}
}
