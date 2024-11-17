\name{Opto-Optical Transfer Function, parameterized}
\alias{general.OOTF}


\title{
Opto-Optical Transfer Function, general
}

\description{
This parameterized OOTF maps from ACES (linear scene) RGB to linear display RGB
(both of these are optical in nature).
\cr
This transform bypasses non-linear signal display RGB (which is electrical in nature).
}


\usage{
general.OOTF( display_pri, Ymin=0.00010, Ymid=7.2, Ymax=108,
                 observerWP=NULL, limiting_pri=NULL,
                 surround='dark', dynrange='SDR', glowmod='1.1', redmod='1.1' )
}

\arguments{
\item{display_pri}{a 4x2 matrix containing the display primaries, or a numeric vector of length 8
that can be converted to such a matrix, by row.
Some built-in matrices are \code{\link{REC709_PRI}}, etc.
This argument cannot be \code{NULL}.
}

\item{Ymin}{the minimum display luminance, in \eqn{cd/m^2}, or nit.}

\item{Ymid}{the middle display luminance, in \eqn{cd/m^2}, or nit.}

\item{Ymax}{the maximum display luminance, in \eqn{cd/m^2}, or nit.}

\item{observerWP}{the xy chromaticity of the assumed observer whitepoint.
This is used to make a Chromatic Adaptation Transform (CAT) from the ACES whitepoint (approximately D60)
to the assumed observer whitepoint.
If \code{observerWP} is \code{NULL}, it is taken from \code{display_pri}.
If ACES and observer whitepoints are the same, there is no CAT.
}

\item{limiting_pri}{a 4x2 matrix containing the limiting primaries, or a numeric vector of length 8
that can be converted to such a matrix, by row.
If \code{limiting_pri} is not \code{NULL}, and not equal to \code{display_pri},
then the output RGB is clamped to the RGB cube that corresponds to \code{limiting_pri}. 
}

\item{surround}{The level of the surround luminance.
Valid values are \code{'dark'} and \code{'dim'}.
If the level is \code{'dark'} there is no special color compensation.
Partial matching is enabled and matching is case-insensitive.
}

\item{dynrange}{the dynamic range of the display system.
Valid values are \code{'SDR'} (standard dynamic range) and \code{'HDR'} (high dynamic range).
If the value is \code{'HDR'} then \code{surround} is ignored.
Matching is partial and case-insensitive.
}

\item{glowmod}{the version of the Glow Modifier to use.  
The only version currently supported is \code{"1.1"}.
\cr
\code{glowmod} can also be \code{NULL}, \code{NA}, or \code{FALSE},
which means to use no Glow Modifier at all.
}

\item{redmod}{the version of the Red Modifier to use.  
The only version currently supported is \code{"1.1"}.
This string can also be \code{"1.1+pinv"}
which means to use a precision inverse; the forward transfer is exactly the same.
This precision inverse uses an iterative root-finder, 
and is slower than the approximate default inverse.
\cr
\code{redmod} can also be \code{NULL}, \code{NA}, or \code{FALSE},
which means to use no Red Modifier at all.
}
}

\value{
\code{general.OOTF()} returns a \code{TransferFunction} of dimension 3
that maps ACES RGB to linear display RGB.
\cr
The domain of the returned \code{TransferFunction} depends on the values of 
\code{Ymin}, \code{Ymid}, and \code{Ymax}.
The range is [0,1]\eqn{^3}, for which clamping may be used.
\cr
The \code{\link{metadata}} contains the display primaries and whitepoint,
which is useful in \code{\link{installRGB}()}.
}

\details{
The transfer is complicated; here is a summary of the steps:
\enumerate{
\item glow module  (see argument \code{glowmod})
\item red modifier (see argument \code{redmod})
\item matrix conversion from AP0 RGB \rarrow AP1 RGB
\item clamp to non-negative RGB
\item global desaturation  (as in \code{\link{RRT.TF}})
\item single-stage tone-scale  (SSTS)  using \code{Ymin}, \code{Ymid}, and \code{Ymax}
\item absolute luminance to linear code-value, in cube [0,1]\eqn{^3}
\item matrix conversion from AP1 RGB to XYZ
\item dim surround compensation (optional, see arguments \code{surround} and \code{dynrange})
\item clamp XYZ to limiting primaries (optional, see argument \code{limiting_pri})
\item adapt XYZ from ACES whitepoint to observer whitepoint (optional, see argument \code{observerWP})
\item matrix conversion from XYZ to linear display RGB (see argument \code{display_pri})
\item scale and roll-white to avoid clipping (optional, only when \code{observerWP} is ACES whitepoint and display whitepoint is D65 or DCI whitepoint)
\item clamp to non-negative RGB
}    
}


\source{
This function was based on source code at:
\url{https://github.com/ampas/aces-core};
especially the file \code{ACESlib.OutputTransforms.ctl}.
This transform is a sub-transform of the function \code{outputTransform()};
it omits the final EOTF\eqn{^{-1}} and optional Full-to-SMPTE range.
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
\link{Standard Primaries}
}
