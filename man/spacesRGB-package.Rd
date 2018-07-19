% File colorSpec-package.Rd
\docType{package}
\name{spacesRGB-package}
\alias{spacesRGB-package}
\title{Calculate with RGB spaces}


\description{

RGB spaces are confusing because for each RGB space, up to three RGB vectors are involved:
\itemize{
\item linear \emph{scene RGB} that corresponds to light from a scene
\item non-linear \emph{signal RGB} encoded from \emph{scene RGB}
\item linear \emph{display RGB} that corresponds to light emitted from a display
}
Without qualification, the non-linear \emph{signal RGB} is usually meant.
This is the case for \code{grDevices::\link{rgb}()},
\code{grDevices::\link{hsv}()},
\code{grDevices::\link{convertColor}()}, etc.
This diagram shows the full processing system.

\if{html}{\figure{diagram.png}{options: width=600}}
\if{latex}{\figure{diagram.pdf}{options: width=5in}}

There are two XYZ vectors.
\emph{Scene XYZ} is converted to \emph{scene RGB} by a 3x3 matrix \eqn{M}.
This matrix is uniquely determined by the xy chromaticities of the RGB primaries
and the XYZ of the whitepoint (3*2 + 3 = 9 total numbers).
\emph{Display RGB} is converted to \emph{display XYZ} through multiplication by \eqn{M^{-1}}.
In this package \eqn{M} and \eqn{M^{-1}} are computed with full precision;
the matrices that may be published in the standard are not used.
% For the pre-installed spaces, the matrices match those computed by \cite{Lindbloom}.
All 4 vectors - \emph{scene RGB}, \emph{scene XYZ},
\emph{display RGB}, and \emph{display XYZ} - 
are called linear because they are all linear functions of the spectral energy distrubution of light.

The terminology below is mostly taken from \cite{BT.709},
\cite{BT.1886}, and \cite{BT.2100} 

In this package, both linear RGB vectors are normalized to the cube [0,1]\eqn{^3};
they are usually thought of as \emph{optical} in nature.
\cr
Signal RGB can be thought of as \emph{electrical} in nature.
It is common to think of this signal as an 8-bit number stored in the frame buffer,
so the cube [0,255]\eqn{^3} is very commmon for signal RGB.
This package allows the user to specify any positive number as the upper limit of the interval.
The relevant function argument is \code{maxSignal}.


The conversion function from optical scene RGB (linear) 
to electrical signal RGB (non-linear) is called the
\emph{Opto-Electronic Transfer Function} (OETF).
Each primary component (R,G, and B) is transformed independently and has the same OETF.
In this package, the OETF is always defined as a map of [0,1] to itself.
The classical OETF is the "\eqn{1/\gamma} power law".
For Adobe RGB, the OETF is classical with \eqn{\gamma = 563/256 \approx 2.2}.
For the most popular RGB space - sRGB - the OETF is a continuous function defined in 2 pieces.
In ICC profiles, it is called a \emph{Tone Response Curve} (TRC), or \emph{shaper curve}.

All the transfer functions in this package are continuous and strictly increasing functions
that take [0,1] to itself.
In fancy terminology they are \emph{orientation preserving homeomorphisms} of [0,1].
All transfer functions are either power laws, or well-approximated by power laws.
The exponent giving the best match is called the
\emph{effective gamma}, the
\emph{approximate gamma}, or the 
\emph{best-fit exponent}
for the transfer function.
The function \code{\link{summaryRGB}()} displays the \emph{best-fit exponent} for
the transfer functions of all the installed RGB spaces.

The conversion function from electrical signal RGB (non-linear) 
to optical display RGB (linear) is called the
\emph{Electro-Optical Transfer Function} (EOTF).
Each primary component (R,G, and B) is transformed independently and has the same EOTF.
In this package, the EOTF is always defined as a map of [0,1] to itself.
The classical EOTF is the "\eqn{\gamma} power law".
That is \eqn{R' = R^\gamma}, where \eqn{R'} is display R and \eqn{R} is signal R.

All the RGB spaces in \cite{Lindbloom} and all the pre-installed RGB spaces in this package
(except for the space \code{'HD+2.4'}),
have EOTF = OETF\eqn{^{-1}}.
This implies that scene RGB and display RGB are the same,
and that scene XYZ and display XYZ are the same.
In \cite{SMPTE 240M}, the OETF is the transfer function of the \emph{reference camera},
and the EOTF is the transfer function of the \emph{reference reproducer},
and they are inverses.
Stated another way the Opto-optical Transfer Function (OOTF),
which is the composition of the OETF and the EOTF,
\deqn{OOTF = EOTF \circ OETF}
is the identity, i.e. trivial.
In some contexts one says that the \emph{system gamma} is 1,
or that the \emph{end-to-end exponent} is 1.
However, due to the \emph{Hunt effect} and the \emph{Stevens effect},
color reproduction is more accurate when the the system gamma is *not* 1;
see \cite{Color Appearance Model}.
See \bold{Examples} for creation of \code{'HD+2.4'} - an RGB space with a non-trivial OOTF.

The following RGB space features are outside the scope of this package:
\itemize{
\item black offset, as in \cite{BT.1886}
\item artistic adjustments, as in \cite{BT.2100}
\item high-dynamic range, as in \cite{BT.2100}
\item Academy Color Encoding System  (ACES)
}

An excellent reference for all this is \cite{Poynton}.
}


\section{The API}{
The API is small.  There are functions to
\itemize{
\item convert between signal RGB and XYZ (both scene and display)
\item convert between signal RGB and linear RGB (both scene and display)
\item install, uninstall, and query the dictionary of RGB spaces
}
The XYZ vector is always relative to the XYZ of the white point of the RGB space;
there are no functions for chromatic adaptation, as in \cite{Pascale}.
% The XYZ of the white-point always has Y=1.

It is also possible to convert between XYZ and linear RGB in one function call
by overriding the OETF and EOTF by setting \code{TF=1} as follows \cr
\code{RGBfromXYZ( XYZ, TF=1 )}\cr
\code{XYZfromRGB( XYZ, TF=1 )}\cr
In the first case, besides being more efficient, there is also no clamping to the unit cube,
so one can tell how far outside the cube that the linear RGB really is.
}



\section{Note}{
All the RGB spaces are stored in a dictionary that initially contains
\code{sRGB}, \code{AdobeRGB},  \code{ProPhotoRGB},
\code{AppleRGB},   \code{BT.709},      \code{BT.2020},   \code{240M},   and \code{HD+2.4}.
\cr
There are some close similarities.
\cr
\code{BT.709} (created 1990) and \code{sRGB} (created 1996)
have the same primaries and whitepoint, but different transfer functions.
\code{sRGB} copied its primaries and whitepoint from \code{BT.709}.
\cr
\code{BT.709} and \code{BT.2020} have very close transfer functions, differing only in precision.
\cr
\code{HD+2.4} has the same primaries and OETF as \code{BT.709}, but a different EOTF.
}

\section{Other Packages}{
Package \pkg{colorspace} has similar functionality,
and is much faster because it is compiled.
But it only supports one RGB space - \bold{sRGB}.
\cr
Package \pkg{colorscience} has very flexible conversion functions,
which include \code{gamma} override and chromatic adaptation,
and it has all the RGB spaces in \cite{Lindbloom}.
But is structured so it only allows one space with non-classical EOTF/OETF, namely \bold{sRGB}.
}


\author{
Glenn Davis <gdavis@gluonics.com>
}

\examples{
# Install an RGB space named 'HD+2.4', with encoding from BT.709 and display from BT.1886.
# the OOTF for this space is non-trivial

prim    = matrix( c(0.64,0.33,  0.30,0.60,  0.15,0.06 ), 3, 2, byrow=TRUE )
white   = c( 0.3127, 0.3290 )
    
installRGB( "HD+2.4", prim, white, OETF='BT.709', EOTF=2.4, overwrite=TRUE )
}



\references{
BT.709.
Parameter values for the HDTV standards for production and international programme exchange.
June 2015.

BT.1886.
Reference electro-optical transfer function for flat panel displays used in HDTV studio production.
March 2011.

BT.2020.
Parameter values for ultra-high definition television systems for production and international programme exchange.
October 2015.

BT.2100.
Image parameter values for high dynamic range television for use in production
and international programme exchange.
June 2017

Lindbloom, Bruce.
\bold{RGB/XYZ Matrices}.
\url{http://brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html}

Pascale, Danny.
A Review of RGB Color Spaces ...from xyY to R'G'B'.
\url{http://www.babelcolor.com/download/A\%20review\%20of\%20RGB\%20color\%20spaces.pdf}.
2003.

Poynton, Charles.
\bold{Digital Video and HD - Algorithms and Interfaces}.
Morgan Kaufmann.
Second Edition. 2012.


Wikipedia. \bold{Adobe RGB (1998)}.
\url{http://en.wikipedia.org/wiki/Adobe_RGB_color_space}.

Wikipedia. \bold{Color Appearance Model}.
https://en.wikipedia.org/wiki/Color_appearance_model

Wikipedia. \bold{ProPhoto RGB}.
\url{http://en.wikipedia.org/wiki/ProPhoto_RGB_color_space}.

Wikipedia. \bold{sRGB}.
\url{http://en.wikipedia.org/wiki/SRGB}.
}

\seealso{
\code{\link{summaryRGB}()},
\code{\link{installRGB}()}
}



\keyword{package}

