% File colorSpec-package.Rd
\docType{package}
\name{spacesRGB-package}
\alias{spacesRGB-package}
\title{Calculate with RGB spaces}



\description{
In the field of color management, RGB spaces are broadly divided into 2 groups:
\emph{display-referred} (or \emph{output-referred}) RGB and \emph{scene-referred} RGB.
The former group is more common and that is what this package covers.
For scene-referred RGB (e.g. RIMM RGB) the package will probably work fine
although some function names are not appropriate.

RGB spaces are confusing because for each space two RGB vectors are involved:
a non-linear RGB vector and a linear RGB vector.
Without qualification, the non-linear vector is usually meant.

The non-linear RGB vector refers to the electrical signal sent from the frame buffer to the display
(in modern LCD displays this signal is often modified in the
display so the light flux response is closer to older CRTs).
Non-linear RGB can be thought of as \emph{electrical} RGB,
and we often call it \emph{display} RGB.
It is common to think of this signal as an 8-bit number stored in the frame buffer,
so the interval [0,255] is very commmon for non-linear RGB.
This package allows the user to specify any positive number as the upper limit of the interval.
The relevant function argument is \code{maxValue}.

The linear RGB vector refers to the \emph{light flux} emitted by the display.
The linear R value is proportional to the light flux emitted by the red display primary.
Linear RGB can be thought of as \emph{optical} RGB.
In this package, linear R is always normalized to the interval [0,1], and the same for G and B.
The light emitted by the display also has a CIE XYZ tristimulus vector.
XYZ and linear RGB are linearly related by two 3x3 matrices - one the inverse of the other.
The two matrices are determined by the xy chromaticities of the primaries and the white point.
In this package the two matrices are computed with full precision;
the matrices that may be published in the standard are not used.
For the pre-installed spaces, the matrices match those computed in \cite{Lindbloom}.
It is tempting to call linear RGB \emph{scene referred} RGB, but this is not quite true.
It has been known for some time that luminance affects human color perception,
and that scenes (especially outdoor scenes) are usually brighter than displays.
So straightforward linear scaling is not appropriate to accurately reproduce the colors
of a scene on a display, and an appropriate correction
(sometimes called an \emph{Opto-Optical Conversion Function} (OOCF))
is built-in to the conversion functions between non-linear and linear RGB.

The conversion function from electrical RGB (non-linear) to optical RGB (linear) is called the
\emph{Electro-Optical Conversion Function} (EOCF).
Each primary component (R,G, and B) has the same EOCF.
In this package, the EOCF is always defined as a map of [0,1] to itself.
The classical EOCF is the "\eqn{\gamma} power law".
That is \eqn{R' = R^\gamma}, where \eqn{R'} is linear R and \eqn{R} is display R.
For Adobe RGB, the EOCF is classical with \eqn{\gamma = 563/256}, which is about 2.2.
For the most popular RGB space - sRGB - the EOCF is a continuous function defined in 2 pieces.
The EOCF is also called the \emph{decoding function}.
In ICC profiles, it is called a \emph{Tone Response Curve} (TRC), or \emph{shaper curve}.

The conversion function from optical RGB (linear) to electrical RGB (non-linear) is called the
\emph{Opto-Electronic Conversion Function} (OECF).
It is also a map of [0,1] to itself, and always the inverse of the EOCF.
The classical OECF is the "\eqn{1/\gamma} power law".
The OECF is also called the \emph{encoding function}.

A good reference for all this is \cite{Poynton}.
}


\section{The API}{
The API is small.  There are functions to
\itemize{
\item convert between XYZ and display RGB
\item convert between linear RGB and display RGB
\item install, uninstall, and query the RGB spaces
}
The XYZ vector is always relative to the XYZ of the white point of the RGB space;
there are no functions for chromatic adaptation, as in \cite{Pascale}.
The XYZ of the white-point always has Y=1.

It is also possible to convert between XYZ and linear RGB in one function call
by overriding the OECF and EOCF by setting \code{gamma=1} as follows \cr
\code{RGBfromXYZ( XYZ, gamma=1 )}\cr
\code{XYZfromRGB( XYZ, gamma=1 )}\cr
In the first case, besides being more efficient, there is also no clamping to the unit cube,
so one can tell how far outside the cube that the linear RGB really is.
}



\section{Note}{
All the RGB spaces are stored in a dictionary that initially contains
\code{'sRGB'}, \code{'AdobeRGB'} and \code{'ProPhotoRGB'}.
See \bold{References} for the corresponding Wikipedia pages on these spaces.
}

\section{Other Packages}{
Package \pkg{colorspace} has similary functionality, but only supports one space - \bold{sRGB}.
\cr
Package \pkg{colorscience} has very flexible conversion functions,
which include \code{gamma} override and chromatic adaptation,
and it has all the RGB spaces in \cite{Lindbloom}.
But is structured so it only allows one space with non-classical EOCF/OECF, namely \bold{sRGB}.
}


\author{
Glenn Davis <gdavis@gluonics.com>
}

\references{
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

Wikipedia. \bold{sRGB}.
\url{http://en.wikipedia.org/wiki/SRGB}.

Wikipedia. \bold{Adobe RGB (1998)}.
\url{http://en.wikipedia.org/wiki/Adobe_RGB_color_space}.


Wikipedia. \bold{ProPhoto RGB}.
\url{http://en.wikipedia.org/wiki/ProPhoto_RGB_color_space}.
}



\keyword{package}

