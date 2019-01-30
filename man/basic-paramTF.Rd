\name{Basic Parameterized TransferFunctions}
\alias{power.OETF}
\alias{power.EOTF}
\alias{power.OOTF}
\alias{affine.TF}       
\alias{BT.1886.EOTF}
\alias{XYZfromRGB.TF}


\title{
Basic Parameterized TransferFunctions
}

\description{
Each of these functions \emph{returns} a \code{TransferFunction} object,
that depends on the argument values passed to it.
The returned object has the parameter values "locked in".
These \code{TransferFunction} objects are a mixture of EOTFs, OETFs, OOTFs, and general-purpose transfer functions.
}


\usage{
power.OETF( gamma )
power.EOTF( gamma )
power.OOTF( gamma )
BT.1886.EOTF( gamma=2.4, Lb=0, Lw=1 )
XYZfromRGB.TF( primaries, white )
affine.TF( y0, y1 )
}

\arguments{
\item{gamma}{the value of \eqn{\gamma}; it must be positive}
\item{Lb}{the black level}
\item{Lw}{the white level}

\item{primaries}{a 3x2 or 4x2 matrix; see \bold{Details}}
\item{white}{a vector of length 1, 2, or 3; see \bold{Details}}

\item{y0}{the number to which 0 maps}
\item{y1}{the number to which 1 maps}
}

\details{
There are 3 valid combinations of \code{primaries} and \code{white}, as given in this table:
\tabular{ccl}{
\code{dim(primaries)}  \tab \code{length(white)} \tab Description \cr
4x2 \tab 1 \tab  \code{primaries} is a 4x2 matrix with CIE xy chromaticities of R,G,B,W in the rows\cr
3x2 \tab 2 \tab  \code{primaries} is a 3x2 matrix with CIE xy chromaticities of R,G,B in the rows\cr
3x2 \tab 3 \tab  \code{primaries} is a 3x2 matrix with CIE xy chromaticities of R,G,B in the rows\cr
}
If \code{length(white)} is 1, then \code{white} is the whitepoint Y.
If \code{length(white)} is 2, then \code{white} is the whitepoint xy (CIE chromaticity); 
the whitepoint Y is taken to be 1.
If \code{length(white)} is 3, \code{white} is the whitepoint XYZ (CIE tristimulus).
\cr
\code{primaries} can also be a plain numeric vector of length 6 or 8,
which is then converted to a 3x2 or 4x2 matrix, by row.
}

\value{
\code{power.OETF()} returns a \code{\link{TransferFunction}} with the classical \eqn{1/\gamma} power law.
\code{power.EOTF()} returns a \code{\link{TransferFunction}} with the classical \eqn{\gamma} power law.
\code{power.OOTF()} is the same as \code{\link{power.EOTF}()},
but having a different name may make the creation of new RGB spaces clearer.
All three of these map [0,1] to [0,1].

\code{BT.1886.EOTF()} returns a \code{\link{TransferFunction}} that maps [0,1] to
[\code{Lb},\code{Lw}], with non-linearity given by \code{gamma}.
The BT.1886 standard has details in Annex 1.

\code{XYZfromRGB.TF()} returns a 3D \code{\link{TransferFunction}} that is linear
and maps RGB=(1,1,1) to the XYZ of white.
The domain is set to the ACES cube  \eqn{[-65504, 65504]^3} and the range is set to the smallest enclosing box.
For the inverse one can use \code{XYZfromRGB.TF()^-1}.

\code{affine.TF()} returns a 1D \code{\link{TransferFunction}} that maps 
0 \rarrow \eqn{y_0} and 1 \rarrow \eqn{y_1} in an affine way.
One must have \eqn{y_0 \ne y_1}, but is is OK to have \eqn{y_0 > y_1}.
No quantities are associated with these values; the function is intended for arbitrary 1D scaling.
}


\references{
BT.1886.
Reference electro-optical transfer function for flat panel displays used in HDTV studio production.
March 2011.
}

\seealso{
\code{\link{TransferFunction}}
}
