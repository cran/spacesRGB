\name{Linear RGB and XYZ Calculation}
\alias{XYZfromRGB}
\alias{LinearRGBfromDisplayRGB}
\title{Convert Display RGB coordinates to XYZ, or Linear RGB}

\description{Convert Display RGB coordinates to XYZ, or linear RGB}

\usage{
XYZfromRGB( RGB, space='sRGB', gamma=NULL, maxValue=1 ) 

LinearRGBfromDisplayRGB( RGB, space='sRGB', gamma=NULL, maxValue=1 ) 
}


\arguments{
\item{RGB}{ a numeric Nx3 matrix with non-linear display RGB coordinates in the rows, or a vector that can be converted to such a matrix, by row.
These should be in the appropriate cube \eqn{[0,maxValue]^3}.
}

\item{space}{the name of an installed RGB space.
The name matching is partial and case-insensitive.
Spaces \code{'sRGB'} and \code{'AdobeRGB'} and \code{'ProPhotoRGB'}  are pre-installed,
and others can be installed with \code{\link{installRGB}()}.}

\item{gamma}{if not NULL, \code{gamma} overrides the gamma of \code{space}.
\code{gamma} can be a positive number, or function that takes the interval [0,1] to itself.
Such a function is an Electro-Optical Conversion Function, or EOCF.
If \code{gamma=1}, the input RGB values taken to be linear. }

\item{maxValue}{maximum value of non-linear \code{RGB}.
Other common values are 100, 255, 1023, 4095, and 65535.
Even when 1, they are still taken to be non-linear display values.}
}

\value{
\code{XYZfromRGB()} returns a \code{data.frame} with N rows and these columns:
\item{XYZ}{the calculated XYZ vectors.
These are for viewing under the white point (with Y=1) of the given RGB space.}

\item{OutOfGamut}{a logical vector. TRUE means the input display RGB is outside the cube \eqn{[0,maxValue]^3}.
If \code{gamma} is not 1, then the display RGBs are clamped to the cube before further calculations.}


\code{LinearRGBfromDisplayRGB()} returns a \code{data.frame} with N rows and these columns:
\item{RGB}{the calculated linear RGB vectors. These are normalized to the unit cube.}

\item{OutOfGamut}{a logical vector. TRUE means the input display RGB is outside the cube \eqn{[0,maxValue]^3}.
If \code{gamma} is not 1, then the display RGBs are clamped to the cube before linearizing.}

In case of error, the functions return NULL.
}

\references{
Wikipedia. \bold{RGB color space}.
\url{https://en.wikipedia.org/wiki/RGB_color_space}
}


\seealso{
\code{\link{RGBfromXYZ}()},
\code{\link{DisplayRGBfromLinearRGB}()},
\code{\link{installRGB}()}
}

\examples{
XYZfromRGB( c(128,200,255, 0,0,0,  255,255,255), max=255 )$XYZ * 100
##              X         Y       Z
##  [1,] 47.59995  53.11406 102.332
##  [2,]  0.00000   0.00000   0.000
##  [3,] 95.04700 100.00000 108.883


XYZfromRGB( c(128,200,255, 0,0,0,  255,255,255), space='Adobe', max=255 )$XYZ   * 100
##              X         Y        Z
##  [1,] 42.36074  50.82687 103.8475
##  [2,]  0.00000   0.00000   0.0000
##  [3,] 95.04700 100.00000 108.8830
}


\keyword{RGB}

