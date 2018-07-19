\name{Linear RGB and XYZ Calculation}
\alias{XYZfromRGB}
\alias{LinearRGBfromSignalRGB}
\title{Convert Signal RGB coordinates to XYZ, or Linear RGB}

\description{Convert signal RGB coordinates to XYZ, or to linear RGB}

\usage{
XYZfromRGB( RGB, space='sRGB', which='scene', TF=NULL, maxSignal=1 ) 

LinearRGBfromSignalRGB( RGB, space='sRGB', which='scene', TF=NULL, maxSignal=1 ) 
}


\arguments{
\item{RGB}{ a numeric Nx3 matrix with non-linear signal RGB coordinates in the rows, or a vector that can be converted to such a matrix, by row.
These should be in the appropriate cube \code{[0,maxSignal]}\eqn{^3}.
}

\item{space}{the name of an installed RGB space.
The name matching is partial and case-insensitive.}

\item{which}{the output XYZ or linear RGB - either \code{'scene'} or \code{'display'}.
Usually the OOTF for the space is the identity and so these two are the same.
}

\item{TF}{if not NULL, \code{TF} overrides the appropriate transfer function of \code{space}.
\code{TF} can be a positive number \eqn{\gamma}, or a transfer function that takes the interval [0,1] to itself.
If \code{TF=1}, the input RGB values taken to be linear and not clamped (see \bold{Value}). 
If \code{TF} is not \code{NULL} in \code{LinearRGBfromSignalRGB()},
then \code{space} and \code{which} are ignored.
}

\item{maxSignal}{maximum value of the input signal \code{RGB}.
Other common values are 100, 255, 1023, 4095, and 65535.
Even when 1, they are still taken to be non-linear signal values.}
}

\value{
\code{XYZfromRGB()} returns a \code{data.frame} with N rows and these columns:
\item{XYZ}{the calculated XYZ vectors.
These are for viewing under the white point of the given RGB space.}

\item{OutOfGamut}{a logical vector. TRUE means the input signal RGB is outside the cube \code{[0,maxSignal]}\eqn{^3}.
If \code{TF} is not 1, then the input signal RGBs are clamped to the cube before further calculations.}


\code{LinearRGBfromSignalRGB()} returns a \code{data.frame} with N rows and these columns:
\item{RGB}{the calculated linear RGB vectors. These are normalized to the unit cube.}

\item{OutOfGamut}{a logical vector. TRUE means the input signal RGB is outside the cube \code{[0,maxSignal]}\eqn{^3}.
If \code{TF} is not 1, then the signal RGBs are clamped to the cube before linearizing.}

In case of error, both functions return NULL.
}

\references{
Wikipedia. \bold{RGB color space}.
\url{https://en.wikipedia.org/wiki/RGB_color_space}
}


\seealso{
\code{\link{RGBfromXYZ}()},
\code{\link{SignalRGBfromLinearRGB}()},
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

