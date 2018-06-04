\name{Display RGB Calculation}
\alias{RGBfromXYZ}
\alias{DisplayRGBfromLinearRGB}
\title{Convert XYZ or Linear RGB to Display RGB}

\description{Convert XYZ or Linear RGB to Display RGB, multiple RGB spaces are available}

\usage{
RGBfromXYZ( XYZ, space='sRGB', gamma=NULL, maxValue=1 ) 

DisplayRGBfromLinearRGB( RGB, space='sRGB', gamma=NULL, maxValue=1 )
}

\arguments{
\item{XYZ}{a numeric Nx3 matrix with CIE XYZ coordinates in the rows, or a vector that can be converted to such a matrix, by row.
The XYZ are assumed to be viewed under the white-point of the given RGB \code{space}
with Y=1.}

\item{RGB}{a numeric Nx3 matrix with linear RGB coordinates in the rows,
or a vector that can be converted to such a matrix, by row.
The RGB are assumed normalized to the unit cube, but are allowed to be outside it.}

\item{space}{the name of an installed RGB space.
The name matching is partial and case-insensitive.
Spaces \code{'sRGB'} and \code{'AdobeRGB'} and \code{'ProPhotoRGB'}  are pre-installed,
and others can be installed with \code{\link{installRGB}()}.}

\item{gamma}{if not NULL, \code{gamma} overrides the gamma of \code{space}.
\code{gamma} can be a positive number, or function that takes the interval [0,1] to itself.
Such a function is an Opto-Electric Conversion Function, or OECF.
If \code{gamma=1}, the output RGB values are actually linear, and not clamped (see \bold{Value}). }

\item{maxValue}{maximum value of non-linear \code{RGB}.
Other common values are 100, 255, 1023, 4095, and 65535.
Even when 1, they are still taken to be non-linear display values.}
}

\value{
a \code{data.frame} with N rows and these columns
\item{RGB}{non-linear RGB coordinates intended for display (display-referred). 
If \code{gamma} is not 1, all linear RGB values are clamped to the unit cube,
which implies that the display RGBs are inside the cube \eqn{[0,maxValue]^3}.
Values are not rounded.}

\item{OutOfGamut}{logical vector, TRUE means the computed linear RGB was outside the unit cube before clamping it}

In case of error, the functions return NULL.
}


\references{
Wikipedia. \bold{RGB color space}.
\url{https://en.wikipedia.org/wiki/RGB_color_space}
}

\author{Glenn Davis}

\seealso{
\code{\link{XYZfromRGB}()},
\code{\link{LinearRGBfromDisplayRGB}()},
\code{\link{installRGB}()}
}

\examples{
RGBfromXYZ( c(80.310897,90.306510,84.613450, 100,100,100)/100, max=255 )
##       RGB.R    RGB.G    RGB.B     OutOfGamut
##  1 230.1514 249.4148 225.2711          FALSE
##  2 255.0000 249.1150 244.4963           TRUE
}

\keyword{RGB}

