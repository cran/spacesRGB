\name{Signal RGB Calculation}
\alias{RGBfromXYZ}
\alias{SignalRGBfromLinearRGB}
\title{Convert XYZ or Linear RGB to Signal RGB}

\description{Convert XYZ or Linear RGB to Signal RGB, multiple RGB spaces are available}

\usage{
RGBfromXYZ( XYZ, space='sRGB', which='scene', TF=NULL, maxSignal=1 ) 

SignalRGBfromLinearRGB( RGB, space='sRGB', which='scene', TF=NULL, maxSignal=1 )
}

\arguments{
\item{XYZ}{a numeric Nx3 matrix with CIE XYZ coordinates in the rows, or a vector that can be converted to such a matrix, by row.
The XYZ are assumed to be viewed under the white-point of the given RGB \code{space}.}

\item{RGB}{a numeric Nx3 matrix with linear RGB coordinates in the rows,
or a vector that can be converted to such a matrix, by row.
The RGB are assumed normalized to the unit cube, but are allowed to be outside it.}

\item{space}{the name of an installed RGB space.
The name matching is partial and case-insensitive.}

\item{which}{the input linear RGB - either \code{'scene'} or \code{'display'}.
Usually the OOTF for the space is the identity and so these two are the same.
}

\item{TF}{if not \code{NULL}, \code{TF} overrides the appropriate transfer function of \code{space}.
\code{TF} can be a positive number, or function that takes the interval [0,1] to itself.
If \code{TF=1}, the output RGB values are actually linear, and not clamped (see \bold{Value}).
If \code{TF}  is not \code{NULL} in \code{SignalRGBfromLinearRGB()},
then \code{space} and \code{which} are ignored.
}

\item{maxSignal}{maximum value of non-linear \code{RGB}.
Other common values are 100, 255, 1023, 4095, and 65535.
Even when 1, they are still taken to be non-linear Signal values.}
}

\value{
a \code{data.frame} with N rows and these columns
\item{RGB}{signal RGB.
If \code{TF} is not 1, all input linear RGB values are clamped to the unit cube,
which implies that the signal RGBs are \emph{inside} the cube \code{[0,maxSignal]}^\eqn{3}.
Values are not rounded.}

\item{OutOfGamut}{logical vector, TRUE means the input linear RGB was outside the unit cube before clamping it}

In case of error, the functions return NULL.
}


\references{
Wikipedia. \bold{RGB color space}.
\url{https://en.wikipedia.org/wiki/RGB_color_space}
}


\seealso{
\code{\link{XYZfromRGB}()},
\code{\link{LinearRGBfromSignalRGB}()},
\code{\link{installRGB}()}
}

\examples{
RGBfromXYZ( c(80.310897,90.306510,84.613450, 100,100,100)/100, max=255 )
##       RGB.R    RGB.G    RGB.B     OutOfGamut
##  1 230.1514 249.4148 225.2711          FALSE
##  2 255.0000 249.1150 244.4963           TRUE
}

\keyword{RGB}

