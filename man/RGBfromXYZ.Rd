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
The RGB may be outside the corresponding domain box (either scene or display), see \bold{Details}.
}

\item{space}{the name of an installed RGB space.
The name matching is partial and case-insensitive.}

\item{which}{either \code{'scene'} or \code{'display'}.
For \code{RGBfromXYZ()} \code{which} describes the input XYZ.
For \code{SignalRGBfromLinearRGB()} \code{which} describes the input RGB.
% Usually the OOTF for the space is the identity and so these two are the same.
}

\item{TF}{if not \code{NULL}, \code{TF} is a \code{TransferFunction}
that overrides the appropriate transfer function of \code{space}.
\code{TF} can also be a positive number.
If \code{TF=1}, then \code{TF} is set to \code{identity.TF},
so the returned RGB values are actually linear, and they are not clamped to the appropriate domain box
(see \bold{Value}).
If \code{TF!=1} it is used to create either
\code{power.EOTF()} or \code{power.OETF()} as approriate.
If \code{TF} is not \code{NULL} in \code{SignalRGBfromLinearRGB()}, then \code{space} is ignored.
}

\item{maxSignal}{maximum value of non-linear \code{RGB}.
Other common values are 100, 255, 1023, 4095, and 65535.
Even when 1, they are still taken to be non-linear Signal values.}
}

\value{
a \code{data.frame} with N rows and these columns
\item{RGB}{signal RGB.
If \code{TF} is not the identity, all input linear RGB values are clamped to the appropriate domain box,
which implies that the signal RGBs are inside the cube \code{[0,maxSignal]}^\eqn{3}.
Values are not rounded.}

\item{OutOfGamut}{logical vector, TRUE means the input linear RGB was outside the domain box before clamping it.}

In case of error, the functions return NULL.
}

\details{
In \code{RGBfromXYZ()}, the conversion is done in 2 steps:
\itemize{
\item XYZ  \rarrow  linear RGB using the appropriate 3x3 matrix for the given \code{space} and \code{which}
\item linear RGB  \rarrow  signal RGB using \code{SignalRGBfromLinearRGB()} and all other the given arguments
}
}

\references{
Wikipedia. \bold{RGB color space}.
\url{https://en.wikipedia.org/wiki/RGB_color_space}
}


\seealso{
\code{\link{XYZfromRGB}()},
\code{\link{LinearRGBfromSignalRGB}()},
\code{\link{installRGB}()},
\code{\link{identity.TF}}
}

\examples{
RGBfromXYZ( c(80.310897,90.306510,84.613450, 100,100,100)/100, max=255 )
##       RGB.R    RGB.G    RGB.B OutOfGamut
##  1 230.1676 249.4122 225.2472      FALSE
##  2 255.0000 249.1125 244.4704       TRUE
}

\keyword{RGB}

