\name{Linear RGB, XYZ, and Lab Calculation}
\alias{XYZfromRGB}
\alias{LabfromRGB}
\alias{LinearRGBfromSignalRGB}
\title{Convert Signal RGB coordinates to Linear RGB, XYZ, or Lab}

\description{Convert signal RGB coordinates to linear RGB, XYZ, or Lab}

\usage{
LinearRGBfromSignalRGB( RGB, space='sRGB', which='scene', TF=NULL, maxSignal=1 ) 

XYZfromRGB( RGB, space='sRGB', which='scene', TF=NULL, maxSignal=1 ) 

LabfromRGB( RGB, space='sRGB', which='scene', TF=NULL, maxSignal=1 ) 
}


\arguments{
\item{RGB}{ a numeric Nx3 matrix with non-linear signal RGB coordinates in the rows, or a vector that can be converted to such a matrix, by row.
These should be in the appropriate cube \code{[0,maxSignal]}\eqn{^3};
but see \code{OutOfGamut} in \bold{Value} below.
}

\item{space}{the name of an installed RGB space.
The name matching is partial and case-insensitive.}

\item{which}{either \code{'scene'} or \code{'display'}.
For \code{LinearRGBfromSignalRGB()}, \code{which} describes the output RGB.
For \code{XYZfromRGB()}, \code{which} describes the output XYZ.
For \code{LabfromRGB()}, \code{which} describes the output Lab.
Usually the OOTF for the space is the identity and so \code{'scene'} and \code{'display'} are the same.
}

\item{TF}{if not \code{NULL}, \code{TF} is a \code{TransferFunction}
that overrides the appropriate transfer function of \code{space}.
\cr
\code{TF} can also be a positive number.
If \code{TF=1}, then \code{TF} is set to \code{identity.TF},
so the returned RGB values are not clamped (see \bold{Value}).
If \code{TF!=1} it is used to create either
\code{power.EOTF()} or \code{power.OETF()} as approriate.
If \code{TF} is not \code{NULL} in \code{LinearRGBfromSignalRGB()}, then \code{space} is ignored.
}


\item{maxSignal}{maximum value of the input signal \code{RGB}.
It is just a simple scaling factor;
other common values are 100, 255, 1023, 4095, and 65535.
Even when \code{maxSignal}=1, the RGBs are still taken to be non-linear signal values.}
}

\details{
In \code{XYZfromRGB()}, the conversion is done in 2 steps:
\itemize{
\item signal RGB  \rarrow  linear RGB using \code{LinearRGBfromSignalRGB()} and all other the given arguments
\item linear RGB  \rarrow  XYZ using the 3x3 matrix for the given \code{space} and \code{which}
}
}


\value{
\code{XYZfromRGB()} returns a \code{data.frame} with N rows and these columns:

\item{OutOfGamut}{a logical vector. TRUE means the input signal RGB is outside the cube \code{[0,maxSignal]}\eqn{^3}.
If \code{TF} is not \code{identity.TF}, 
then the input signal RGBs are clamped to this cube before further calculations.}

\item{XYZ}{the calculated XYZ vectors, using the appropriate 3x3 matrix for the RGB space.
These are for viewing under the white point (either \code{'scene'} or \code{'display'})
of the given RGB space.}



\code{LabfromRGB()} returns a \code{data.frame} with N rows and these columns:

\item{OutOfGamut}{a logical vector. TRUE means the input signal RGB is outside the cube \code{[0,maxSignal]}\eqn{^3}.
If \code{TF} is not \code{identity.TF}, 
then the input signal RGBs are clamped to this cube before further calculations.}

\item{Lab}{the calculated Lab vectors.
First XYZ is calculated using \code{XYZfromRGB()}.
Then XYZ \rarrow Lab
is calculated using
\code{\link[spacesXYZ:fromXYZ]{spacesXYZ::LabfromXYZ}()},
with the white point
(either \code{'scene'} or \code{'display'}) of the given RGB space.
If R=G=B exactly, then a=b=0 exactly.}



\code{LinearRGBfromSignalRGB()} returns a \code{data.frame} with N rows and these columns:

\item{OutOfGamut}{a logical vector. TRUE means the input signal RGB is outside the cube \code{[0,maxSignal]}\eqn{^3}.
If \code{TF} is not \code{identity.TF}, then the signal RGBs are clamped to this cube before linearizing.}

\item{RGB}{the calculated linear RGB vectors, either scene linear or display linear, as specified by \code{which}.}

In case of error, all functions return NULL.
}

\references{
Wikipedia. \bold{RGB color space}.
\url{https://en.wikipedia.org/wiki/RGB_color_space}
}


\seealso{
\code{\link{RGBfromXYZ}()},
\code{\link{SignalRGBfromLinearRGB}()},
\code{\link{installRGB}()},
\code{\link[spacesXYZ:fromXYZ]{spacesXYZ::LabfromXYZ}()}
}

\examples{
XYZfromRGB( c(128,200,255, 0,0,0,  255,255,255), max=255 )$XYZ * 100
##              X         Y        Z
##  [1,] 47.60334  53.11601 102.3549
##  [2,]  0.00000   0.00000   0.0000
##  [3,] 95.04559 100.00000 108.9058


XYZfromRGB( c(128,200,255, 0,0,0,  255,255,255), space='Adobe', max=255 )$XYZ   * 100
##              X         Y        Z
##  [1,] 42.36398  50.82876 103.8704
##  [2,]  0.00000   0.00000   0.0000
##  [3,] 95.04559 100.00000 108.9058
}


\keyword{RGB}

