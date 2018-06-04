\name{RGB Space Management}
\alias{installRGB}
\alias{uninstallRGB}
\alias{getRGB}
\alias{summaryRGB}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Manage RGB Spaces}

\description{
Install user-defined RGB Spaces, and query and summarize all installed RGB spaces
}

\usage{
installRGB( space, primaries, white, gamma, peakRGB=1, overwrite=FALSE )
uninstallRGB( space )

getRGB( space )

summaryRGB( verbosity=1 )
}

\arguments{
\item{space}{name of the RGB space to install or uninstall or query.
After the RGB space is installed, 
\code{space} it can be used in the conversion functions,
e.g. \code{\link{RGBfromXYZ}()} and \code{\link{XYZfromRGB}()}.}

\item{primaries}{3x2 matrix with the CIE xy chromaticities of R,G,B in the rows, in that order}

\item{white}{a numeric 2-vector or 3-vector defining the whitepoint of the display.
If \code{white} is a 2-vector, it is interpreted as the CIE xy chromaticity of the whitepoint.
If \code{white} is a 3-vector, it is interpreted as the CIE XYZ of the whitepoint,
and the Y must be 1.
For better numeric compatibility with other software, XYZ is recommended.
}

\item{gamma}{
If \code{gamma} is a positive number, it is taken to be the \emph{decoding gamma},
i.e. from display to linear RGB  (the EOCF).
In all the usual spaces, \code{gamma} \eqn{\ge} 1.
\cr
\code{gamma} can also be a list with 2 functions.
The 1st function is display to linear (the EOCF), and the 2nd is linear to display (the OECF).
Of course, they should be inverses.
Both of them should operate on the interval [0,1] and take
0 to 0 and 1 to 1 (this is checked).
Both functions must also be able to take a matrix as argument
and return a matrix with the same dimensions (this is checked).
\cr
\code{gamma} can also be the string \code{'sRGB'} or \code{'ProPhotoRGB'} which then installs
the appropriate special OECF/EOCF conversion functions.
}

\item{peakRGB}{
A number in the interval (0,1], and the neutral RGB that maps to the XYZ of the white point.
The default \code{peakRGB=1} is the conventional value,
but values less than 1 allow for some "white headroom".
The number is replicated 3 times to get a neutral gray.
Use of this parameter should be very rare.
}

\item{overwrite}{
in \code{installRGB()}, \code{space} is compared with previously installed
RGB space names in case-insensitive way.  
If there is a match, and \code{overwrite} is \code{FALSE},
then the installation fails.
If \code{overwrite} is \code{TRUE}, then the existing space is overwritten.
}

\item{verbosity}{
an integer that controls the return value of \code{summaryRGB()}, see \bold{Value}.
}

}


\value{
\code{installRGB()} and \code{uninstallRGB()} return \code{TRUE} or \code{FALSE}.

If \code{verbosity=0}, \code{summaryRGB()} returns the names of all installed RGB spaces.
If \code{verbosity=1}, it returns a \code{data.frame}
with rows that contain  primaries and gamma information for each space.

\code{getRGB()} returns a list with these items:

\describe{
\item{\code{space}}{the full and original name of the RGB space}
\item{\code{primaries}}{4x2 matrix with the xy chromaticities of the RGB primaries plus white}
\item{\code{whiteXYZ}}{XYZ of the display white point, which maps to \code{peakRGB}}
\item{\code{peakRGB}}{the peak linear RGB, which maps to \code{whiteXYZ}.  It is always a neutral gray.}
\item{\code{RGB2XYZ}}{3x3 matrix taking RGB to XYZ}
\item{\code{XYZ2RGB}}{3x3 matrix taking XYZ to RGB}
\item{\code{EOCF}}{Electro-Optical Conversion Function, taking [0,1] to itself}
\item{\code{OECF}}{Opto-Electronic Conversion Function, taking [0,1] to itself}
}

In case of error, \code{getRGB()} returns NULL.
}



\section{Warning}{
All the RGB spaces are stored in a dictionary that initially contains
\code{'sRGB'}, \code{'AdobeRGB'} and \code{'ProPhotoRGB'}.
The dictionary is initially locked, 
but is unlocked in the function \code{\link{.onAttach}()}.
If package \pkg{spacesRGB} is imported using
\code{\link{requireNamespace}()} or  \code{\link{loadNamespace}()},
then \code{\link{installRGB}()} will fail because \code{\link{.onAttach}()} is never called
(\code{\link{.onLoad}()} was tried but failed).
If this is a problem, use \code{\link{attachNamespace}()} instead.

If \code{installRGB()} is successful, the installed space is only active until the end of the \R session.
To make it persist, please put the function call in an \R script that is executed on startup.
}

\details{
Both  \code{installRGB()} and \code{uninstallRGB()} check for matches with existing names.
The matching is full (not partial) and case-insensitive.
So it is impossible to have 2 spaces that differ only in case.

In \code{getRGB()} the name matching is partial and case-insensitive.
}



\seealso{
\code{\link{RGBfromXYZ}()},
\code{\link{XYZfromRGB}()}
}

\examples{
#  install native RGB space for NEC PA242W display
prim = matrix( c(0.675,0.316, 0.199,0.715, 0.157,0.026), 3, 2, byrow=TRUE )
installRGB( 'PA242W', prim, c(0.95047,1,1.08883), gamma=2 )

#  install a linear version of sRGB  (gamma=1)
prim    = matrix( c(0.64,0.33,  0.30,0.60, 0.15,0.06), 3, 2, byrow=TRUE )
installRGB( 'linear-sRGB', prim,  c(0.95047,1,1.08883),  gamma=1 )    

# make plot comparing three EOCFs
plot( c(0,1), c(0,1), type='n', xlab='Electronic', ylab='Optical', las=1 )
grid(lty=1)
x = seq(0,1,by=1/64)
lines( x, getRGB('sRGB')$EOCF(x), col='black' )
lines( x, getRGB('linear')$EOCF(x), col='red' )
lines( x, getRGB('PA242W')$EOCF(x), col='blue' )
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ RGB }% use one of  RShowDoc("KEYWORDS")

