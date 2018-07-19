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
installRGB( space, primaries, white, OETF, EOTF=NULL, overwrite=FALSE )
uninstallRGB( space )

getRGB( space, full=TRUE )

summaryRGB( verbosity=1 )
}

\arguments{
\item{space}{name of the RGB space to install or uninstall or query.
After the RGB space is installed,
the string \code{space} can be used in the conversion functions,
e.g. \code{\link{RGBfromXYZ}()} and \code{\link{XYZfromRGB}()}.}

\item{primaries}{3x2 matrix with the CIE xy chromaticities of R,G,B in the rows, in that order}

\item{white}{a numeric 2-vector or 3-vector defining the whitepoint of the display.
If \code{white} is a 2-vector, it is interpreted as the CIE xy chromaticity of the whitepoint.
If \code{white} is a 3-vector, it is interpreted as the CIE XYZ of the whitepoint.
All numbers must be positive, but it is not necessary for Y to be 1.
The whitepoint is linearly transformed to RGB=(1,1,1).
For better numeric compatibility with standards, xy is recommended.
For better numeric compatibility with \cite{Lindbloom}, XYZ is recommended.
}

\item{OETF}{
If \code{OETF} is a positive number \eqn{\gamma},
the Opto-Electronic Transfer function is the classical \eqn{1/\gamma} power law.
\cr
\code{OETF} can also be a list with 2 transfer functions.
The 1st function is the OETF, and the 2nd function is the inverse OETF\eqn{^{-1}}.
Both of them should be monotone on the interval [0,1] and take
0 to 0 and 1 to 1 (this is checked).
They must be inverses of each other (this is checked).
Both functions must also be able to take a matrix as argument
and return a matrix with the same dimensions (this is checked).
\cr
\code{OETF} can also be a single transfer function
and then a suitable inverse is computed using \code{\link{splinefun}()}.
\cr
\code{OETF} can also be one of these strings:
\code{'sRGB'}, \code{'ProPhotoRGB'}, \code{'BT.709'},  \code{'BT.2020'}, or  \code{'240M'},
which then installs
the appropriate special OETF/OETF\eqn{^{-1}} conversion functions.
}

\item{EOTF}{
If \code{EOTF} is \code{NULL} (the default),
then an appropriate EOTF and EOTF\eqn{^{-1}} are assigned so that the composition
OOTF = EOTF\eqn{\circ}OETF is the identity.
\cr
\code{EOTF} can also be one the values availabe for \code{OETF}
(positive number, list of 2 functions, single function, or character string)
which works just like \code{OETF}.
}

\item{overwrite}{
in \code{installRGB()}, \code{space} is compared with previously installed
RGB space names in case-insensitive way.  
If there is a match, and \code{overwrite} is \code{FALSE},
then the installation fails.
If \code{overwrite} is \code{TRUE}, then the existing space is overwritten.
}

\item{full}{
a logical that controls the return value of \code{getRGB()}, see \bold{Value}.
}

\item{verbosity}{
an integer that controls the return value of \code{summaryRGB()}, see \bold{Value}.
}

}


\value{
\code{installRGB()} and \code{uninstallRGB()} return \code{TRUE} or \code{FALSE}.

\code{summaryRGB()}, with the default \code{verbosity=0}, returns the names of all installed RGB spaces.
If \code{verbosity=1}, it returns a \code{data.frame}
with rows that contain  primary, whitepoint, and transfer function information for each space.
\cr
If the OETF is classical (pure 1/\eqn{\gamma} power law), the string is 1/\eqn{\gamma}.
If the OETF is not classical, the string is 1/~\eqn{\gamma},
where \eqn{\gamma} is the best-fit (or approximate or effective) \eqn{\gamma} to the OETF in the \eqn{L^1}-norm.
\cr
Similarly, if the EOTF is classical (pure \eqn{\gamma} power law) the string is \eqn{\gamma},
and if the EOTF is not classical the string is ~\eqn{\gamma}.
\cr
The OOTF is the quotient (to 2 decimal places) of the gammas of EOTF and OETF
(either true gamma or best-fit gamma).
If either gamma is best-fit then the string is preceede by a \code{'~'},
which means \emph{effective}.


\code{getRGB()}, with the default \code{full=TRUE}, returns a list with these items:

\describe{
\item{\code{space}}{the full and original name of the RGB space}
\item{\code{primaries}}{4x2 matrix with the xy chromaticities of the RGB primaries plus white}
\item{\code{whiteXYZ}}{XYZ of the display white point, which maps to RGB=(1,1,1)}
\item{\code{RGB2XYZ}}{3x3 matrix taking RGB to XYZ}
\item{\code{XYZ2RGB}}{3x3 matrix taking XYZ to RGB}
\item{\code{EOTF}}{Electro-Optical Transfer Function, taking [0,1] to itself}
\item{\code{OETF}}{Opto-Electronic Transfer Function, taking [0,1] to itself}
\item{\code{OOTF}}{Opto-Optical Transfer Function, taking [0,1] to itself,
and equal to EOTF\eqn{\circ}OETF}
}
All transfer functions are actual functions, and not the numerical exponent.
They are suitable for plotting.
If  \code{full=FALSE}, the transfer functions are omitted for extra speed.
In case of error, \code{getRGB()} returns NULL.
}



\section{Warning}{
All the RGB spaces are stored in a dictionary.
When package \pkg{spacesRGB} is loaded the dictionary is locked, 
but is unlocked in the function \code{\link{.onAttach}()}.
If the package is imported using
\code{\link{requireNamespace}()} or  \code{\link{loadNamespace}()},
then \code{\link{installRGB}()} will fail because \code{\link{.onAttach}()} is never called
(\code{\link{.onLoad}()} was tried but it did not work).
If this is a problem, use \code{\link{attachNamespace}()} instead.

If \code{installRGB()} is successful, the installed space is only in
the dictionary until the end of the \R session.
To make it persist, please put the function call in an \R script that is executed
after the package is attached.
}

\details{
Both  \code{installRGB()} and \code{uninstallRGB()} check for matches with existing names.
The matching is full (not partial) and case-insensitive.
So it is impossible to have 2 spaces that differ only in case.

In \code{getRGB()} the name matching is partial and case-insensitive.
}


\references{
Lindbloom, Bruce.
\bold{RGB/XYZ Matrices}.
\url{http://brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html}
}


\seealso{
\code{\link{RGBfromXYZ}()},
\code{\link{XYZfromRGB}()}
}

\examples{
#  install native RGB space for NEC PA242W display
prim = matrix( c(0.675,0.316, 0.199,0.715, 0.157,0.026), 3, 2, byrow=TRUE )
installRGB( 'PA242W', prim, c(0.95047,1,1.08883), OETF=2 )

#  install a linear version of sRGB  (OETF=1)
prim    = matrix( c(0.64,0.33,  0.30,0.60, 0.15,0.06), 3, 2, byrow=TRUE )
installRGB( 'linear-sRGB', prim,  c(0.3127,0.3290),  OETF=1 )    

# make plot comparing three EOTFs
plot( c(0,1), c(0,1), type='n', xlab='Electronic', ylab='Optical', las=1 )
grid(lty=1)
x = seq(0,1,by=1/64)
lines( x, getRGB('sRGB')$EOTF(x), col='black' )
lines( x, getRGB('linear')$EOTF(x), col='red' )
lines( x, getRGB('PA242W')$EOTF(x), col='blue' )

# Install an RGB space named 'HD+2.4', with encoding from BT.709 and display from BT.1886.
# the OOTF for this space is non-trivial
prim    = matrix( c(0.64,0.33,  0.30,0.60,  0.15,0.06 ), 3, 2, byrow=TRUE )
white   = c( 0.3127, 0.3290 )
installRGB( "HD+2.4", prim, white, OETF='BT.709', EOTF=2.4, over=TRUE )

# make plot comparing two OOTFs
plot( c(0,1), c(0,1), type='n', xlab='Optical in', ylab='Optical out', las=1 )
grid(lty=1)
lines( x, getRGB('sRGB')$OOTF(x), col='black' )
lines( x, getRGB('HD+2.4')$OOTF(x), col='red' )
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ RGB }% use one of  RShowDoc("KEYWORDS")

