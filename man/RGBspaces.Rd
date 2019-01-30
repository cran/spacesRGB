\name{RGB Space Query}
\alias{getRGB}
\alias{getWhiteXYZ}
\alias{summaryRGB}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Query RGB Spaces}

\description{
Query and summarize the installed RGB spaces.
The RGB spaces are stored in a dictionary, which comes with 8 RGB spaces pre-installed.
}

\usage{
summaryRGB( verbosity=1 )

getRGB( space )
getWhiteXYZ( space, which='scene' )
}

\arguments{
\item{space}{name of the RGB space to query.
The name matching is partial and case-insensitive.
}

\item{verbosity}{
an integer that controls the return value of \code{summaryRGB()}, see \bold{Value}.
}

\item{which}{the source of the whitepoint, either \code{'scene'} or \code{'display'}.
Matching is partial and case-insensitive.
}
}

\details{
The function \code{getWhiteXYZ()} is provided because some applications only need the whitepoint
for chromatic adaptation purposes, and this function is faster than \code{getRGB()}.
}


\value{
\code{summaryRGB()}, with the default \code{verbosity=1}, returns a \code{data.frame}
with a row for each RGB space.
The row contains  primary, whitepoint, and transfer function information for each space.
The primary/whitepoint data is for both scene and display; all the data is numerical and the columns are labeled.
There are 22 columns so the display is very wide.
\cr
The transfer function data is a very short string.
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
\cr
If the \code{TransferFunction} has dimension 1, but the domain and range are not the interval [0,1],
the string is \code{'1D'}.
If the \code{TransferFunction} has dimension 3, the string is \code{'3D'}.
\cr
If \code{verbosity=0}, \code{summaryRGB()} returns the names of all the RGB spaces.


\code{getRGB()} returns a list with these items:

\tabular{ll}{
\code{space}   \tab   the full and original name of the RGB space  \cr
\code{scene} 	 \tab   a list with items \code{primaries}, \code{whiteXYZ}, \code{RGB2XYZ}, and \code{XYZ2RGB}  \cr
\code{display} \tab   a list with items \code{primaries}, \code{whiteXYZ}, \code{RGB2XYZ}, and \code{XYZ2RGB}  \cr
\code{EOTF}    \tab   Electro-Optical Transfer Function  \cr
\code{OETF}    \tab   Opto-Electronic Transfer Function  \cr
\code{OOTF}    \tab   Opto-Optical Transfer Function, and numerically equal to \code{OETF*EOTF}  \cr
}

The items in the lists \code{scene} and \code{display} are
\tabular{ll}{
\code{primaries} \tab   4x2 matrix with the xy chromaticities of the RGB primaries and white \cr
\code{whiteXYZ}  \tab   XYZ of the display white point, which maps to RGB=(1,1,1) \cr
\code{RGB2XYZ}   \tab   3x3 matrix taking RGB to XYZ \cr
\code{XYZ2RGB}   \tab   3x3 matrix taking XYZ to RGB \cr
}
All transfer functions are actual \code{TransferFunctions} objects, 
and not the numerical exponent or character string name.
They are suitable for plotting with \code{\link{plot.TransferFunction}()}; see the \bold{Examples}.
In case of error, \code{getRGB()} returns NULL.


\code{getWhiteXYZ()} 
returns a numeric 3-vector with the XYZ of the whitepoint of the scene or the display.
In case of error it returns \code{NULL}.
}




\references{
Lindbloom, Bruce.
\bold{RGB/XYZ Matrices}.
\url{http://brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html}
}


\seealso{
\code{\link{installRGB}()},
\code{\link{plot.TransferFunction}()}
}

\examples{
# make plot comparing three EOTFs
plot( getRGB('sRGB')$EOTF, col='black' )
plot( getRGB('BT.709')$EOTF, col='blue', add=TRUE )
plot( getRGB('ProPhotoRGB')$EOTF, col='red', add=TRUE )
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ RGB }% use one of  RShowDoc("KEYWORDS")

