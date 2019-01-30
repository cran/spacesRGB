\name{RGB Space Management}
\alias{installRGB}
\alias{uninstallRGB}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{Manage RGB Spaces}

\description{
Install/uninstall RGB spaces in a dictionary. 
The dictionary comes with 8 RGB spaces pre-installed.
To query the dictionary, use \code{\link{getRGB}()} and \code{\link{summaryRGB}()}.
}

\usage{
installRGB( space, scene, display=NULL, OETF=NULL, EOTF=NULL, OOTF=NULL, overwrite=FALSE )
uninstallRGB( space )
}

\arguments{
\item{space}{name of the RGB space to install or uninstall or query.
After the RGB space is installed,
the string \code{space} can be used in the conversion functions,
e.g. \code{\link{RGBfromXYZ}()} and \code{\link{XYZfromRGB}()}.}

\item{scene}{
the specification of the scene primaries and whitepoint.
There are many options here.
The 1st option is a 4x2 matrix with the CIE xy chromaticities of R,G,B,W in the rows, in that order.
The 2nd option is a list with 2 items: the primaries data and the whitepoint data.
These are described in the section \bold{Primaries and Whitepoint Details} below.
If \code{scene} is \code{NULL}, it will duplicate the data from argument \code{display}.
}

\item{display}{
the specification of the display primaries and whitepoint.
The options are the same as for argument \code{scene}.
If this is \code{NULL} (the default), the function will first look at the metadata of the transfer functions.
These built-in transfer functions
- \code{general.OOTF()} and \code{general.PODT()}  - 
already have this metadata assigned.
If the metadata is not found, it will duplicate the data from argument \code{scene}.
}


\item{OETF}{
a \code{\link{TransferFunction}} of dimension 1 or 3.
\cr
\code{OETF} can also be a positive number \eqn{\gamma},
which is then passed to \code{\link{power.OETF}()} to create the \code{\link{TransferFunction}}.
This is the classical \eqn{1/\gamma} power law.
\cr
\code{OETF} can also be \code{NULL}; see section \bold{Transfer Function Details} for valid combinations.
}

\item{EOTF}{
a \code{\link{TransferFunction}} of dimension 1 or 3.
\cr
\code{EOTF} can also be a positive number \eqn{\gamma},
which is then passed to \code{\link{power.EOTF}()} to create the \code{\link{TransferFunction}}.
This is the classical \eqn{\gamma} power law.
\cr
\code{EOTF} can also be one of these strings:
\code{'sRGB'}, \code{'ProPhotoRGB'}, \code{'BT.709'},  \code{'BT.2020'}, or  \code{'240M'},
which then installs the appropriate special EOTF function.
\cr
\code{EOTF} can also be \code{NULL}; see section \bold{Transfer Function Details} for valid combinations.
}

\item{OOTF}{
a \code{\link{TransferFunction}} of dimension 1 or 3.
\cr
\code{OOTF} can also be a positive number \eqn{\gamma},
which is then passed to \code{\link{power.OOTF}()} to create the \code{\link{TransferFunction}}.
This is the classical \eqn{\gamma} power law.
\cr
\code{EOTF} can also be \code{NULL}; see section \bold{Transfer Function Details} for valid combinations.
}

\item{overwrite}{
in \code{installRGB()}, \code{space} is compared with previously installed
RGB space names in case-insensitive way.  
If there is a match, and \code{overwrite} is \code{FALSE},
then the installation fails.
If \code{overwrite} is \code{TRUE}, then the existing space is overwritten.
}

}


\value{
\code{installRGB()} and \code{uninstallRGB()} return \code{TRUE} or \code{FALSE}.
}

\section{Primaries and Whitepoint Details}{
The arguments \code{scene} and \code{display} can be a list with 2 items: 
\code{primaries} and \code{white} in that order.
There are 3 options for this list, as given in this table:
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
The whitepoint is linearly transformed to RGB=(1,1,1).
For better numeric compatibility with standards, xy is recommended.
For better numeric compatibility with \cite{Lindbloom}, XYZ is recommended.
See the \bold{Examples} below.
}


\section{Transfer Function Details}{
The 3 transfer functions - \code{OETF}, \code{EOTF}, \code{OOTF} - can be \code{NULL} (the default) or given.
This yields 8 combinations, but only 6 are valid, as given in this table:
\tabular{cccl}{
\code{OETF} \tab \code{EOTF} \tab \code{OOTF}         \tab Description \cr
 given      \tab   given     \tab   given             \tab INVALID \cr
 given      \tab   given     \tab   \code{OETF*EOTF}  \tab \code{OOTF} is the composition \code{OETF} followed by \code{EOTF} \cr
 given      \tab  \code{OETF^-1*OOTF}      \tab given \tab \code{EOTF} is the composition \code{OETF^-1} followed by \code{OOTF} \cr
 \code{OOTF*EOTF^-1} \tab  given    \tab given        \tab \code{OETF} is the composition \code{OOTF} followed by \code{EOTF^-1} \cr 
 
given \tab  \code{OETF^-1}  \tab  \code{identity.TF}  \tab \code{EOTF} is set to \code{OETF^-1}, and \code{OOTF} is set to the identity \cr 
  
 \code{EOTF^-1}  \tab  given  \tab  \code{identity.TF}  \tab \code{OETF} is set to \code{EOTF^-1}, and \code{OOTF} is set to the identity \cr 

\code{NULL} \tab \code{NULL} \tab  given              \tab INVALID \cr
\code{NULL} \tab \code{NULL} \tab \code{NULL}         \tab all 3 transfer functions are set to \code{identity.TF}.  \cr
}
Think of these 3 functions as forming a triangle.
If all 3 are given, the transfers may be ambiguous, i.e. the triangle may not commute.
If 2 functions are given, the 3rd is computed from those 2.
If only 1 function is given, and it is EOTF or OETF, then it makes sense to make the other one the inverse
of the given one, so that the OOTF is the identity.
If only the OOTF is given, there is no well-defined way to define the other 2.
If none are given, as in the last row,
this might be useful for testing conversion between RGB and XYZ.
}



\section{Warning}{
All the RGB spaces are stored in a dictionary.
If \code{installRGB()} is successful, the installed space is only in
the dictionary until the end of the \R session.
To make it persist, please put the function call in an \R script that is executed
after the package is loaded.
\cr
The dictionary comes with 8 RGB spaces pre-installed.
}

\details{
Both  \code{installRGB()} and \code{uninstallRGB()} check for matches with existing names.
The matching is full (not partial) and case-insensitive.
So it is impossible to have 2 spaces that differ only in case.
}


\references{
Lindbloom, Bruce.
\bold{RGB/XYZ Matrices}.
\url{http://brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html}
}


\seealso{
\code{\link{getRGB}()},
\code{\link{summaryRGB}()}
\code{\link{RGBfromXYZ}()},
\code{\link{XYZfromRGB}()},
\code{\link{TransferFunction}},
\code{\link{power.OETF}()},
\code{\link{power.EOTF}()},
\code{\link{power.OOTF}()}
}

\examples{
#  install native RGB space for NEC PA242W display
prim = matrix( c(0.675,0.316, 0.199,0.715, 0.157,0.026), 3, 2, byrow=TRUE )
installRGB( 'PA242W', scene=NULL, display=list(primaries=prim,white=c(0.95047,1,1.08883)), OETF=2 )

#  install a linear version of sRGB  (OETF=1)
prim    = matrix( c(0.64,0.33,  0.30,0.60, 0.15,0.06), 3, 2, byrow=TRUE )
installRGB( 'linear-sRGB', scene=NULL, display=list(prim,c(0.3127,0.3290)), OETF=1 )

# make plot comparing three EOTFs
plot( getRGB('sRGB')$EOTF, col='black' )
plot( getRGB('linear')$EOTF, col='red', add=TRUE )
plot( getRGB('PA242W')$EOTF, col='blue', add=TRUE )

# Install an RGB space named 'HD+2.4', with encoding from BT.709 and display from BT.1886.
# the OOTF for this space is non-trivial
prim    = matrix( c(0.64,0.33,  0.30,0.60,  0.15,0.06 ), 3, 2, byrow=TRUE )
white   = c( 0.3127, 0.3290 )
installRGB( "HD+2.4", scene=NULL, display=list(prim,white),
                          OETF=(BT.709.EOTF)^-1, EOTF=BT.1886.EOTF(), over=TRUE )

# make plot comparing two OOTFs
plot( getRGB('HD+2.4')$OOTF, col='red')
plot( getRGB('sRGB')$OOTF, col='black', add=TRUE  )
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ RGB }% use one of  RShowDoc("KEYWORDS")

