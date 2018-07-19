% File plotPatchesRGB.Rd

\name{plotPatchesRGB}
\title{Plot Patches defined by RGB}
\alias{plotPatchesRGB}


\description{
RGB patches are a very common way of comparing color renderings.
This function draws rectangular patches,
and can also draw triangles formed by omitting one vertex from the rectangle.
}

\usage{
plotPatchesRGB( obj, space='sRGB', which='signal', maxColorValue=1,
                     background='gray50', shape='full', add=FALSE, labels=FALSE, ... )
}


\arguments{

\item{obj}{
an Nx3 matrix of RGBs for N patches, preferably with assigned \code{rownames}.
\code{obj} can also be a \code{data.frame} containing a unique matrix column whose name
starts with the string \code{'RGB'}.
If \code{obj} has columns \code{LEFT,TOP,WIDTH,HEIGHT} then these are used to place the patches,
with the Y coordinate increasing going \emph{down} the page.
If \code{obj} has columns \code{LEFT,BOTTOM,WIDTH,HEIGHT} then these are used to place the patches,
with the Y coordinate increasing going \emph{up} the page.
If there are no columns defining the location and size of the patches,
then defaults are supplied, see \bold{Details}.
}

\item{space}{the name of an installed RGB space.
When the input RGB is linear, a transfer function of this RGB space is used
to convert linear RGB to signal RGB, see \bold{Details}.}

\item{which}{the meaning of the RGB values in \code{obj}.
Valid values are \code{'signal'}, \code{'scene'}, and \code{'display'}.
See the Figure on page 2.
Partial matching is used.
For the RGB processing, see \bold{Details}.
}

\item{maxColorValue}{a positive number used for input RGB scaling, see \bold{Details}}

\item{background}{
the color for the background behind all the patches.
If it is a character string, it is passed directly to \code{par()} as parameter \code{bg}.
If it is a numeric vector of length 3, it is processed just like the input RGB in \code{obj},
see  \bold{Details}.
If it is a number, it is interpreted as graylevel,
replicated to length 3, and treated as in the previous sentence.}
  

\item{shape}{
If \code{shape='full'} (the default) then the full rectangle is drawn.
\cr
If \code{shape='half'} then the rectangle is shrunk to 1/2 size, and with the same center.
\cr
If \code{shape} is one of \code{'left'}, \code{'right'}, \code{'bottom'}, or \code{'top'}
then only a  half-rectangle is drawn, and keeping the specified side.
\cr
If \code{shape} is one of \code{'topleft'}, \code{'topright'}, \code{'bottomleft'}, or \code{'bottomright'},
then only a triangular half of the rectangle is drawn, 
and keeping the specified vertex.
\cr
If \code{shape='hhex'} then a hexagon is drawn inscribed in the rectangle with 2 horizontal
opposite sides (in contact with the rectangle sides).
And if the aspect ratio of the rectangle is \eqn{ 2:\sqrt{3} } the hexagon is regular.
If \code{shape='vhex'} then the inscribed rectangle has 2 vertical opposite sides.
}
               
\item{add}{if \code{TRUE} then the patches are added to an existing plot}

\item{labels}{controls how the patches are labeled, using \code{rownames(obj)},
or 1:N if \code{rownames(obj)} is \code{NULL}.
The function used is \code{graphics::\link{text}()}.
If \code{labels=FALSE} then no labels are plotted.
If \code{labels=TRUE} then labels are plotted in the center of the patch when there are
columns defining the location and size of the patches, and to the right of the patch otherwise.
\cr
\code{labels} can also be a character string defining the location where the labels are drawn.
It can be the side of the patches, i.e. \code{left}, \code{right}, \code{top}, or \code{bottom},
or the corner of the patches, i.e. \code{bottomleft}, \code{bottomright}, \code{topleft}, or \code{topright}.
}

\item{...}{additional arguments passed to \code{graphics::\link{text}()}.
For example: \code{adj}, \code{cex}, etc.
}


}


\value{ \code{TRUE} if successful, and \code{FALSE} otherwise}


\details{
If \code{which='signal'} then the input RGBs are converted to hex codes using
\code{\link{rgb}()} using the \code{maxColorValue} argument,
and the \code{space} argument is ignored.
\cr
If \code{which='scene'} or  \code{which='display'} then the input linear RGBs are normalized
by division by \code{maxColorValue}, and then converted to signal RGB using
\code{\link{SignalRGBfromLinearRGB}()} with the \code{space} argument.
The signal RGB is then converted to hex codes using \code{\link{rgb}()}.

If \code{obj} is a matrix, or a data.frame without columns \code{LEFT,TOP,WIDTH,HEIGHT},
then the patches are drawn vertically stacked and abutting from top to bottom.
}

\seealso{
\code{\link{SignalRGBfromLinearRGB}()},
\code{\link{installRGB}()},
\code{\link{rgb}()}
}

\examples{
set.seed(0)
RGB = round( 255 * matrix( runif(6*3), 6, 3 ) )
plotPatchesRGB( RGB, max=255 )
}


\keyword{RGB}
