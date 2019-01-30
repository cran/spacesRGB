% File plot.Rd

\name{plot}
\alias{plot.TransferFunction}

\title{plot a TransferFunction}

\description{
plot a \code{TransferFunction} of dimension 1, 2, or 3.
}

\usage{
\S3method{plot}{TransferFunction}( x, color='red', main=TRUE, add=FALSE, ... )
}

\arguments{
\item{x}{a \code{TransferFunction} object with dimension N = 1, 2, or 3.}

\item{color}{
Any value acceptable as the \code{col} argument to \code{graphics::lines()}.
If N=3 this argument is currently ignored.
}

\item{main}{If \code{main=TRUE} then a main title is generated from the object \code{x}.
If \code{main=FALSE} then no main title is displayed.
And if \code{main} is a character string then that string is used as the main title.
If N=3 this argument is currently ignored.
}

\item{add}{
If \code{add=TRUE} then the lines are added to an existing plot.
If N=3 this argument is currently ignored.
}

\item{...}{other graphical parameters, see \bold{Details}}

}

\value{
\code{TRUE} or \code{FALSE}
}

\details{
If N=1 a conventional plot is drawn using \code{graphics::lines()}.
Commonly used graphical parameters applicable when N=1 are:

\tabular{ll}{
\code{log}  \tab  passed on to \code{\link{plot.default}()}. 
Care must be taken because many transfer functions have 0 in their domains. \cr
}
If N=2 a grid is generated in the domain box, and the image of that grid is plotted using using \code{graphics::lines()}.
\cr
If N=3 a grid is generated in the domain box, and the image of that grid is plotted in 3D
using \code{rgl::lines3d()}.
}



\seealso{
\code{\link[graphics:lines]{graphics::lines}()},
\code{\link[rgl:3dobjects]{rgl::lines3d}()}
}


