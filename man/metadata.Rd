% File metadata.Rd

\name{metadata}
\title{metadata of a TransferFunction object}
\alias{metadata}
\alias{metadata<-}
\alias{metadata.TransferFunction}
\alias{metadata<-.TransferFunction}

\description{
Retrieve or set the metadata of a \bold{TransferFunction} object.
The user is free to set this as he/she wishes.
}

\usage{
\S3method{metadata}{TransferFunction}( x, ...)

\S3method{metadata}{TransferFunction}( x, add=FALSE ) <- value
}

\arguments{
  \item{x}{a \code{TransferFunction} \R object}
  
  \item{...}{optional names of metadata to return}
  
  \item{value}{a named \code{list}.  
  If \code{add} is \code{FALSE}, \code{value} replaces any existing \code{metadata}.   
  If \code{add} is \code{TRUE},  \code{value} is appended to the existing list of metadata.
  If a name already exists, its value is updated using \code{\link{modifyList}()}.
  Unnamed items in  \code{value} are ignored.}
  
  \item{add}{if \code{add=FALSE}, any existing metadata is discarded.
  If \code{add=TRUE} then existing metadata is preserved, using \code{modifyList()}.
}
}

\details{
The metadata list is stored as \code{attr(x,'metadata')}.
After construction this list is empty.
}

\value{
\code{metadata(x)} with no additional arguments returns the complete named list of metadata.
If arguments are present, then only those \code{metadata} items are returned.

}



\seealso{
\code{\link{modifyList}}
}

\examples{

\dontrun{
# get list of *all* metadata
metadata(TF)

# get just the number 'gamma'
metadata( TF, 'gamma' )

# alternative method to get just the number 'gamma'
metadata( TF )$gamma

# set the 'date'
metadata( TF ) = list( date="2016-04-01" )
}
}

\keyword{TransferFunction}
