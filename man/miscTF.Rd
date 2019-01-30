% File miscTF.Rd

\name{miscTF}
\title{Miscellaneous TransferFunction Methods}
\alias{dimension}
\alias{dimension.TransferFunction}
\alias{domain}
\alias{domain.TransferFunction}
\alias{orientation}
\alias{orientation.TransferFunction}



\description{
Miscellaneous \code{TransferFunction} methods
}


\usage{
\S3method{dimension}{TransferFunction}( TF )

\S3method{domain}{TransferFunction}( TF )

\S3method{orientation}{TransferFunction}( TF )
}

\arguments{
\item{TF}{a \code{TransferFunction} object}
}


\value{
\code{dimension()} returns a positive integer - the dimension of the domain and range of \code{TF}.
If \code{TF} is a universal identity, it returns \code{NA}.

\code{domain()} returns a 2xN matrix with the domain box of \code{TF},
where N is \code{dimension(TF)}.
If \code{TF} is a universal identity, it returns \code{NA}.

\code{orientation()} returns a real number.
If the value is positive it means that \code{TF} preserves orientation.
If the value is negative it means that \code{TF} reverses orientation.
When \code{dimension(TF)}=1, this simply corresponds to
the function being monotone increasing or decreasing, respectively.
In case of ERROR, the function returns \code{NA}.
}

\examples{
TF = affine.TF( 1, 108 )

dimension(TF)					# 1

orientation(TF)					# 107

orientation( affine.TF( 100, 1 ) )		# -99

domain(TF)
##      AU
##  min  0
##  max  1
}



\seealso{
\code{\link{TransferFunction}}
\code{\link{identity.TF}}
}

