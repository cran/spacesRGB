% File inverse.Rd

\name{inverse}
\title{The inverse of a TransferFunction Object}
\alias{inverse}
\alias{inverse.TransferFunction}
\alias{^.TransferFunction}
\alias{is.invertible}
\alias{is.invertible.TransferFunction}



\description{
The function \code{inverse()} returns a \code{TransferFunction} that is the inverse
of the argument (if the argument is invertible).

\code{is.invertible()} tests whether a \bold{TransferFunction} object has an inverse.\cr
}


\usage{
\S3method{inverse}{TransferFunction}( TF )

\S3method{^}{TransferFunction}( TF, n )

\S3method{is.invertible}{TransferFunction}( TF )
}

\arguments{
\item{TF}{a \code{TransferFunction} object}

\item{n}{an integer exponent; valid values are -1, 0, and 1}
}


\value{
\code{inverse()} returns a \code{TransferFunction} object obtained by 
swapping \code{fun} and \code{funinv} and 
swapping \code{domain} and \code{range}.
The names of the elements composing \code{TF} are changed appropriately.
If \code{TF} is not invertible, it returns \code{NULL}.

If \code{n=-1}, \code{TF^n} returns \code{inverse(TF)}.
If \code{n=1}, it returns \code{TF}.
If \code{n=0}, it returns the universal \code{\link{identity.TF}}.
For any other value of \code{n} it returns \code{NULL}.

\code{is.invertible()} returns \code{TRUE} or \code{FALSE}.
}    

\seealso{
\code{\link{identity.TF}}
}

