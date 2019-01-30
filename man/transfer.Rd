% File transfer.Rd

\name{transfer}
\title{Apply TransferFunction to a Vector or an Array}
\alias{transfer}
\alias{transfer.TransferFunction}




\description{
The function \code{transfer()} applies the given \code{TransferFunction} to the
given vector or array \code{x} and returns a numeric object of the same dimensions.
}


\usage{
\S3method{transfer}{TransferFunction}( TF, x, domaincheck=TRUE )
}

\arguments{
\item{TF}{a \code{TransferFunction} object, with dimension N}

\item{x}{a numeric vector or array.
If N\eqn{\ge}2 then \code{x} must be an MxN matrix, or a vector that can be converted to such
a matrix, by row.
}

\item{domaincheck}{check whether numbers or rows of \code{x} are in the domain box of \code{TF} before application}
}


\value{
Let N := \code{dimension(TF)}.
\cr
If N=1 then \code{x} can have any length or dimension;
the function is applied to each number in \code{x} in a vectorized way,
and the returned object is then assigned the same dimensions as \code{x}.  
If \code{x} is a matrix then the returned object is assigned the same \code{rownames}.
If a number is \code{NA} then the returned number is also \code{NA}.
\cr
If N\eqn{\ge}2 and \code{x} is an MxN matrix,
then the function is applied to each row of \code{x} individually and the returned object
is a matrix with the same dimensions and \code{rownames} as \code{x}.
If any number in a row is \code{NA} then the returned row is all \code{NA}s.
\cr
If \code{TF} is a universal identity (e.g. \code{identity.TF}), the function returns \code{x}
with no checking.
\cr
In case of a global error (e.g. dimension mismatch) the function returns \code{NULL}.
}    

\seealso{
\code{\link{TransferFunction}},
\code{\link{identity.TF}}
}

