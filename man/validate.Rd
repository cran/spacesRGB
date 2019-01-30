% File validate.Rd

\name{validate}
\title{Validate a TransferFunction by applying some simple Tests}
\alias{validate}
\alias{validate.TransferFunction}


\description{
Each \code{TransferFunction} object is actually a list of so-called \emph{elementary} transfer functions;
for details on this see \code{\link{composition}()}.
This \code{validate()} applies an internal \code{validate()} function to each elementary function
individually.
The internal \code{validate()} function generates some points in the domain of the function
and checks that all points are transfered into the range of the function.
If the function is also invertible, it checks that the inverse transfers back to the original point.
}

\usage{
\S3method{validate}{TransferFunction}( TF, points=1300, tol=5.e-7, domain=NULL )
}

\arguments{
\item{TF}{a \code{TransferFunction} object with dimension N, and consisting of M elementary transfer functions}

\item{points}{the number of points to test, in each elementary function}

\item{tol}{the numerical tolerance for the inversion test - this is relative to the length of the corresponding side of the domain box}

\item{domain}{a 2xN matrix to use as an alternate domain, for the first elementary function in the list only.
\code{domain} can also be a vector of length 2, which is then replicated to a 2xN matrix.
}
}


\value{
The function returns a logical vector of length M.
The value of the i'th element is the validation status of the i'th elementary function.
The returned vector has the attribute \code{'message'} which is a list of length M
with explanatory text.
For nicely formatted text see \code{\link{print}()}.
}    

\seealso{
\code{\link{TransferFunction}},
\code{\link{identity.TF}},
\code{\link{composition}()},
\code{\link{print.TransferFunction}()}
}

