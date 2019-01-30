% File TransferFunction.Rd

\name{TransferFunction}
\title{Constructing and Testing TransferFunction Objects}
\alias{TransferFunction}
\alias{is.TransferFunction}
\alias{as.TransferFunction}
\alias{as.TransferFunction.default}

\description{
The function \code{TransferFunction()} is the constructor for \bold{TransferFunction} objects.

\code{is.TransferFunction()} tests whether an object is a valid \bold{TransferFunction} object.\cr
\code{as.TransferFunction()} converts other variables to a \bold{TransferFunction} object,
and is designed to be overridden by other packages.
}



\usage{
TransferFunction( fun, funinv, domain, range, id=NULL )

is.TransferFunction(x)

\S3method{as.TransferFunction}{default}( \ldots )
}

\arguments{
  \item{fun}{a function that accepts a numeric argument, and returns one of the same length.
  The dimension of \code{fun} is determined by the arguments \code{domain} and \code{range}.
  The function must be \emph{injective} and this is checked if the function is univariate.  
  The requirements for univariate and multivariate functions are very different, see \bold{Details}.  }

  \item{funinv}{a function that the inverse for \code{fun}.
  If \code{fun} is univariate and \code{funinv=NULL}, then an approximation for the inverse is computed
  using \code{stats::splinefun()}.
  If \code{fun} is multivariate and \code{funinv=NULL}, then it is an ERROR. }

  \item{domain}{
  a 2xN matrix, or a numeric vector that can be converted to such a matrix, by column.
  In each column, the entry in row 1 must be strictly less than the entry in row 2.  
  The columns of \code{domain} define N intervals whose product is a box in \eqn{R^N}
  that is the domain of \code{fun}.
  The box must be finite.
  If N=1 then the box is just an interval, and \code{fun} is univariate.
  Otherwise, \code{fun} is multivariate with dimension N.
  }
  
  \item{range}{
  a 2xN matrix, or a numeric vector that can be converted to such a matrix, by column.
  The N here must be equal to the N for \code{domain}.
  The matrix defines a box that encloses the image of \code{domain} under \code{fun}.
  The box must be finite.
  }
  
  \item{id}{
  a character string that is helpful when printing the object, and in logging messages.
  If \code{id=NULL} then an appropriate string is created from the function call.
  }
  
  \item{x}{an \R object to test for being a valid \bold{TransferFunction} object.}
  
  \item{...}{arguments for use in other packages.}  
}
    
\value{
\code{TransferFunction()} returns a \bold{TransferFunction} object, or \code{NULL} in case of ERROR.

\code{is.TransferFunction()} returns \code{TRUE} or \code{FALSE}.
It only checks the \code{class}, using \code{base::inherits()}.

\code{as.TransferFunction.default()} issues an ERROR message and returns \code{NULL}..
}    


\details{
If \code{fun} is univariate, then it must be able to accept a numeric vector of any length,
and apply the function to each number in the vector; i.e. \code{fun} must be vectorized.
If a number in the vector is NA, then the function must silently return NA for that number;
usually this is not a problem.
The function is *not* required to test whether the number is in the \code{domain} interval;
this is handled by the \code{TransferFunction} code.
\cr
If \code{fun} is multivariate with dimension N, then it must be able to accept a vector
of length N and return a vector of length N.
It is *not* required to accept an MxN matrix.
It is *not* required to test whether the vector is in the domain box.
\cr
The function \code{funinv} has the same requirements as \code{fun}.
}


\seealso{
\code{\link{dimension}()},
\code{\link{composition}()},
\code{\link{is.invertible}()},  
\code{\link{metadata}()},
\code{\link{inverse}()},
\code{\link{transfer}()},
\code{\link{orientation}()},
\code{\link{validate}()},  
\code{\link{print.TransferFunction}()},
\code{\link{plot.TransferFunction}()}
}

\examples{
#  make a test TransferFunction

myfun = function(x) {x*x}

test  = TransferFunction( myfun, sqrt, domain=c(0,3), range=c(0,9), id='test.TF' )

#  print it
test
#--------------------   test.TF    ---------------------#
## test.TF is a univariate TransferFunction.
## domain:      [0,3]  (x)
## range:       [0,9]  (y)
## invertible:  Yes
## orientation: preserving
## range-test points = 1300, max(distance)=0.
## validation:  Passed

# and now plot it
plot( test )
}


\keyword{TransferFunction}
