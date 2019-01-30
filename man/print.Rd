% File validate.Rd

\name{print}
\title{Print Basic Facts about a TransferFunction}
\alias{print}
\alias{print.TransferFunction}


\description{
Each \code{TransferFunction} object is actually a list of so-called \emph{elementary} transfer functions;
for details on this see \code{\link{composition}()}.
This \code{print()} calls an internal \code{print()} function for each elementary function individually.
The internal \code{print()} also calls an 
internal \code{validate()} (with default arguments) which runs some basic tests
and formats the results nicely for printing, see \code{\link{validate}()}.
}

\usage{
\S3method{print}{TransferFunction}( x, ... )
}

\arguments{
\item{x}{a \code{TransferFunction} object consisting of M \emph{elementary} transfer functions}

\item{...}{further arguments ignored, but required by the generic \code{print()}}
}


\value{
The function returns \code{TRUE} or \code{FALSE}.
}    

\seealso{
\code{\link{TransferFunction}},
\code{\link{validate}()},
\code{\link{composition}()}
}

\examples{
tf = sRGB.EOTF^-1 * power.EOTF(2.5)
tf
##  #--------------------   [sRGB.EOTF]^-1    ---------------------#
##  [sRGB.EOTF]^-1 is a univariate TransferFunction.
##  domain:      [0,1]  (linear display)
##  range:       [0,1]  (non-linear signal)
##  invertible:  Yes
##  orientation: preserving
##  range-test points = 1300, max(distance)=0.
##  validation:  Passed
##  #--------------------   power.EOTF(2.5)    ---------------------#
##  power.EOTF(2.5) is a univariate TransferFunction.
##  domain:      [0,1]  (non-linear signal)
##  range:       [0,1]  (linear display)
##  invertible:  Yes
##  orientation: preserving
##  range-test points = 1300, max(distance)=0.
##  validation:  Passed
}

