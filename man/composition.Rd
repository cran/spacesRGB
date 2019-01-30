% File composition.Rd

\name{composition}
\title{The composition of TransferFunction objects}
\alias{composition}
\alias{composition.TransferFunction}
\alias{*.TransferFunction}
\alias{\%X\%.TransferFunction}
\alias{\%;\%.TransferFunction}
\alias{\%O\%.TransferFunction}
\alias{\%X\%}
\alias{\%;\%}
\alias{\%O\%}
\alias{is.identity}
\alias{is.identity.TransferFunction}
\alias{identity.TF}

\description{
The function \code{composition(TF1,TF2)} returns a \code{TransferFunction} that is 
\code{TF1} followed by \code{TF2}.
Four equivalent infix operators are also available.

}


\usage{
\S3method{composition}{TransferFunction}( TF1, TF2 )

\S3method{*}{TransferFunction}( TF1, TF2 )
\S3method{\%;\%}{TransferFunction}( TF1, TF2 )
\S3method{\%X\%}{TransferFunction}( TF1, TF2 )
\S3method{\%O\%}{TransferFunction}( TF2, TF1 )

identity.TF

\S3method{is.identity}{TransferFunction}( TF )
}

\arguments{
\item{TF1}{a \code{TransferFunction} object}

\item{TF2}{a \code{TransferFunction} object}

\item{TF}{a \code{TransferFunction} object}
}

\details{
In order to be composed, the dimensions of \code{TF1} and \code{TF2} must be equal,
or the dimension of one of them must be 1.
In the latter case, the function is applied to each coordinate in exactly the same way.

All the above represent the function \code{TF1} followed by \code{TF2}.
In mathematics this operation is usually called \emph{composition of functions}
(and \emph{composition of morphisms} in category theory),
and in computer science and BT.2100 and BT.2390 it is called the \emph{concatenation}.
In BT.2390 it is also called the \emph{cascade}.

The ACES literature uses infix notation with the symbol \code{'+'} 
which is unfortunate because in mathematics
the plus symbol is only used for commutative operations, which composition certainly is not.
The symbol \code{'*'} is offered here as an alternative,
since \code{'*'} does not imply commutativity (e.g. as in MATLAB's matrix multiplication).
In computer science the symbol \code{';'} is common, and so \code{\%;\%} is offered as an alternative.
In BT.2100 and BT.2390 the symbol \otimes is used, and so \code{\%X\%} is offered as an alternative.
And finally, in mathematics \circ is used but in the opposite order,
so that \code{TF2 \%O\% TF1} is identical to \code{composition(TF1,TF2)}.
\cr
Each \code{TransferFunction} object is actually a list of so-called \emph{elementary} transfer functions.
If \code{TF1} has \eqn{M_1} elementary functions and \code{TF2} has \eqn{M_2} elementary functions,
then \code{composition(TF1,TF2)} has \eqn{\le M_1 + M_2} elementary functions.
It can be strictly less if there is cancellation of elementary functions at the end of \code{TF1}
and the beginning of \code{TF2}.
}

\value{
\code{composition(TF1,TF2)} returns a \code{TransferFunction} object,
which applies \code{TF1} followed by \code{TF2}.
The individual objects \code{TF1} and \code{TF2} are stored  inside the returned object.
In case of ERROR it returns \code{NULL}.
The 4 infix operators above all invoke \code{composition()}.

\code{identity.TF} is a built-in global \code{TransferFunction} object which is a 
universal identity for composition.
This means that for any \code{TransferFunction}  \code{TF},  
\code{TF*identity.TF = identity.TF*TF = TF}.
Moreover, \code{TF*TF^-1 = TF^1*TF = identity.TF}.
This is *not* the same as \code{base::identity()}.

\code{is.identity(TF)} tests whether \code{TF} is the universal identity,
and returns \code{TRUE} or \code{FALSE}.
}


\references{
Technical Bulletin.
TB-2018-002.
ACES Output Transform Details.
June 2018 (draft).

ACES Retrospective and Enhancements
March 2017.

BT.2100.
Image parameter values for high dynamic range television for use in production
and international programme exchange.
June 2017.

BT.2390.
High dynamic range television for production and international programme exchange.
April 2018.
}

\examples{
comp = power.OOTF(2.2) * power.OOTF(1.4)
x = 0:100 / 100
max( abs( transfer(comp,x) - transfer(power.OOTF(2.2*1.4),x) ) )  #  1.110223e-16

comp * comp^-1
## This is a universal identity TransferFunction.

is.identity(comp * comp^-1)          # TRUE

identical( comp * identity.TF, comp )  # TRUE
}


\seealso{
\code{\link{TransferFunction}},
\code{\link{transfer}()},
\code{\link{inverse}()}
}

