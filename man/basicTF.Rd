\name{Basic TransferFunctions}
\alias{sRGB.EOTF}
\alias{BT.709.EOTF}
\alias{BT.2020.EOTF}
\alias{ProPhotoRGB.EOTF}
\alias{SMPTE.240M.EOTF}

\docType{data}

\title{
Basic TransferFunctions
}

\description{
\tabular{ll}{
\code{sRGB.EOTF}        \tab the standardized sRGB transfer function \cr
\code{BT.709.EOTF}      \tab the standardized BT.709 transfer function \cr
\code{BT.2020.EOTF}     \tab the standardized BT.2020 transfer function \cr
\code{ProPhotoRGB.EOTF} \tab the standardized ProPhotoRGB transfer function \cr
\code{SMPTE.240M.EOTF}  \tab the standardized SMPTE-240M transfer function \cr
}
}



\details{
All of these are built-in \code{\link{TransferFunction}} objects;
they have no parameters and are ready-to-go.
All are EOTFs and have domain and range the interval [0,1], and all are monotone increasing.
All are defined in 2 pieces, with a linear segment near 0.
All are easily inverted.
}



\references{
Wikipedia. \bold{sRGB}.
\url{http://en.wikipedia.org/wiki/SRGB}.

BT.709.
Parameter values for the HDTV standards for production and international programme exchange.
June 2015.

BT.2020.
Parameter values for ultra-high definition television systems for production and international programme exchange.
October 2015.

Wikipedia. \bold{ProPhoto RGB}.
\url{http://en.wikipedia.org/wiki/ProPhoto_RGB_color_space}.


ANSI/SMPTE 240M-1995.
SMPTE STANDARD for Television Signal Parameters 1125-Line High-Definition Production Systems.
}


\seealso{
\code{\link{TransferFunction}} 
}

\examples{
# make plot comparing 5 EOTFs
colvec	= c('black','red','blue','green','orange')
plot( sRGB.EOTF, color=colvec[1], main="The Basic 5 EOTFs" )
plot( BT.709.EOTF, color=colvec[2], add=TRUE )
plot( BT.2020.EOTF, color=colvec[3], add=TRUE )
plot( ProPhotoRGB.EOTF, color=colvec[4], add=TRUE )
plot( SMPTE.240M.EOTF, color=colvec[5], add=TRUE )
legend( 'topleft', legend=c('sRGB','BT.709','BT.2020','ProPhotoRGB','SMPTE.240M'), 
             col=colvec, bty='n', lty=1, lwd=2 )
}

\keyword{datasets}
