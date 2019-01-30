\name{Reference Rendering Transform}
\alias{RRT.TF}
\alias{general.RRT}


\title{
Reference Rendering Transform
}

\description{
The fixed \code{RRT.TF} corresponds to the RRT in aces-dev 1.1.
\cr
A parameterized version \code{general.RRT()} is also provided - for experimentation.
This one returns a \code{TransferFunction} with the argument values "locked-in".
}


\usage{
RRT.TF
general.RRT( glowmod="1.1", redmod="1.1" )
}

\arguments{
\item{glowmod}{the version of the Glow Modifier to use.  
The only version currently supported is \code{"1.1"}.
\cr
\code{glowmod} can also be \code{NULL}, \code{NA}, or \code{FALSE},
which means to use no Glow Modifier at all.
}

\item{redmod}{the version of the Red Modifier to use.  
The only version currently supported is \code{"1.1"}.
This string can also be \code{"1.1+pinv"}
which means to use a precision inverse; the forward transfer is exactly the same.
This precision inverse uses an iterative root-finder, 
and is slower than the approximate default inverse.
\cr
\code{redmod} can also be \code{NULL}, \code{NA}, or \code{FALSE},
which means to use no Red Modifier at all.
}
}

\details{
\code{RRT.TF} is a \code{Transferfunction} that maps ACES RGB to OCES RGB.
Both spaces are relative to the AP1 primaries.
\code{RRT.TF} is constructed by calling \code{general.RRT()} with its default arguments.
The transfer is complicated; here is a summary of the steps
starting with ACES RGB as input:
\enumerate{
\item glow module (see argument \code{glowmod})
\item red modifier (see argument \code{redmod})
\item matrix conversion from AP0 RGB \rarrow AP1 RGB
\item clamp to non-negative RGB
\item global desaturation
\item segmented spline, applied to each channel separately
\item matrix conversion from AP1 \rarrow  AP0 (now OCES RGB)
}
}

\value{
\code{general.RRT()} returns a \code{Transferfunction} that maps ACES RGB to OCES RGB.
The domain is [0,47000]\eqn{^3} and the range is [0,10000]\eqn{^3}.
}


\references{
ST 2065-1:2012.
SMPTE Standard - Academy Color Encoding Specification (ACES).
2013.
}

\seealso{
\code{\link{TransferFunction}}
}
