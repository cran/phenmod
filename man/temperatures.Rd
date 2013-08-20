\name{temperatures}
\docType{data}
\alias{temperatures}
\title{Temperature data}
\description{
This dataset contains interpolated temperature data based on data obtained from WebWerdis, 
the Web-based weather request and distribution system of the German Weather Service (DWD).
The interpolation was done via external drift kriging provided by the \sQuote{gstat} package.
}
\usage{extractedObs}
\format{A matrix containing 4 years (366 days) of temperature data at 9 different locations given by \code{\link{relatedGrid}}.}
\source{WebWerdis (DWD)}
\references{Krige, D.G., 1951. A statistical approach to some basic mine valuation problems on the witwatersrand. Journal of the Chemical, Metallurgical and Mining Society of South Africa 52, 119-139.
Pebesma, E.J., 2004. Multivariable geostatistics in S: the gstat package. Computers & Geosciences 30, 683-691.}
\seealso{\code{\link{relatedGrid}}}
\keyword{datasets}
