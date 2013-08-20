\name{util.geoco2gk}
\alias{util.geoco2gk}

\title{Geographic coordinates to Gauss-Krueger-Coordinates}
\description{Converts geographic coordinates to Gauss-Krueger-Coordinates.}
\usage{util.geoco2gk(x,y, meridian=4)}
\arguments{
  \item{x}{A vector of longitudes to transform.}
  \item{y}{A vector of latitudes to transform.}
  \item{meridian}{The referenced meridian for Gauss-Krueger-Coordinates.}
 }
\details{Converts geographic coordinates (longitude, latitude) to Gauss-Krueger-Coordinates (\sQuote{Rechtswert}, \sQuote{Hochwert}) referenced by a given meridian.}
\value{A matrix containing the the \sQuote{Rechtswert} of the Gauss-Krueger-Coordinates in its first column and the \sQuote{Hochwert} in its second column.}
\author{Daniel Doktor, Maximilian Lange}

\examples{
 x <- c(51.3, 54.7)
 y <- c(12.3, 13.5)
 gk.xy <- util.geoco2gk(x,y,4)
}

