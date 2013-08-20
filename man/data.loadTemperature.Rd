\name{data.loadTemperature}
\alias{data.loadTemperature}

\title{Loads a temperature vector}
\description{Loads a temperature vector of given length from a given day of a given year.}
\usage{data.loadTemperature(year, temperature.year, temperature.previous.year, 
				from.previous.year.doy, length, position, 
				scale.factor=0.1)}
\arguments{
  \item{year}{The year of the budburst doy which should be calculated.}
  \item{temperature.year}{A vector of temperature data (with length 365 or 366) of year given by \sQuote{year}.}
  \item{temperature.previous.year}{A vector of temperature data (with length 365 or 366) of previous year.}
  \item{from.previous.year.doy}{The day the modelling starts (leafcolouring day of previous year).}
  \item{length}{The length of the resulting temperature vector.}
  \item{position}{The number of the grid-cell for which the temperatures should be loaded.}
  \item{scale.factor}{The down-scaling factor for the temperature data (needed if the data is scaled).}
 }
\details{Loads a temperature vector of given length from temperature data. The day of leafcolouring of the 
previous year should be given as starting day. The vector contains one temperature datapoint per day and 
should have a length of 300 or more to ensure a stable modelling process. Internal function used by method \code{\link{data.addTemperatures}}.}
\value{A vector with temperature data for a given location and a given period.}
\author{Daniel Doktor, Maximilian Lange}
\seealso{\code{\link{data.addTemperatures}},\code{\link{data.coordinates2gridcellnumber}}}

\examples{
  t.year <- as.matrix(rnorm(365, mean=10, sd=5))
  t.p.year <- as.matrix(rnorm(365, mean=9, sd=6))
  temperatures <- data.loadTemperature(year=2006,temperature.year=t.year, 
			temperature.previous.year=t.p.year, 
			from.previous.year.doy=210, length=365, 
			position=1, scale.factor=0.1)
}

