\name{util.isLeapYear}
\alias{util.isLeapYear}

\title{Leap year check}
\description{Checks wether a given year is a leap year or not.}
\usage{util.isLeapYear(year)}
\arguments{
  \item{year}{A year or vector of years to check.}
 }
\details{This function checks wether a year (yyyy) is a leap year or not.}
\value{Returns a boolean vector determining wether the given years are leap years (TRUE) or not (FALSE).}
\author{Daniel Doktor, Maximilian Lange}

\examples{
  util.isLeapYear(2000:2012)
}

