\name{result.rsquare}
\alias{result.rsquare}

\title{R-square calculation.}
\description{Calculates the r-square value of a model.}
\usage{result.rsquare(values, type="cod")}
\arguments{
  \item{values}{A list containing the modelled (doy.model) and the observed (doy.observed) values.}
  \item{type}{The value of type (either \sQuote{cod} or \sQuote{pearson}) determines whether the coefficient of 
			determination or the squared pearson correlation coefficient is calculated.}
 }
\details{Calculates the r-square value of modelled values with given observed values.}
\value{The r-square as numeric value.}
\author{Daniel Doktor, Maximilian Lange}
\seealso{\code{\link{result.main}}}

\examples{
 	modelled <- c(100,102,98,97,96)
	observed <- rep(100,5)
	values <- list(doy.model=modelled, doy.observed=observed)
	result.rsquare(values, type="cod")
}

