\name{data.addTemperatures}
\alias{data.addTemperatures}

\title{Add temperature data to a dataset.}
\description{Adds a temperature vector to a dataset containing spatial and seasonal phenological information.}
\usage{data.addTemperatures(dataset, grid.related.to.Temperatures,
		temperature.filenames, 
		temperature.matrix, temperature.scale.factor, 
		out2File=FALSE, silent=FALSE)}
\arguments{
  \item{dataset}{A dataset containing spatial and seasonal phenological information. Can be created by using function \code{\link{data.extract}}.}
  \item{grid.related.to.Temperatures}{A grid containing spatial informations for the temperature data.}
  \item{temperature.filenames}{The full name of temperature files. The filenames have to contain the year (YYYY) of the temperature observation. 
	The files should be stored as RData-Files and have to contain a matrix called \sQuote{edk.one.year} with 366 columns (one per day) 
	and the number of rows equal to the number of rows in \sQuote{grid.related.to.Temperatures}. 
	Should have the value NULL if a temperature matrix should be used instead.}
  \item{temperature.matrix}{An array containing temperature data. The year of the observation should given as rowname, the columns should equal the julian day 
	of the observation and the third dimension of the array should equal the location given in \sQuote{grid.related.to.Temperatures}}.
	The matrix will be used instead of temperature files if \sQuote{temperature.filenames} is NULL.
  \item{temperature.scale.factor}{The down-scaling factor for the temperature data (needed if the data is scaled).}
  \item{out2File}{A boolean value determining wether the output will be stored in log-files.}
  \item{silent}{A boolean value determining wether the function should generate output messages or not.}
 }
\details{Adds a temperature vector to a dataset containing spatial (The coordinates of the station as Gauss-Krueger coordinates) and seasonal (The year for which the budburst day should be modelled and the leafcolouring day of the previous year) phenological information.}
\value{A dataset containing spatial and seasonal phenological data and the added temperatures.}
\author{Daniel Doktor, Maximilian Lange}
\seealso{\code{\link{data.loadTemperature}},\code{\link{data.extract}},\code{\link{relatedGrid}},\code{\link{temperatures}}}

\examples{
  ## load extracted observations as created by 'data.extract'
  data(extractedObs)
  ## load temperature data
  data(temperatures)
  ## load grid with spatial informations for temperature files
  data(relatedGrid)	

  ## add temperatures from files in working directory
  dataset <- data.addTemperatures(dataset=extractedObs, 
				out2File=FALSE,
				grid.related.to.Temperatures=relatedGrid,
				temperature.filenames=NULL, 
				temperature.matrix=temperatures,
				temperature.scale.factor=0.1)
 
  ## resulting dataset should equal 'dataFinal'
  data(dataFinal)
  dataFinal
  dataset
}

