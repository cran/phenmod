\name{data.main}
\alias{data.main}

\title{Main function for data conversion.}
\description{This function extracts data from tsv-Files, creates combined timeseries out of them if needed and adds temperature data.}
\usage{data.main(extraction.done=TRUE, out2File=FALSE, 
			grid.related.to.temperaturefiles, 
			valid.years=1952:2009,
			combine.time.series=TRUE, range=10000, 
			alt.range=50, temperature.scale.factor=1,
			dataPath=getwd(),temperature.filenames,
			temperature.matrix, pathForTmpFiles=getwd(),
			pathToSave=getwd(), plant="beech")}
\arguments{
  \item{extraction.done}{If data extraction was already done, turn this value to TRUE and the 
	extracted values will be loaded from \sQuote{pathForTmpFiles}.}
  \item{out2File}{A boolean value determining wether the output will be stored in log-files.}
  \item{grid.related.to.temperaturefiles}{A grid containing spatial informations for the temperature files.}
  \item{valid.years}{The years that should be included in extraction process 
			(may be needed if some data in datasets are insufficient).}
  \item{combine.time.series}{A boolean value determining wether the extracted timeseries should be combined or not.}
  \item{range}{The maximum distance between two stations that should be connected in the station station-net 
	needed for the timeserie-combination.}
  \item{alt.range}{The maximum altitude difference between two stations that should be connected in the 
	station-net needed for the timeserie-combination.}
  \item{temperature.scale.factor}{The down-scaling factor for the temperature data (needed if the data is scaled).}
  \item{dataPath}{The path where the tsv-tables are stored. These tsv-tables should have the following name: 
	\sQuote{<PLANT>_budburst.tsv} for budburst-data or \sQuote{<PLANT>_leafcolouring.tsv} for leafcolouring-data, 
	where <PLANT> stands for the name of the plant, for example \sQuote{beech}.}
   \item{temperature.filenames}{The full name of temperature files. The filenames have to contain the year (YYYY) of the temperature observation. 
	The files should be stored as RData-Files and have to contain a matrix called \sQuote{edk.one.year} with 366 columns (one per day) 
	and the number of rows equal to the number of rows in \sQuote{grid.related.to.Temperatures}.
	Should have the value NULL if a temperature matrix should be used instead.}
  \item{temperature.matrix}{An array containing temperature data. The year of the observation should given as rowname, the columns should equal the julian day 
	of the observation and the third dimension of the array should equal the location given in \sQuote{grid.related.to.Temperatures}}.
	The matrix will be used instead of temperature files if \sQuote{temperature.filenames} is NULL.
  \item{pathForTmpFiles}{The path where the extracted data should be temporarily stored as RData-files.}
  \item{pathToSave}{The path where the resulting dataframe should be stored as RData-file with the name 
	\sQuote{<PLANT>-dataset.RData} for non-combined data or \sQuote{<PLANT>-dataset-cts.RData} for combined data.}
  \item{plant}{The name of the plant which should be processed.}
 }
\details{This function is the main method for data extraction. It exctracts budburst and leafcolouring data from tsv-Files, 
	creates combined timeseries out of them if needed and adds temperature data from RData files containing daily mean temperatures.}
\value{Returns nothing, but saves the resulting dataset in \sQuote{pathToSave}.}
\author{Daniel Doktor, Maximilian Lange}
\seealso{\code{\link{data.extract}},\code{\link{data.combine}},\code{\link{data.addTemperatures}}}

\examples{
  ## load temperature data
  data(temperatures)
  ## load grid with spatial informations for temperature files
  data(relatedGrid)
  
  ## extract dataset from datafiles in working directory
  \dontrun{dataset <- data.main(extraction.done=FALSE, out2File=FALSE, 
			grid.related.to.temperaturefiles=relatedGrid,
			combine.timeseries=TRUE, range=5000, alt.range=50, 
			temperature.scale.factor=0.1, dataPath=getwd(), 
			temperature.filenames=NULL,
			temperature.matrix=temperatures,
			pathForTmpFiles=getwd(), pathToSave=getwd(),
			plant="beech")}
}

