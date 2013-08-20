\name{data.combine}
\alias{data.combine}

\title{Main function to combine timeseries}
\description{This function creates a station net and builds clusters of stations out of it. These clusters are used to create combined timeseries.}
\usage{data.combine(dataset, range=5000, alt.range=50, shuffle=TRUE, 
				tries=100, silent=FALSE, out2File=FALSE,
				clusters.tmp.file="tmpcluster.RData")}
\arguments{
  \item{dataset}{A dataset created by \code{\link{data.extract}} containing the information that 
			should be used to generate the combined timeseries.}
  \item{range}{The maximum distance between two stations that should be connected in the station net.}
  \item{alt.range}{The maximum altitude difference between two stations that should be connected in the station net.}
  \item{shuffle}{A boolean value determining wether the stations should be processed in their order (value: FALSE) 
			or if they should be shuffled befor processing (value: TRUE). 
			Shuffled stations produce different results each run.}
  \item{tries}{If value of \sQuote{shuffle} is true, the integer value \sQuote{tries} determines how 
			much cluster-lists should be created. The cluster-list with the lowest 
			number of entries will be returned (this will increase the size of the clusters).}
  \item{silent}{A boolean value determining wether the function should generate output messages or not.}
  \item{out2File}{A boolean value determining wether the output will be stored in log-files.}
  \item{clusters.tmp.file}{A file where the clusters are saved for evaluation. If the value is NULL, no file will be created.}
 }
\details{This function joins the functions \code{\link{data.combine.stationNet}}, \code{\link{data.combine.clusters}} and \code{\link{data.combine.timeseries}}.}
\value{A dataset containing the combined timeseries as a data.frame with same columns like a data.frame created by \code{\link{data.extract}}.}
\references{Schaber J., Badeck F. (2002). Evaluation of methods for the combination of phenological time series and outlier detection. Tree Physiology, 22:973-982}
\author{Daniel Doktor, Maximilian Lange}
\seealso{\code{\link{data.combine.clusters}},\code{\link{data.combine.stationNet}},\code{\link{data.combine.timeseries}}}

\examples{
  ## load extracted observations as created by 'data.extract'
  data(extractedObs)
  
  ## combine timeseries
  data.combined <- data.combine(dataset=extractedObs, range=5000, 
				alt.range=50, shuffle=TRUE, tries=3, 
				silent=FALSE, out2File=FALSE, 
				clusters.tmp.file=NULL)
}

