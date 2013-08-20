\name{data.combine.clusters}
\alias{data.combine.clusters}

\title{Station cluster creation.}
\description{Creates a list of station clusters.}
\usage{data.combine.clusters(dataset, stations.net, 
			shuffle=TRUE, tries=100, 
			silent=FALSE, out2File=FALSE)}
\arguments{
  \item{dataset}{The dataset with the stations to cluster and their coordinates.}
  \item{stations.net}{A list containing the information which stations are neighbours. 
			Created by function \code{\link{data.combine.stationNet}}.}
  \item{shuffle}{A boolean value determining wether the stations should be processed in their order (value: FALSE) 
			or if they should be shuffled befor processing (value: TRUE). 
			Shuffled stations produce different results each run.}
  \item{tries}{If value of \sQuote{shuffle} is true, the integer value \sQuote{tries} determines how 
			much cluster-lists should be created. The cluster-list with the lowest 
			number of entries will be returned (this will increase the size of the clusters).}
  \item{silent}{A boolean value determining wether the function should generate output messages or not.}
  \item{out2File}{A boolean value determining wether the output will be stored in log-files.}
 }
\details{Creates a list of station clusters by randomly creating clusters (if value of \sQuote{shuffle} is TRUE) and 
		chosing the list with the lowest number of clusters.}
\value{Returns a list of station clusters (a station cluster is a vector containing related stations).}
\author{Daniel Doktor, Maximilian Lange}
\seealso{\code{\link{data.combine.clusters.search}},\code{\link{data.combine.stationNet}}}

\examples{
  ## load extracted observations as created by 'data.extract'
  data(extractedObs)

  ## create station net
  stations.net <- data.combine.stationNet(extractedObs, 
				range=5000, alt.range=50, silent=FALSE, 
				out2File=FALSE)

  ## search clusters in station net
  clusters <- data.combine.clusters(extractedObs, stations.net, 
				shuffle=TRUE, tries=3, silent=FALSE, 
				out2File=FALSE)
}

