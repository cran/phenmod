\name{data.combine.clusters.search}
\alias{data.combine.clusters.search}

\title{Cluster search}
\description{Searches clusters in a station net.}
\usage{data.combine.clusters.search(stations, stations.net, 
			shuffle=TRUE)}
\arguments{
  \item{stations}{A list of all stations.}
  \item{stations.net}{A list of neighbourstations (neighbours of the station at same list position in list \sQuote{stations}).}
  \item{shuffle}{A boolean value determining wether the stations should be processed in their order (value: FALSE) 
	or if they should be shuffled befor processing (value: TRUE).}
 }
\details{Searches related stations in a station net. If value of \sQuote{shuffle} is TRUE, 
		each run of this method produces different results.}
\value{Returns a list containing related stations.}
\author{Daniel Doktor, Maximilian Lange}

\examples{
  ## load extracted observations as created by 'data.extract'
  data(extractedObs)

  ## create station net
  stations.net <- data.combine.stationNet(extractedObs, range=5000, 
				alt.range=50, silent=FALSE, out2File=FALSE)

  ## search clusters in station net
  stations <- as.list(unique(extractedObs$STAT_ID))
  clusters <- data.combine.clusters.search(stations, 
			stations.net, shuffle=TRUE)

  
}

