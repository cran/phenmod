\name{data.combine.stationNet}
\alias{data.combine.stationNet}

\title{Station net creator}
\description{Creates a net of the stations of dataset by checking the distances between the stations.}
\usage{data.combine.stationNet(dataset, range, 
					alt.range, silent=FALSE, 
					out2File=FALSE)}
\arguments{
  \item{dataset}{A dataset created by \code{\link{data.extract}} containing the information that 
			should be used to generate the station net.}
  \item{range}{The maximum distance between two stations that should be connected in the net.}
  \item{alt.range}{The maximum altitude difference between two stations that should be connected in the net.}
  \item{silent}{A boolean value determining wether the function should generate output messages or not.}
  \item{out2File}{A boolean value determining wether the output will be stored in log-files.}
 }

\details{Only stations with a maximal distance of 'range' and maximal altitude difference of \sQuote{alt.range} are connected.}
\value{Returns a list with entries for all stations of dataset. Each list entry stands for a station and contains all neighbours of that station.}
\author{Daniel Doktor, Maximilian Lange}

\examples{
  ## load extracted observations as created by 'data.extract'
  data(extractedObs)

  ## create station net
  stations.net <- data.combine.stationNet(extractedObs, 
				range=5000, alt.range=50, silent=FALSE, 
				out2File=FALSE)
}

