\name{data.extract}
\alias{data.extract}

\title{Essential data extraction.}
\description{Extracts the essential data out of budburst and leafcolouring data and combines them.}
\usage{data.extract(data.budburst, data.leafcolouring, 
			valid.years=1952:2009, out2File=FALSE, 
			silent=FALSE)}
\arguments{
  \item{data.budburst}{A dataset extracted by function \code{\link{data.read.phase}} containing budburst 
				information of currently processed plant.}
  \item{data.leafcolouring}{A dataset extracted by function \code{\link{data.read.phase}} containing leafcolouring 
				information of currently processed plant.}
  \item{valid.years}{The years that should be included in extraction process 
			(may be needed if some data in datasets are insufficient).}
  \item{out2File}{A boolean value determining wether the output will be stored in log-files.}
  \item{silent}{A boolean value determining wether the function should generate output messages or not.}
 }
\details{This function extracts the essential data out of budburst and leafcolouring data and combines them. 
	Therefor it extracts the stations contained in both datasets and searches the leafcolouring dates and 
	the related budburst date (at the following year).}
\value{A dataset containing the station ID, the station geographic and the Gauss-Krueger-coordinates of the station,
	the altitude of the station, the year and the julian day of the budburst and the leafcolouring and 
	outlier information for budburst and leafcolouring.}
\author{Daniel Doktor, Maximilian Lange}
\seealso{\code{\link{data.read.phase}},\code{\link{lcObs}},\code{\link{bbObs}},\code{\link{extractedObs}}}

\examples{
  ## load budburst observations (extracted via 'data.read.phase')
  data(bbObs)
  ## load leafcolouring observations (extracted via 'data.read.phase')
  data(lcObs)

  ## extract essential data
  data.extracted <- data.extract(data.budburst=bbObs, 
				data.leafcolouring=lcObs, 
				out2File=FALSE)
}

