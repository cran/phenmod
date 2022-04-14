\name{data.read.phase}
\alias{data.read.phase}

\title{Read data from TSV-table}
\description{This function reads phenological data out of a tsv-table.}
\usage{data.read.phase(path, filename)}
\arguments{
  \item{path}{The path of the tsv-file to read.}
  \item{filename}{The filename of the tsv-table.}
 }
\details{Reads phenological data out of a tsv-table containing the columns 
	\sQuote{DWD_STAT_ID} (the ID of the phenological station), 
	\sQuote{STAT_NAME} (the name of the phenological station),
	\sQuote{STAT_LON} (the longitude of the station),
	\sQuote{STAT_LAT} (the latitude of the station),
	\sQuote{STAT_ALT} (the altitude of the station),
	\sQuote{BEGIN_OBS} (the year the station started the observation),
	\sQuote{END_OBS} (the year the station stopped the observation),
	\sQuote{NATURRAUM_ID} (the id of the region),
	\sQuote{PHASE_ID} (the ID of the observed phase),
	\sQuote{OBS_DAY} (the julian day the observed phase starts),
	\sQuote{OBS_YEAR} (the year at which \sQuote{OBS_DAY} was observed),
	\sQuote{CHECKED} (a value determining wether the result was checked or not),
	\sQuote{outlier} (if result was checked, this value determines wether it was an outlier or not).
	Such a dataset can be obtained at \url{http://www.phenology.de}.}
\value{Returns the table as R-Dataframe.}
\references{Dierenbach, J., Badeck, F.W., Schaber, J., 2013. The plant phenological online database (PPODB): an online database for long-term phenological data. International Journal of Biometeorology , 1-8.}
\author{Daniel Doktor, Maximilian Lange}
\seealso{\code{\link{lcObs}},\code{\link{bbObs}}}
\examples{
  ##read file 'beech_budburst.tsv' out of working directory
  \dontrun{data.budburst <- data.read.phase(path=getwd(), 
				filename="beech_budburst.tsv")}

  ## resulting data looks like the following:
  data(lcObs)
  data(bbObs)

  lcObs
  bbObs
}

