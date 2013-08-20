\name{tsm.solve}
\alias{tsm.solve}

\title{Temperature-Sum-Model}
\description{Applies a temperature-sum-model to a given dataset.}
\usage{tsm.solve(params, data, silent=FALSE, out2File=FALSE)}
\arguments{
  \item{params}{The parameters for the temperature-sum-model as list or vector with following order: \eqn{T_b}{T_b}, \eqn{F^*}{F^*}. \eqn{T_b}{T_b} is the threshold temperature. Above this value, temperatures are summed for calculating the development rate. \eqn{F^*}{F^*} is the threshold development rate. If the development rate reached this threshold, the phase occurs.}
  \item{data}{A dataset containing the station ID, the station geographic and the Gauss-Krueger-coordinates of the station,
	the altitude of the station, the year and the julian day of the budburst and the leafcolouring,
	outlier information for budburst and leafcolouring and temperature data for modelling. Can be created by using the function \code{\link{data.main}}.}
  \item{silent}{A boolean value determining wether the function should generate output messages or not.}
  \item{out2File}{A boolean value determining wether the output will be stored in log-files.}
 }
\details{Applies a temperature-sum-model with given parameters to a given dataset.}
\value{A dataset containing the values of the origin dataset and additionally the modelled budburst days.}
\references{Menzel, A. (1997). Phaenologie von Waldbaeumen unter sich aendernden Klimabedingungen - Auswertung der Beobachtungen in den Internationalen Phaenologischen Gaerten und Moeglichkeiten der Modellierung von Phaenodaten. Thesis. Forstwissenschaftliche Fakultaet der Uni Muenchen. Muenchen, Universitaet Muenchen.}
\author{Daniel Doktor, Maximilian Lange}
\seealso{\code{\link{data.main}}}

\examples{
  ## load preprocessed data
  data(dataFinal)

  ## set or load params
  params <- c(0, 100)

  ## apply model
  result <- tsm.solve(params, dataFinal, 
		silent=FALSE, out2File=FALSE)

}

