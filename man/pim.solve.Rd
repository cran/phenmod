\name{pim.solve}
\alias{pim.solve}

\title{Promoter-Inhibitor-Model}
\description{Applies a promoter-inhibitor-model to a given dataset.}
\usage{pim.solve(params, data, model.no=1, 
	silent=FALSE, out2File=FALSE)}
\arguments{
  \item{params}{The parameters for the promoter-inhibitor-model as list or vector with following order: a1, a2, a3, a4, T.min.i, T.opt.i, T.max.i, T.min.p, T.opt.p, T.max.p}
  \item{data}{A dataset containing the station ID, the station geographic and the Gauss-Krueger-coordinates of the station,
	the altitude of the station, the year and the julian day of the budburst and the leafcolouring,
	outlier information for budburst and leafcolouring and temperature data for modelling. Can be created by using the function \code{\link{data.main}}.}
  \item{model.no}{The promoter-inhibitor-model to use. See references for more details.}
  \item{silent}{A boolean value determining wether the function should generate output messages or not.}
  \item{out2File}{A boolean value determining wether the output will be stored in log-files.}
 }
\details{Applies a promoter-inhibitor-model with given parameters to a given dataset.}
\value{A dataset containing the values of the origin dataset and additionally the modelled budburst days.}
\references{Schaber, J. and Badeck, F.-W. (2003). Physiology-based phenology models for forest tree species in Germany
. International Journal of Biometeorology 47:193-201}
\author{Daniel Doktor, Maximilian Lange}
\seealso{\code{\link{data.main}}}

\examples{
  ## load preprocessed data
  data(dataFinal)

  ## set or load params
  params <- c(0, 0.058326, 0.109494, 0.039178, 
			-10.34, -0.89, 18.11,-10.03, 
			28.61, 44.49)

  ## apply model
  result <- pim.solve(params, dataFinal, model.no=11, 
			silent=FALSE, out2File=FALSE)

}

