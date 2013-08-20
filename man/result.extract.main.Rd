\name{result.extract.main}
\alias{result.extract.main}

\title{Essential result extraction}
\description{Extracts essential results from model.}
\usage{result.extract.main(mask.grid, result.grid, model="pim", 
	interpolate=TRUE, silent=FALSE, withOutliers=FALSE)}
\arguments{
  \item{mask.grid}{A grid containing spatial information related to values (which should be extracted) and their Gauss-Krueger-Coordinates. The values in the resulting data.frame will be ordered related to values in mask.grid.}
  \item{result.grid}{The grid created by \code{\link{pim.solve}} containing observated and modelled values.}
  \item{model}{A character value determining which model was used to create the results (either \sQuote{pim} or \sQuote{tsm}).}	
  \item{interpolate}{A boolean value determining wether the results should be interpolated (with spatial information of mask.grid) or not.}
  \item{silent}{A boolean value determining wether the function should generate output messages or not.}
  \item{withOutliers}{A boolean value determining wether outliers should be included in extraction and transformation or not.}
 }
\details{Extracts essential results (observed and modelled results and their difference and coordinates) from used model.}
\value{A data.frame containing the modelled values (\sQuote{doy.model}), the observed values (\sQuote{doy.observed}), their difference (\sQuote{doy.dif}) and related coordinates (\sQuote{x},\sQuote{y})}
\author{Daniel Doktor, Maximilian Lange}
\seealso{\code{\link{result.extract.sub}},\code{\link{result.extract.mask}},\code{\link{result.extract.interpolate}}}
\examples{
  ## load preprocessed data
  data(dataFinal)
  ## load spatial information
  data(relatedGrid)

  ## set or load params
  params <- c(0, 0.058326, 0.109494, 0.039178, 
			-10.34, -0.89, 18.11,-10.03, 
			28.61, 44.49)

  ## apply model
  result <- pim.solve(params, dataFinal, model.no=11, 
			silent=FALSE, out2File=FALSE)
 
  ## extract essential results
  result.values <- result.extract.main(
				mask.grid=relatedGrid, 
				result.grid=result, model="pim", 
				interpolate=FALSE, silent=FALSE)
}

