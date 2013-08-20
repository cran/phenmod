\name{result.pic.maps}
\alias{result.pic.maps}

\title{Map creation}
\description{Creates maps out of results.}
\usage{result.pic.maps(values, picPath=getwd(), 
	picName="beech-budburst", silent=FALSE,
	createFile=TRUE)}
\arguments{
  \item{values}{The data.frame with result values (can be created with \code{\link{result.extract.main}}).}
  \item{picPath}{The path where the png-files should be stored.}
  \item{picName}{The name of the created files with the maps.}
  \item{silent}{A boolean value determining wether the function should generate output messages or not.}
  \item{createFile}{A boolean value determining wether a png-File will be created or not.}
 }
\details{Creates maps out of results.}
\value{Returns nothing but stores maps as png-files with given path and filename.}
\author{Daniel Doktor, Maximilian Lange}
\seealso{\code{\link{result.extract.main}},\code{\link{result.main}}}

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
			interpolate=TRUE, silent=FALSE)

  ## not enough successfully calculated budbursts
  ## replace with examples
  result.values$doy.model <- round(rnorm(n=9,mean=100,sd=20))
  result.values$doy.observed <- round(rnorm(n=9,mean=100,sd=20))

  ## create maps
  result.pic.maps(values=result.values, 
		picPath=getwd(), picName="beech_budburst", 
		silent=FALSE, createFile=FALSE)
}

