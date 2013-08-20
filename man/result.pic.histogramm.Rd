\name{result.pic.histogramm}
\alias{result.pic.histogramm}

\title{Histogramm Creation}
\description{Creates histogramms of results.}
\usage{result.pic.histogramm(values, picPath=getwd(), 
		picName="budburst-beech", silent=FALSE,
		createFile=TRUE)}
\arguments{
  \item{values}{The data.frame with result values (can be created with \code{\link{result.extract.main}}).}
  \item{picPath}{The path where the png-files should be stored.}
  \item{picName}{The name of the created files with the histogramms.}
  \item{silent}{A boolean value determining wether the function should generate output messages or not.}
  \item{createFile}{A boolean value determining wether a png-File will be created or not.}
 }
\details{Creates histogramms of results (observed value, modelled value and difference of them).}
\value{Returns nothing but creates histogramms as png-files with given path and filename.}
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
			interpolate=FALSE, silent=FALSE)

  ## create histogramm
  result.pic.histogramm(values=result.values, 
		picPath=getwd(), picName="beech_budburst", 
		silent=FALSE, createFile=FALSE)
}

