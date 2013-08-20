\name{result.pic.scatterplot}
\alias{result.pic.scatterplot}

\title{Scatterplot creation}
\description{Creates a scatterplot out of result values.}
\usage{result.pic.scatterplot(values, picPath=getwd(), 
			picName="beech-budburst", createFile=TRUE)}
\arguments{
  \item{values}{The values for which the scatterplot should be created.}
  \item{picPath}{The path where the png-files should be stored.}
  \item{picName}{The name of the created file with the scatterplot.}
  \item{createFile}{A boolean value determining wether a png-File will be created or not.}
 }
\details{Creates a scatterplot out of result values. The observed values are shown at the abscissa and the modelled values at the ordinate.}
\value{Returns nothing but stores the scatterplot as png-file with given path and filename.}
\author{Daniel Doktor, Maximilian Lange}
\seealso{\code{\link{result.main}}}

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

  ## create scatterplot
  result.pic.scatterplot(values=result.values, 
		picPath=getwd(), picName="beech_budburst",
		createFile=FALSE)
}

