\name{result.main}
\alias{result.main}

\title{Result evaluation}
\description{Main function for result evaluation.}
\usage{result.main(mask.grid, result.grid, plant="beech", model="pim", 
			year=1954, picPath=getwd(), picName="beech-budburst", 
			createFiles=TRUE, rsquarePath=getwd(), 
			rsquareFile="rsquare.RData", rsquare.type="cod", 
			silent=FALSE, withOutliers=FALSE)}
\arguments{
  \item{mask.grid}{A grid with spatial information related to the resulting grid of modelling.}
  \item{result.grid}{The resulting grid of modelling. Can be obtained by using \code{\link{pim.solve}} or \code{\link{tsm.solve}}.}
  \item{plant}{The plant name for which the values of 'result.grid' are modelled.}
  \item{model}{A character value determining which model was used to create the results (either \sQuote{pim} or \sQuote{tsm}).}
  \item{year}{The processed year.}
  \item{picPath}{The path where the created png-files should be stored.}
  \item{picName}{The filename of the created png-files.}
  \item{createFiles}{A boolean flag determining wether the results should be stored in files or not.}
  \item{rsquarePath}{The path where the RData-file with the rsquare-dataset should be stored.}
  \item{rsquareFile}{The filename of the RData-file with the rsquare-dataset.}
  \item{rsquare.type}{The value of type (either \sQuote{cod} or \sQuote{pearson}) determines whether the coefficient of 
			determination or the squared pearson correlation coefficient is used as rsquare.}
  \item{silent}{A boolean value determining wether the function should generate output messages or not.}
  \item{withOutliers}{A boolean value determining wether outliers should be included in extraction and transformation or not.}
 }
\details{Main function for result evaluation. Extracts essential values from resulting grid, interpolates and masks values, creates histogramms, scatterplot and maps and calculates r-square.}
\author{Daniel Doktor, Maximilian Lange}
\seealso{\code{\link{result.extract.main}},\code{\link{result.extract.interpolate}},\code{\link{result.pic.histogramm}},\code{\link{result.pic.maps}},\code{\link{result.pic.scatterplot}},\code{\link{result.rsquare}},\code{\link{pim.solve}},\code{\link{tsm.solve}}}

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
  ## evaluate resuts
  result.main(mask.grid=relatedGrid, 
		result.grid=result, plant="beech", model="pim", 
		year="1952-2009", picPath=getwd(), picName="beech-budburst",
		createFiles=FALSE, rsquarePath=getwd(), 
		rsquareFile="rsquare.RData", rsquare.type="cod", 
		silent=FALSE, withOutliers=FALSE)
}

