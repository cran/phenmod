\name{result.extract.sub}
\alias{result.extract.sub}

\title{Extract values}
\description{Extracts and transformates given values to a data.frame with same order as a given masking grid.}
\usage{result.extract.sub(mask.grid, values, gk4.x, gk4.y, 
	outliers, silent=FALSE, withOutliers=FALSE)}
\arguments{
  \item{mask.grid}{A grid containing spatial information related to the values and their Gauss-Krueger-Coordinates. The values in the resulting data.frame will be ordered related to values in \sQuote{mask.grid}.}
  \item{values}{The values that should be extracted and transformed.}
  \item{gk4.x}{Gauss-Krueger-Rechtswert related to values.}
  \item{gk4.y}{Gauss-Krueger-Hochwert related to values.}
  \item{outliers}{Outlier information related to values.}
  \item{silent}{A boolean value determining wether the function should generate output messages or not.}
  \item{withOutliers}{A boolean value determining wether outliers should be included in extraction and transformation or not.}
 }
\details{Extracts and transformates given values (with related Gauss-Krueger-Coordinates) to a data.frame with same order as a given masking grid.}
\value{A grid containing the values and their coordinates (as given in \sQuote{mask.grid}).}
\author{Daniel Doktor, Maximilian Lange}
\seealso{\code{\link{result.extract.main}}}
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
 
  ## resolve outlier information
  outliers <- result$outlier.bb + result$outlier.lc
  outliers.na <- which(is.na(outliers)==TRUE)
  outliers[outliers.na] <- rep(0, length(outliers.na))

  mask.grid <- relatedGrid

  ## extract valid modelled values
  values.model <- result.extract.sub(mask.grid=mask.grid, 
				result$doy.bb.pim, result$gk4.x, 
				result$gk4.y, outliers=outliers, 
				silent=FALSE, withOutliers=FALSE)$values
}

