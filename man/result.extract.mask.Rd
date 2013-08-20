\name{result.extract.mask}
\alias{result.extract.mask}

\title{Result masking}
\description{Masks values with the spatial informations of a given grid.}
\usage{result.extract.mask(mask.grid, values)}
\arguments{
  \item{mask.grid}{The grid with the spatial informations which should be used for masking.}
  \item{values}{The values to mask (in same order as \sQuote{mask.grid}).}
 }
\details{Masks values with the spatial informations of a given grid. Gridcells with NA value are masked in the resulting data.frame by setting them to \sQuote{-9999}.}
\value{A vector of masked values.}
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

  ## mask result values with spatial informations of mask.grid
  values.model <- result.extract.mask(relatedGrid, 
				values.model)
}

