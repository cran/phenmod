\name{result.extract.interpolate}
\alias{result.extract.interpolate}

\title{Result interpotion}
\description{Interpolates result values with given spatial information.}
\usage{result.extract.interpolate(mask.grid, values, alt, x, y)}
\arguments{
  \item{mask.grid}{The grid with spatial information the values are ordered by.}
  \item{values}{The values which should be interpolated.}
  \item{alt}{The related altitude for the gridcells of \sQuote{mask.grid}.}
  \item{x}{The related Rechtswert (Gauss-Krueger-coordinates) for the gridcells of \sQuote{mask.grid}.}
  \item{y}{The related Hochwert (Gauss-Krueger-coordinates) for the gridcells of \sQuote{mask.grid}.}
 }
\details{Interpolates result values with given spatial information by external drift kriging.}
\value{A vector with the interpolated values.}
\author{Daniel Doktor, Maximilian Lange}
\references{Krige, D.G., 1951. A statistical approach to some basic mine valuation problems on the witwatersrand. Journal of the Chemical, Metallurgical and Mining Society of South Africa 52, 119-139.
Pebesma, E.J., 2004. Multivariable geostatistics in S: the gstat package. Computers & Geosciences 30, 683-691.}
\seealso{\code{\link{result.extract.main}},\code{\link{data.main}}}

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

  ## interpolate result values with spatial informations of mask.grid
  values.model <- result.extract.interpolate(mask.grid=mask.grid, 
				values=values.model, alt=mask.grid$alt, 
				x=mask.grid$x, y=mask.grid$y)
}

