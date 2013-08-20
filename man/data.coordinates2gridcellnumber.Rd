\name{data.coordinates2gridcellnumber}
\alias{data.coordinates2gridcellnumber}

\title{Search number of gridcell for given coordinates.}
\description{Searches the number of the cell of a given spatial grid by given coordinates (Gauss-Krueger-Coordinates).}
\usage{data.coordinates2gridcellnumber(grid, x,y)}
\arguments{
  \item{grid}{The grid in which the cell should be searched.}
  \item{x}{The \sQuote{Rechtswert} of the Gauss-Krueger-Coordinates.}
  \item{y}{The \sQuote{Hochwert} of the Gauss-Krueger-Coordinates.}
 }
\details{Searches the number of the cell of a given spatial grid by given coordinates (Gauss-Krueger-Coordinates).}
\value{Returns the cell-number of the grid-cell which is next to the given coordinates.}
\author{Daniel Doktor, Maximilian Lange}

\examples{
  ## load grid containing spatial information
  data(relatedGrid)

  ## search position in grid
  xy <- util.geoco2gk(54.12,10.17,4)
  gridposition <- data.coordinates2gridcellnumber(grid=relatedGrid, 
	x=xy[,1],y=xy[,2])

  gridposition
}

