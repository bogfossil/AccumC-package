\name{AccumC-package}
\alias{AccumC-package}
\alias{AccumC}
\docType{package}
\title{
\packageTitle{AccumC}
}
\description{
\packageDescription{AccumC}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{AccumC}
\packageIndices{AccumC}

This package contains all of the functions and data to recreate the analyses from the publication, "A probabilistic method of assessing carbon accumulation rate at Imnavait Creek Peatland, Arctic LTER, Alaska." These include all the radiocarbon and ash-free bulk density measurements reported in the paper, the functions \code{AccumC()} and \code{plot.AccumC()} that are used to analyze the data and plot the results of the analyses. For convenince, the package also includes the IntCal13 and SHCal13 radiocarbon calibration datasets.
}
\author{
\packageAuthor{AccumC}

Maintainer: \packageMaintainer{AccumC}
}
\references{
Nichols et al., JQS, 2016
}
\keyword{ package }
\seealso{
}
\examples{
data(imv)
x <- AccumC(imv)
plot(x, bigplot = F)
}
