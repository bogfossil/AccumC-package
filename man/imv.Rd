\name{imv}
\alias{imv}
\docType{data}
\title{
Imnavait Creek Radiocarbon Data
}
\description{
Radiocarbon dates with ash-free bulk density measurements from the Imnavait Creek Peatland cores IMV11 and IMV12-B.
}
\usage{data("imv")}
\format{
  A data frame with 31 observations on the following 5 variables.
  \describe{
    \item{\code{labcode}}{the lab codes for the radiocarbon dates}
    \item{\code{depth}}{the depths of the radiocarbon dates in centimeters}
    \item{\code{bp}}{radiocarbon years bp}
    \item{\code{std}}{lab-reported standard deviation for radiocarbon measurements}
    \item{\code{afbd}}{ash-free bulk density (organic matter density) in g/cc}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
Nichols et al., 2016
}
\examples{
x <- AccumC(imv)
}
\keyword{datasets}
