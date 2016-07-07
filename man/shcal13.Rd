\name{shcal13}
\alias{shcal13}
\docType{data}
\title{
SHCal13}
\description{
The SHCal13 radiocarbon calibration curve
}
\usage{data("shcal13")}
\format{
  A data frame with 5141 observations on the following 5 variables.
  \describe{
    \item{\code{cal}}{calibrated years BP}
    \item{\code{bp}}{radiocarbon years BP}
    \item{\code{error}}{calibration error in years BP}
    \item{\code{delta}}{delta 14C in per mille}
    \item{\code{sigma}}{per mille}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
http://www.radiocarbon.org/IntCal13%20files/shcal13.14c
}
\references{
Alan G Hogg, Quan Hua, Paul G Blackwell, Caitlin E Buck, Thomas P Guilderson Timothy J  Heaton, Mu Niu, Jonathan G Palmer, Paula J Reimer, Ron W Reimer, Christian S M Turney, Susan R H Zimmerman Radiocarbon 55(4). DOI: 10.2458/azu_js_rc.55.16783}
\examples{
data(shcal13)
View(shcal13)
}
\keyword{datasets}
