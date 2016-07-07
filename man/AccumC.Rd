\name{AccumC}
\alias{AccumC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Accum(ulate) C(arbon)
}
\description{
AccumC is a method for estimating the rate of carbon accumulation in peatlands, especially those where many reworked macrofossils have been dated. 
}
\usage{
AccumC(corename = NA, bin_size = NA, age_bw = NA, kernel = "r", calcurve = intcal13, it_power = 17, pct_C = 0.48)
}

\arguments{
  \item{x}{a \code{data.frame} where \code{names(x) == c("labcode", "depth", "bp", "std", "afbd")} }
  \item{corename}{
The case-sensitive name of the directory and filenames i.e. \code{"corename/corename.dat"} and \code{"corename/corename.afbd"}
}
  \item{bin_size}{
size of age bandwith in years. If \code{NA} (default), it is set automatically.
}
  \item{age_bw}{
size of age bandwith in years. If \code{NA} (default), it is set automatically.
}
  \item{kernel}{
kernel shape for PDF smoothing, as in \code{stats::density()}
}
  \item{calcurve}{
The radiocarbon calibration curve to be used. Built-in choices are \code{intcal13} or \code{shcal13}, but any \code{data.frame} with columns \code{cal} and \code{bp} can be used.
}
  \item{it_power}{
The Monte Carlo simulation of accumulation rates will be run \code{2^it_power} times.
}
  \item{pct_C}{
Fraction carbon in ash-free bulk density. A default value of \code{0.48} is used. Replace this value with \code{1} if fraction carbon is used instead of ash-free bulk density.
}
}
\details{
If \code{x} is used, \code{corename} is merely a text string that is used as the name of the core, but if \code{x} is not specified, \code{corename} is assumed to be the name of a folder inside R's working directory that contains a comma-separated file called "\code{corename.dat}".
}
\value{
A \code{list} is returned with the following components:
%%  This function returns a \code{list} of class \code{AccumC} with the following items:
%%  If it is a LIST, use
	\item{$metadata}{settings and descriptive information about the dataset}
	\item{$calibrations}{a \code{data.frame} with columns \code{labcodes}, \code{depths}, and \code{maxLikeAge} }
	\item{$rates}{a \code{data.frame} with the age of the midpoint, top, and bottom of each depth bin, its sediment accumulation rate, }
	\item{$total_prob}{a \code{list} of two vectors describing the PDF of the age of sediment in the core}
	\item{$age_mid_quantiles}{Quantiles of the probable age of depth bin midpoints, where \code{ p = c(0, 0.25, 0.5, 0.75, 1)}}
	\item{$age_test}{an \code{array} with the results of the Monte Carlo simulation of the age of the top and bottom of each depth bin.}
	\item{$rate_test}{an \code{array} with the results of the Monte Carlo simulation of the sedimentation rate of each depth bin.}
	\item{$crate_test}{an \code{array} with the results of the Monte Carlo simulation of the carbon accumulation rate of each depth bin.}
	\item{$all_pdfs}{a \code{list} containing the PDF of the calibrated age of each radiocarbon date}
	\item{$age_quantiles}{Quantiles of the probable age of depth ends, where \code{ p = c(0, 0.25, 0.5, 0.75, 1)}}
	\item{$rate_quantiles}{Quantiles of the probable carbon accumulation rate of each depth bin, where \code{ p = c(0, 0.25, 0.5, 0.75, 1)}}
	\item{$binned_pdfs}{the PDF of the age of each depth bin}
%% ...
}
\references{
Nichols et al., 2016
}
\author{
Jonathan E. Nichols\cr
Lamont-Doherty Earth Observatory\cr
Columbia University

jnichols@ldeo.columbia.edu
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{plot.AccumC}}
}
\examples{
data(imv)
x <- AccumC(imv)
plot(x, bigplot = F)
}

\keyword{ paleoecology }
\keyword{ carbon cycle }