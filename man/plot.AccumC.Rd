\name{plot.AccumC}
\alias{plot.AccumC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot \code{"AccumC"} Objects
}
\description{
This function plots objects of class \code{"AccumC"}. It makes small plots in the current graphics window and/or a big plot saved as a \code{.pdf} file.}
\usage{
plot.AccumC(x, plot = T, bigplot = F, box_col = "gray", prettytitle = NA, filename = "AccumC_plot", ratelabel = T, log_y_axis = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class \code{"AccumC"}
}
  \item{plot}{
Logical. Make plots in the current graphics window?
}
  \item{bigplot}{
Logical. Make a complex plot saved as a \code{".pdf"} file?
}
  \item{box_col}{
the fill color of boxes and bars; default is \code{"gray"}
}
  \item{prettytitle}{
A title for the \code{bigplot} that is different from \code{corename}
}
  \item{filename}{
A filename or path for the \code{bigplot}
}
  \item{ratelabel}{
Logical. Annotate the age/depth model plot with (text) sedimentation rates?
}
  \item{log_y_axis}{
Logical. Should the y axes on the carbon accumulation rate plots be log?
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
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
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
plot(
  AccumC(imv), bigplot = F
  )
}
