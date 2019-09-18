\name{tumor}
\alias{tumor}
\docType{data}
\title{
  Growth of colon carcinoma tumors in mice.
}
\description{
  Growth of colon carcinoma tumors in mice.
}
\usage{data("tumor")}
\format{
  A data frame with 330 observations on the following 4 variables.
  \describe{
    \item{\code{group}}{Factor for group, 3 levels}
    \item{\code{mouse}}{Mouse id number}
    \item{\code{day}}{Day}
    \item{\code{size}}{Tumor size}
  }
}
\details{
  Three therapies were tested against colon cancer in mice.
  Mice were injected with mouse colon carcinoma cells and tumors were
  measured at times from 7 to 21 days.

  Some mice died before the end of the experiment.

  The tumoz size was measured in cubic millimeters.
}
\source{
  James A. Koziol, Donna A. Maxwell, Matsuro Fukushima,
  M. E. M. Colmerauer and Yosef H. Pilch. (1981).
  A Distribution-Free Test for Tumor-Growth Curve Analyses with
  Application to an Animal Tumor Immunotherapy Experiment.
  \emph{Biometrics}, 37, 383-390.
  http://www.jstor.org/stable/2530427

  Used with permission of James Koziol.
}
\references{
  Also appears in Dawson 1997 (with different mouse numbering).
}
\examples{
data(tumor)
library(lattice)
xyplot(size~day|group, data=tumor, groups=mouse, type='l',
       main="Tumor growth over time for each mouse",
       xlab="Day", ylab="Tumor size",
       layout=c(3,1))
}
\keyword{datasets}
