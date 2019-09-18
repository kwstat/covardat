\name{alzheimers}
\alias{alzheimers}
\docType{data}
\title{
  Disease progression in an alzheimer's trial
}

\description{
  Disease progression in an alzheimer's trial
}

\usage{data("alzheimers")}
\format{
  A data frame with 480 observations on the following 4 variables.
  \describe{
    \item{\code{trt}}{Treatment factor: Low, High, Placebo}
    \item{\code{patient}}{Patient id number}
    \item{\code{month}}{Month number: 2,4,6,8,10,12}
    \item{\code{score}}{Score for disease}
  }
}

\details{
  Patents had two daily doses of placebo, low dose, or high dose of a
  treatment for Alzheimer's disease.

  Patients were assessed every 2 months for 1 year.

  Disease progression was measuered with the Alzheimer's Disease
  Assessment Scale, which evaluates memory, language, and praxis
  function.  Scores range from 0 to 70, with higher scores indicating
  greater severity of disease.

  Key question: Do the treatments have different response profiles over
  time?
  
}
\source{
  Glenn Walker and Jack Shostack. 2010.
  \emph{Common statistical methods for clinical research with {SAS}
    examples}.
  SAS publishing.
}

\references{
  None.
}
\examples{
data(alzheimers)
require("lattice")
xyplot(score ~ month|trt , data=alzheimers, groups=patient,
       xlab="month", ylab="ADAS cognitive score", type='l',layout=c(3,1),
       scales=list(x=list(rot=90, cex=.75, alternating=FALSE)))
}

\keyword{datasets}
