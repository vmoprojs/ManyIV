\name{SJIVEfitF}
\alias{SJIVEfitF}
\encoding{UTF-8}
\title{IV Estimation and Inference}
\description{
  It will compute the SJIVE/SJEF estimator and the associated variance covariance matrix
}
\usage{
SJIVEfitF(y=NULL,X1=NULL,X2=NULL,Z1=NULL,Z2=NULL,full=2)
}

\arguments{
\item{y}{Blah. The default is \code{NULL}}
\item{X1}{Blah. The default is \code{NULL}}
\item{X2}{Blah. The default is \code{NULL}}
\item{Z1}{Blah. The default is \code{NULL}}
\item{Z2}{Blah. The default is \code{NULL}}
\item{full}{Blah. The default is \code{NULL}. See \bold{Details}.}
}

\details{
  The optimization method is \code{Blah}.

  The \code{full} parameter allows to... 


  The \code{Z1} parameter allows to...
}

\value{
  Returns an object of class \code{IVEstim}.
  An object of class \code{IVEstim} is a list containing
  at most the following components:
  \item{param}{The vector of  parameters' estimates;}
}


\references{
  P. A. Bekker & F. Crudu (2015). "Jackknife instrumental variable estimation with heteroskedasticity",\emph{Journal of Econometrics}, \bold{185}, 332--342.
}



\author{Federico Crudu, \email{federico.crudu@gmail.com},\url{https://sites.google.com/site/federicocrudu/},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/}
}
  
\examples{
library(ManyIV)
# ST: Data Simulation:
y<-rnorm(500)#log(packs)  # dependent variable
X1<- matrix(rnorm(500),500,1)#cbind(log(rprice),log(rincome))# matrix of endogenous variables  
Z1<- matrix(rnorm(10000),500,20)#cbind(log(rincome) , tdiff , I(tax/cpi))# non included instruments
Z2<- matrix(1,length(y),1)# included instruments
# END: Data Simulation:


# ST: Estimation
sol = SJIVEfitF(y,X1,X2,Z1,Z2,full=2)
# END: Estimation
str(sol)
plot(density(sol$beta))
}
\keyword{SJIVEfit}
