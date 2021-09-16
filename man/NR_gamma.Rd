\name{NR_gamma}
\alias{NR_gamma}
\title{
Newton-Raphson method for estimation of gamma
}
\description{
This function estimates gamma with Newton-Raphson method
}
\usage{
NR_gamma(Y, X, beta0, gamma0, n_iter)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
Observation matrix
}
  \item{X}{
Design matrix
} 
  \item{beta0}{
Initial beta vector
}
  \item{gamma0}{
Initial gamma vector
}
  \item{n_iter}{
Number of iterations of the algorithm. Default=100
}
}
\value{
\item{gamma}{
Estimated gamma vector
}
}
\references{
M. Gomtsyan et al. "Variable selection in sparse GLARMA models", arXiv:2007.08623v1
}
\author{
Marina Gomtsyan, Celine Levy-Leduc, Sarah Ouadah, Laure Sansonnet

Maintainer: Marina Gomtsyan <marina.gomtsyan@agroparistech.fr>
}
\examples{
n=50
p=30
X = matrix(NA,(p+1),n)
f = 1/0.7
for(t in 1:n){X[,t]<-c(1,cos(2*pi*(1:(p/2))*t*f/n),sin(2*pi*(1:(p/2))*t*f/n))}
gamma0 = c(0)
data(Y)
glm_pois<-glm(Y~t(X)[,2:(p+1)],family = poisson)
beta0<-as.numeric(glm_pois$coefficients)
gamma_est = NR_gamma(Y, X, beta0, gamma0, n_iter=100)
}
