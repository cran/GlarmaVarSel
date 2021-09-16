\name{grad_hess_beta}
\alias{grad_hess_beta}
\title{
Gradient and Hessian of the log-likelihood with respect to beta
}
\description{
This function calculates the gradient and Hessian of the log-likelihood with
respect to beta.
}
\usage{
grad_hess_beta(Y, X, beta0, gamma0)
}
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
}

\value{
\item{grad_L_beta}{
Vector of the gradient of L with respect to beta
}
\item{hess_L_beta}{
Matrix of the Hessian of L with respect to beta
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
result = grad_hess_beta(Y, X, beta0, gamma0)
grad = result$grad_L_beta
Hessian = result$hess_L_beta
}

