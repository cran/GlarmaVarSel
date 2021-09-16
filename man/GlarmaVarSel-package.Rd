\name{GlarmaVarSel-package}
\alias{GlarmaVarSel-package}
\alias{GlarmaVarSel}
\docType{package}
\title{
\packageTitle{GlarmaVarSel}
}
\description{
GlarmaVarSel consists of four functions: "variable_selection.R", "grad_hess_beta.R", "grad_hess_gamma.R" and "NR_gamma.R"
For further information on how to use these functions, 
we refer the reader to the vignette of the package.
}
\details{
GlarmaVarSel consists of four functions: "variable_selection.R", "grad_hess_beta.R", "grad_hess_gamma.R" and "NR_gamma.R"
For further information on how to use these functions, 
we refer the reader to the vignette of the package.
}
\author{
Marina Gomtsyan, Celine Levy-Leduc, Sarah Ouadah, Laure Sansonnet

Maintainer: Marina Gomtsyan <marina.gomtsyan@agroparistech.fr>
}
\references{
M. Gomtsyan et al. "Variable selection in sparse GLARMA models", arXiv:2007.08623v1
}
\keyword{package}
\examples{
n=50
p=30
X = matrix(NA,(p+1),n)
f = 1/0.7
for(t in 1:n){X[,t]<-c(1,cos(2*pi*(1:(p/2))*t*f/n),sin(2*pi*(1:(p/2))*t*f/n))}
gamma0 = c(0)
data(Y)
result = variable_selection(Y, X, gamma0, k_max=2, n_iter=100, method="min",
nb_rep_ss=1000, threshold=0.7, parallel=FALSE, nb.cores=1)
beta_est = result$beta_est
Estim_active = result$estim_active
gamma_est = result$gamma_est    
}
