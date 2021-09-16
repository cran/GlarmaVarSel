\name{variable_selection}
\alias{variable_selection}
\title{
Variable selection
}
\description{
This function performs variable selection, estimates a new vector beta and a new vector gamma
}
\usage{
variable_selection(Y, X, gamma0, k_max = 2, n_iter = 100, method = "min", 
  nb_rep_ss = 1000, threshold = 0.8, parallel = FALSE, nb.cores = 1)
}
\arguments{
  \item{Y}{
Observation matrix
}
  \item{X}{
Design matrix
}
  \item{gamma0}{
Initial gamma vector
}
  \item{k_max}{
Number of iteration to repeat the whole algorithm
}
  \item{n_iter}{
Number of iteration for Newton-Raphson algorithm
}
  \item{method}{
Stability selection method: "fast", "min" or "cv". In "min" the smallest lambda is chosen, in "cv" cross-validation lambda is chosen for stability selection. "fast" is a fater stability selection approach. The default is "min"
}
  \item{nb_rep_ss}{
Number of replications in stability selection step. The default is 1000
}
  \item{threshold}{
Threshold for stability selection. The default is 0.9
}
  \item{parallel}{
Whether to parallelize stability selection step or not. The default is FALSE
}
  \item{nb.cores}{
Number of cores for parallelization. The default is 1
}
}
\value{
\item{estim_active }{Estimated active coefficients}
\item{beta_est}{Vector of estimated beta values}
\item{gamma_est}{Vector of estimated gamma values}
%% ...
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
result = variable_selection(Y, X, gamma0, k_max=2, n_iter=100, method="min",
nb_rep_ss=1000, threshold=0.7, parallel=FALSE, nb.cores=1)
beta_est = result$beta_est
Estim_active = result$estim_active
gamma_est = result$gamma_est    
}