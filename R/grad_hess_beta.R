grad_hess_beta <-
function(Y, X, beta0, gamma0)
{
  mu<-numeric(n) 
  E<-numeric(n) 
  Z<-numeric(n) 
  W<-numeric(n) 
  q<-length(gamma0)
  p<-length(X[,1])-1
  n<-length(Y)
  
  grad_W_beta=matrix(NA,nrow=(p+1),ncol=n)
  hess_W_beta_liste=list()
  hess_L_beta=matrix(0,ncol=(p+1),nrow=(p+1))
  
  Z[1]=0
  W[1]=beta0%*%X[,1]
  mu[1]=exp(W[1])
  E[1]=Y[1]/mu[1]-1
  hess_W_beta_liste[[1]]=matrix(0,ncol=(p+1),nrow=(p+1))
  grad_W_beta[,1]=X[,1]
  
  for (t in 2:n)
  {
    hess_W_beta_liste[[t]]=matrix(NA,ncol=(p+1),nrow=(p+1))
    jsup=min(q,(t-1))
    W[t]=beta0%*%X[,t]+sum(gamma0[1:jsup]*E[(t-1):(t-jsup)])
    for (k in 1:(p+1))
    {
      grad_W_beta[k,t]=X[k,t]-sum(gamma0[1:jsup]*(1+E[(t-1):(t-jsup)])*grad_W_beta[k,(t-1):(t-jsup)])
      
      for (j in k:(p+1))
      {
        B=(1+E[(t-1):(t-jsup)])*grad_W_beta[k,(t-1):(t-jsup)]*grad_W_beta[j,(t-1):(t-jsup)]
        C=-(1+E[(t-1):(t-jsup)])*as.numeric((lapply(hess_W_beta_liste[((t-1):(t-jsup))],'[',k,j)))
        hess_W_beta_liste[[t]][k,j]=sum(gamma0[1:jsup]*(B+C))
      }
    }
    mu[t]=exp(W[t])
    E[t]=Y[t]/mu[t]-1
    
    terme1_hess=hess_W_beta_liste[[t]]*(Y[t]-exp(W[t]))
    terme2_hess=exp(W[t])*(grad_W_beta[,t]%*%t(grad_W_beta[,t]))
    hess_L_beta=hess_L_beta+terme1_hess-terme2_hess
  }
  hess_L_beta[lower.tri(hess_L_beta)]=hess_L_beta[upper.tri(hess_L_beta)]
  grad_L_beta=grad_W_beta%*%(Y-exp(W))
  return(list(grad_L_beta=grad_L_beta, hess_L_beta=hess_L_beta))
}
