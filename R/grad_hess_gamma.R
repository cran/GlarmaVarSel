  grad_hess_gamma <-
  function(Y, X, beta0, gamma0)
  {
    mu<-numeric(n) 
    E<-numeric(n) 
    Z<-numeric(n) 
    W<-numeric(n) 
    q<-length(gamma0)
    p<-length(X[,1])-1
    n<-length(Y)
    
    grad_W=matrix(0,nrow=q,ncol=n)
    hess_W_liste=list()
    hess_L_gamma=matrix(0,ncol=q,nrow=q)
    
    Z[1]=0
    W[1]=beta0%*%X[,1]
    mu[1]=exp(W[1])
    E[1]=Y[1]/mu[1]-1
    hess_W_liste[[1]]=matrix(0,ncol=q,nrow=q)
    
    for (t in 2:n)
    {
      hess_W_liste[[t]]=matrix(NA,ncol=q,nrow=q)
      jsup=min(q,(t-1))
      W[t]=beta0%*%X[,t]+sum(gamma0[1:jsup]*E[(t-1):(t-jsup)])
      for (l in 1:q)
      {
        if (l<=jsup) 
        {
          grad_W[l,t]=E[(t-l)]-sum(gamma0[1:jsup]*(1+E[(t-1):(t-jsup)])*grad_W[l,(t-1):(t-jsup)])
        }
        for (m in l:q)
        {
          if (m+l<t) 
          {
            A1=-(1+E[(t-l)])*grad_W[m,(t-l)]
            A2=-(1+E[(t-m)])*grad_W[l,(t-m)]
          }
          else 
          {
            A1=0
            A2=0
          }
          A=A1+A2
          B=-(1+E[(t-1):(t-jsup)])*grad_W[l,(t-1):(t-jsup)]*grad_W[m,(t-1):(t-jsup)]
          C=(1+E[(t-1):(t-jsup)])*as.numeric((lapply(hess_W_liste[((t-1):(t-jsup))],'[',l,m)))
          hess_W_liste[[t]][l,m]=A-sum(gamma0[1:jsup]*(B+C))
        }
      }
      mu[t]=exp(W[t])
      E[t]=Y[t]/mu[t]-1
      
      terme1_hess=hess_W_liste[[t]]*(Y[t]-exp(W[t]))
      terme2_hess=exp(W[t])*(grad_W[,t]%*%t(grad_W[,t]))
      hess_L_gamma=hess_L_gamma+terme1_hess-terme2_hess
    }
    hess_L_gamma[lower.tri(hess_L_gamma)]=hess_L_gamma[upper.tri(hess_L_gamma)]
    grad_L_gamma=grad_W%*%(Y-exp(W))
    return(list(grad_L_gamma = grad_L_gamma, hess_L_gamma=hess_L_gamma))
  }
