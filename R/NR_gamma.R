NR_gamma <-
function(Y, X, beta0, gamma0, n_iter=100){
  gamma<-gamma0
  new.params<-gamma
  n<-length(Y)
  iter<-1
  eps<-10^-6
  erreur<-1+eps
  
  while((erreur>eps)&&(iter<=n_iter)){
    old.params<-new.params
    res_grad_hess=grad_hess_gamma(Y, X, beta0, gamma)
    hessienne=res_grad_hess[[2]]
    grad=res_grad_hess[[1]]
    res_svd=svd(hessienne)
    U=res_svd$u
    V=res_svd$v
    ind_vp_not_null=which(round(res_svd$d,digits=6)!=0) 
    Lambda=diag(1/res_svd$d[ind_vp_not_null],length(ind_vp_not_null))
    hess_inv=V[,ind_vp_not_null]%*%Lambda%*%t(U[,ind_vp_not_null])
    gamma<-gamma-as.numeric(hess_inv%*%grad)
    
    new.params<-gamma
    
    erreur<-max(abs(new.params-old.params))
    iter<-iter+1
  }
  return(gamma=gamma)
}
