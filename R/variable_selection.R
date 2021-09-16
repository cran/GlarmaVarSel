variable_selection <-
function(Y, X, gamma0, k_max=2, n_iter=100, method="min", nb_rep_ss=1000,threshold=0.8, parallel=FALSE, nb.cores=1)
{
  search_non_null=function(x)
  {
    ind=as.numeric(which(x!=0))
    return(ind)
  }
  n=length(Y)
  p=dim(X)[1]-1
  q=length(gamma0)
  
  if((p+1)<n){
    glm_pois<-glm(Y~t(X)[,2:(p+1)],family = poisson)
    beta0<-as.numeric(glm_pois$coefficients)
  }else{
    cv_lambda<-cv.glmnet(t(X), Y, family = "poisson", alpha = 1)$lambda.min
    glm_pois<-glmnet(t(X), Y, family = "poisson", alpha = 1, lambda=cv_lambda)
    beta0<-as.numeric(glm_pois$beta)
    beta0[1] = as.numeric(glm_pois$a0)
  }
  gamma_est=NR_gamma(Y, X, beta0, gamma0, n_iter)
  
  for(k in 1:k_max)
  {
    grad_hess_res_est=grad_hess_beta(Y, X, beta0, gamma_est)
    Gradient<-grad_hess_res_est[[1]][1:(p+1)]
    Hessienne<-grad_hess_res_est[[2]]
    
    res_svd=svd(-Hessienne) 
    U=res_svd$u
    ind_vp_not_null=which(round(res_svd$d,digits=6)!=0)
    
    Lambda_rac=diag(sqrt(res_svd$d[ind_vp_not_null])) 
    Lambda_rac_inv=diag(1/sqrt(res_svd$d[ind_vp_not_null]))
    Y_eta=Lambda_rac_inv%*%t(U[,ind_vp_not_null])%*%Gradient+Lambda_rac%*%t(U[,ind_vp_not_null])%*%beta0
    X_eta=Lambda_rac%*%t(U[,ind_vp_not_null])
    
    if(method=='fast'){
      beta_all_lambda=glmnet(X_eta,Y_eta,alpha=1,family = "gaussian")$beta
      
      mat_ind_active=apply(beta_all_lambda,2,FUN = search_non_null)
      position_active=rep(0,(p+1))
      for (i in 1:length(mat_ind_active)){
        position_active[mat_ind_active[[i]]]=position_active[mat_ind_active[[i]]]+1
      }
      freq=position_active/(length(mat_ind_active)-1)
    }else if(method=='cv'){
      lambda_cv=cv.glmnet(X_eta,Y_eta,family="gaussian",alpha=1,parallel=TRUE)$lambda.min
      res.cum = rep(0,(p+1))
      b_sort_matrix = matrix(NA, nrow = nb_rep_ss, ncol = floor((p+1)/2))
      for(j in 1:nb_rep_ss){
        b_sort_matrix[j,] <- sort(sample(1:(p+1),floor((p+1)/2)))
      }
      for(j in 1:nb_rep_ss){
        b_sort = b_sort_matrix[j,]
        resultat_glmnet=glmnet(X_eta[b_sort,],Y_eta[b_sort],family="gaussian",alpha=1,lambda=lambda_cv)
        ind_glmnet=which(as.numeric(resultat_glmnet$beta)!=0)
        res.cum = res.cum + tabulate(ind_glmnet,(p+1))
      }
      freq=res.cum/nb_rep_ss
    }else{
      all_lambda=cv.glmnet(X_eta,Y_eta,family="gaussian",alpha=1,parallel=parallel)$lambda
      lambda_min = c(min(all_lambda))
      res.cum = rep(0,(p+1))
      b_sort_matrix = matrix(NA, nrow = nb_rep_ss, ncol = floor((p+1)/2))
      for(j in 1:nb_rep_ss){
        b_sort_matrix[j,] <- sort(sample(1:(p+1),floor((p+1)/2)))
      }
      for(j in 1:nb_rep_ss){
        b_sort = b_sort_matrix[j,]
        resultat_glmnet=glmnet(X_eta[b_sort,],Y_eta[b_sort],family="gaussian",alpha=1,lambda=lambda_min)
        ind_glmnet=which(as.numeric(resultat_glmnet$beta)!=0)
        res.cum = res.cum + tabulate(ind_glmnet,(p+1))
      }
      freq=res.cum/nb_rep_ss
    }
    
    Estim_active=which(freq>=threshold)
    
    beta_est=rep(0,(p+1))
    if(length(Estim_active) == 1){
      glm_pois<-glm(Y~t(X)[,Estim_active]-1,family = poisson)
      beta_est[Estim_active]<-as.numeric(glm_pois$coefficients)
    }else if((length(Estim_active)<n) & (rankMatrix(t(X)[,Estim_active])[1]==length(Estim_active))){
      glm_pois<-glm(Y~t(X)[,Estim_active]-1,family = poisson)
      beta_est[Estim_active]<-as.numeric(glm_pois$coefficients)
    }else{
      lambda_cv<-cv.glmnet(t(X)[,Estim_active], Y, family = "poisson", alpha = 0)$lambda.min
      glm_pois<-glmnet(t(X)[,Estim_active], Y, family = "poisson", alpha = 0, lambda=lambda_cv)
      beta_est[Estim_active]<-as.numeric(glm_pois$beta)
      if((1%in% Estim_active)==TRUE){
        beta_est[1]=as.numeric(glm_pois$a0)
      }
    }
    gamma_est=NR_gamma(Y, X, beta_est, gamma0, n_iter)
    gamma0 = gamma_est
    beta0=beta_est
    
  }
  return(list(estim_active=Estim_active,beta_est=beta_est,gamma_est=gamma_est))
}
