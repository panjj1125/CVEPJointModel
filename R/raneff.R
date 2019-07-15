# A function to obtain the estimates of the random effects.
#
#
raneff <- function(res){
  theta <<- res$para$theta
  Sigmab <<- res$para$Sigmab
  TT_Sigmab <<- res$para$TT_Sigmab
  DS_Sigmab <<- res$para$TT_Sigmab
  rho <<- res$para$rho
  naLast_TT_data = res$compdata$TT
  naLast_DS_data = res$compdata$DS
  TT_DS_init = res$compdata$init
  TT_DS_obs = res$compdata$obs
  TT_DS_X = res$compdata$X
  TT_DS_Z = res$compdata$Z
  DS_mu = res$compdata$DS_mu; DS_st = res$compdata$DS_st
  TT_Tobs = res$compdata$TT_Tobs
  DS_Tobs = res$compdata$DS_Tobs
  TT_fixeff <<- res$model$TT$fixed ## fixed effects of target trait
  TT_raneff <<- res$model$TT$random ## random effects of target trait
  DS_fixeff <<- res$model$DS$fixed ## random effects of decision score
  DS_raneff <<- res$model$DS$random ## random effects of decision score
  Sigmab_all = mapply(SigmaB,naLast_TT_data,naLast_DS_data)
  sigma2_err_all = mapply(function(TT_dat,DS_dat){return(diag(rep(sigm2, c(nrow(TT_dat),nrow(DS_dat)))))},naLast_TT_data, naLast_DS_data)

  u <- mapply(slice_samp, TT_DS_init, TT_DS_obs, naLast_TT_data, naLast_DS_data, TT_DS_X, TT_DS_Z, DS_mu, DS_st, Sigmab_all, SIMPLIFY=FALSE)

  DS_u <- mapply(function(x,n1){
    if(is.null(x)){
      ## no missing
      return(NULL)
    } else {
      return(x[, (dim(x)[2]-sum(is.na(n1$DS_Trait)) + 1):dim(x)[2],drop=F])
    }
  }, u, naLast_DS_data)
  TT_u <- mapply(function(x, n1){
    #  a list of nsamp x nmiss (only T) resamples
    if(is.null(x)) {
      return(NULL)
    } else {
      return(x[, 1:sum(is.na(n1$TT_Trait)),drop=F])
    }
  }, u, naLast_TT_data)
  ## Update Sigmab
  TT_Ttrait <- mapply(function(t,u1){
    ## this gives you a matrix for target trait: observed value repeating nsample times and missing values imputed
    if (length(u1) == 0){
      return(t)
    }else{
      return(rbind(t, t(u1)))
    }
  }, TT_Tobs, TT_u, SIMPLIFY=FALSE)
  DS_Ttrait <- mapply(function(t,u1){
    ## this gives you a matrix for DS trait: observed value repeating nsample times and missing values imputed
    if (length(u1) == 0){
      return(t)
    }else{
      return(rbind(t, t(u1)))
    }
  }, DS_Tobs, DS_u, SIMPLIFY=FALSE)
  TT_DS_Ttrait <- mapply(function(tt, DS){return(rbind(tt, DS))}, TT_Ttrait, DS_Ttrait, SIMPLIFY=FALSE)

  b1_all <- mapply(
    function(Ttrait, x, z, Sigmab_TT_DS,sigma2_err){
      b1 <- matrix(apply(Ttrait, 2, function(t1){
        Sigmab_TT_DS %*% t(z) %*% solve(z %*% Sigmab_TT_DS %*% t(z) + sigma2_err) %*% (t1-x %*% theta)
      }),nrow=dim(z)[2]); rownames(b1) = colnames(Sigmab_TT_DS)
      return(b1)
    },
    TT_DS_Ttrait,TT_DS_X,TT_DS_Z,Sigmab_all,sigma2_err_all)
  b1 = sapply(b1_all,function(x){apply(x,1,mean)})
  return(b1)
}
