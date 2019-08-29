# Joint Modelling Procedure in Analyzing Highly Unbalanced CVEP Data
#
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

CVEP_JM_Year <- function(MET_dat,factors,TT_mm=list(fixed=c("Year","Loc"),random=c("Variety","Variety:Year","Variety:Loc")),DS_mm=list(fixed=c("Year","Loc"),random=c("Variety")),window=Inf,yes_full=F,converg_control=list(nsamp=500,max.iter=1000,err=10^(-7),err1=10^(-4),seed=20190421)){
  
  library(lme4)
  library("dplyr")
  library(parallel)
  library(mvtnorm)
  
  TT_fixeff <<- TT_mm$fixed ## fixed effects of target trait
  TT_raneff <<- TT_mm$random ## random effects of target trait
  DS_fixeff <<- DS_mm$fixed ## random effects of decision score
  DS_raneff <<- DS_mm$random ## random effects of decision score
  nsamp <<- converg_control$nsamp
  max.iter = converg_control$max.iter
  err = converg_control$err
  err1 = converg_control$err1
  seed1 <<- converg_control$seed
  
  dat0 = MET_dat
  dat0[,factors] = apply(dat0[,factors],2,as.character)
  yes.exist = is.element(factors,colnames(dat0))
  if (!all(yes.exist)) stop(paste("Error: There are no factors:",paste(factors[!yes.exist],"in the data",collapse=",")))
  alleff <- unique(c(TT_fixeff,TT_raneff,DS_fixeff,DS_raneff))
  index_inter = grep(":",alleff) # the index of interaction
  if (length(index_inter)>0){
    ## if there are interaction terms
    eff_int = unlist(strsplit(alleff[index_inter],split=":"))
    eff_model = unique(alleff[-index_inter],eff_int)
  } else eff_model = alleff
  yes.exist = is.element(eff_model,factors)
  if (!all(yes.exist)) stop(paste("Error: There are no factors:",paste(eff_model[!yes.exist],collapse=",")))
  
  if (!is.element("TT_Trait",colnames(dat0))) stop("Error: There is no TT_Trait in the data.")
  if (!is.element("DS_Trait",colnames(dat0))) stop("Error: There is no DS_Trait in the data.")
  
  dat = dat0[,c(factors,"TT_Trait","DS_Trait","Checks")]
  
  factors_list <<- sapply(dat[factors], function(x){sort(unique(x))}, simplify = FALSE)
  factors_length = lapply(factors_list,length)
  
  
  TT_Sigmab_nm = unlist(sapply(TT_raneff,function(r){
    if (r=="Variety") return("TT_Variety") else return(paste("TT_",r,"_",factors_list[[strsplit(r,split=":")[[1]][2]]],sep=""))
  },simplify = F),use.names = FALSE)
  DS_Sigmab_nm = unlist(sapply(DS_raneff,function(r){
    if (r=="Variety") return("DS_Variety") else return(paste("DS_",r,"_",factors_list[[strsplit(r,split=":")[[1]][2]]],sep=""))
  },simplify = F),use.names = FALSE)
  namesAllRan = c(TT_Sigmab_nm,DS_Sigmab_nm)
  
  ## split standard and non-standard varieties
  if (all(dat$Checks=="No")) {
    data = dat
  } else {
    datControlV <- dat[(dat$Variety %in% unique(dat$Variety[dat$Checks=="Yes"])), ]
    data <- dat[!(dat$Variety %in% unique(dat$Variety[dat$Checks=="Yes"])), ]
  }
  
  ## construct working balanced data
  AllLoc_Year = lapply(sort(unique(dat$Year)),function(x){return(sort(unique(dat$Loc)))}); names(AllLoc_Year) = sort(unique(dat$Year))
  
  ## non-standard varieties
  data$Variety = as.factor(data$Variety)
  dataVariety= split(data,data$Variety)
  split_data <- lapply(dataVariety, completeData, window=window, factors_list,AllLoc_Year)
  complete.data <- do.call(rbind, split_data)
  
  if (!is.null(datControlV)){
    ## standard varieties
    complete.data.control = datControlV %>% mutate(VarietybyYear=paste(Variety,Year,sep="-"), P_tilde=1,New=1*(!duplicated(Variety)),E=1, rep1 = 1, weights1=1,Checks=1)
    complete.data.control$New[which(complete.data.control$VarietybyYear %in% complete.data.control$VarietybyYear[which(!is.na(complete.data.control$New))])] <- 1
    ## Combine the non-standard and standard varieties
    complete.data.all = rbind(complete.data[,match(colnames(complete.data.control),colnames(complete.data))],complete.data.control)
  } else complete.data.all = complete.data
  
  ## Imputation for missing values at locations/replicates
  complete.data.all$VbyLbyYbyR <- paste(complete.data.all$Variety, complete.data.all$Loc, complete.data.all$Year, complete.data.all$Rep, sep = '-')
  im1dat <- complete.data.all[which(complete.data.all$P_tilde==1),]
  TT_im1dat_na <- im1dat[is.na(im1dat$TT_Trait), ];
  DS_im1dat_na <- im1dat[is.na(im1dat$DS_Trait), ];
  TT_imlm <- lm(TT_Trait ~ Year + Variety + Loc + Rep, data = im1dat)
  TT_imLoc <- predict(TT_imlm, newdata = TT_im1dat_na)
  complete.data.all[complete.data.all$VbyLbyYbyR %in% TT_im1dat_na$VbyLbyYbyR, 'TT_Trait'] <- TT_imLoc
  DS_imlm <- lm(DS_Trait ~ Year + Variety + Loc + Rep, data = im1dat)
  DS_imLoc <- predict(DS_imlm, newdata = DS_im1dat_na)
  complete.data.all[complete.data.all$VbyLbyYbyR %in% DS_im1dat_na$VbyLbyYbyR, 'DS_Trait'] <- DS_imLoc
  
  ## Select the post-entry data
  if (!yes_full) {complete.data.all = complete.data.all[complete.data.all$E==1,]}
  
  
  complete.data.all$Year <- factor(complete.data.all$Year, levels = sort(unique(as.character(complete.data.all$Year))))
  complete.data.all$Loc <- factor(complete.data.all$Loc, levels = sort(unique(as.character(complete.data.all$Loc))))
  
  split_complete_data <- split(complete.data.all, complete.data.all$Variety)
  naLast_TT_data <- lapply(split_complete_data, naLast,nm="TT_Trait")
  TT_Num <- sum(sapply(naLast_TT_data,nrow))
  naLast_DS_data <- lapply(split_complete_data,naLast,nm="DS_Trait")
  DS_Num <- sum(sapply(naLast_DS_data,nrow))
  
  TT_X <- lapply(naLast_TT_data, design_matrix_X, TT_fixeff) # design matrix of fixed effects for TT traits
  DS_X <- lapply(naLast_DS_data, design_matrix_X, DS_fixeff) # design matrix of fixed effects for DS trait
  TT_DS_X <- mapply(FUN=function(x1, x2){as.matrix(bdiag(x1, x2))}, TT_X, DS_X, SIMPLIFY=F) # Overall design matrix of fixed effects
  
  TT_Z <- lapply(naLast_TT_data, design_matrix_Z, TT_raneff) # deisgn matrix of random effects for
  
  DS_Z <- lapply(naLast_DS_data, design_matrix_Z, DS_raneff)# design matrix of random effects
  TT_DS_Z <- mapply(function(z1, z2){as.matrix(bdiag(z1, z2))}, TT_Z, DS_Z, SIMPLIFY = F) #
  
  complete.data.all[,factors] = apply(complete.data.all[,factors],2,function(x){factor(as.character(x))})
  
  fmla.fixeff = paste(TT_fixeff,collapse=" + ")
  fmla.raneff = paste(paste("(1|",TT_raneff,")",sep=""),collapse=" + ")
  TT_fmla = as.formula(paste("TT_Trait ~",fmla.fixeff,"+",fmla.raneff))
  TT_fitCC = lmer(TT_fmla,data=complete.data.all[!is.na(complete.data.all$TT_Trait),]) # complete case analysis for target traits
  
  fmla.fixeff = ifelse(is.null(DS_fixeff),"",paste(DS_fixeff,collapse="+"))
  fmla.raneff = paste(paste("(1|",DS_raneff,")",sep=""),collapse=" + ")
  DS_fmla = as.formula(paste("DS_Trait ~",fmla.fixeff,"+",fmla.raneff))
  DS_fitCC = lmer(DS_fmla, data=complete.data.all[!is.na(complete.data.all$DS_Trait),])# complete case analysis for decision-making traits
  
  ## initial parameter values for target trait
  TT_theta = summary(TT_fitCC)$coef[,1]; names(TT_theta) = paste("TT_",c("int",names(TT_theta)[-1]),sep="")
  TT_para_raneff = as.data.frame(VarCorr(TT_fitCC))[,'vcov']; names(TT_para_raneff) = as.data.frame(VarCorr(TT_fitCC))[,'grp']
  gamma_TT_Sigmab = TT_para_raneff[TT_raneff]; names(gamma_TT_Sigmab) = paste("TT_",TT_raneff,sep="")
  TT_Sigmab <<- unlist(lapply(TT_raneff,function(r){
    if (r=="Variety") return(TT_para_raneff[r]) else {
      r1 = strsplit(r,split=":")[[1]][2]
      out = rep(TT_para_raneff[r],factors_length[[r1]])
      names(out) = paste(r,factors_list[[r1]],sep="_")
      return(out)
    }
  }));
  names(TT_Sigmab) = paste("TT_",names(TT_Sigmab),sep="")
  TT_Sigmab = TT_Sigmab[match(TT_Sigmab_nm,names(TT_Sigmab))]
  ## initial parameter values for decision score
  DS_theta = summary(DS_fitCC)$coef[,1]; names(DS_theta) = paste("DS_",c("int",names(DS_theta)[-1]),sep="")
  DS_para_raneff = as.data.frame(VarCorr(DS_fitCC))[,'vcov']; names(DS_para_raneff) = as.data.frame(VarCorr(DS_fitCC))[,'grp']
  gamma_DS_Sigmab = DS_para_raneff[DS_raneff]; names(gamma_DS_Sigmab) = paste("DS_",DS_raneff,sep="")
  DS_Sigmab <<- unlist(lapply(DS_raneff,function(r){
    if (r=="Variety") return(DS_para_raneff[r]) else {
      r1 = strsplit(r,split=":")[[1]][2]
      out = rep(DS_para_raneff[r],factors_length[[r1]])
      names(out) = paste(r,factors_list[[r1]],sep="_")
      return(out)
    }
  }));
  names(DS_Sigmab) = paste("DS_",names(DS_Sigmab),sep="")
  DS_Sigmab = DS_Sigmab[match(DS_Sigmab_nm,names(DS_Sigmab))]
  ## Combine initial parameter values of target trait and decision score
  theta <<- c(TT_theta, DS_theta)
  Sigmab <<- c(TT_Sigmab, DS_Sigmab)
  rho <<- 0.2
  sigm2 <<- c(TT_para_raneff['Residual'], DS_para_raneff['Residual']); names(sigm2)=c("TT_sigm2","DS_sigm2") # error variance
  ## initial values values of logistic models and
  set.seed(2019)
  bet <<- runif(factors_length$Year)
  
  
  gamma_now <- list(theta=theta, Sigmab=c(gamma_TT_Sigmab,gamma_DS_Sigmab), sigm2=sigm2, bet=bet, rho = rho) # all the parameters
  XtX <- mapply(function(x){t(x) %*% x}, TT_DS_X, SIMPLIFY=FALSE) ## Xi'Xi
  sumXtX <- apply(array(unlist(XtX), dim = c(dim(XtX[[1]])[1], dim(XtX[[1]])[1], factors_length$Variety)), c(1,2), sum) ## sum_i Xi'Xi
  TT_Tobs <- lapply(naLast_TT_data, function(x){matrix(rep(x$TT_Trait[!is.na(x$TT_Trait)], nsamp), nrow=sum(!is.na(x$TT_Trait)), ncol=nsamp, byrow=F)}) ## create a matrix with the observed value of target trait repeated nsamp columns
  DS_Tobs <- lapply(naLast_DS_data, function(x){matrix(rep(x$DS_Trait[!is.na(x$DS_Trait)], nsamp), nrow=sum(!is.na(x$DS_Trait)), ncol=nsamp, byrow=F)}) ## create a matrix with the observed value of decision-making trait repeated nsamp columns
  
  TT_obs <- lapply(naLast_TT_data, function(x){x$TT_Trait[!is.na(x$TT_Trait)]})
  DS_obs <- lapply(naLast_DS_data, function(x){x$DS_Trait[!is.na(x$DS_Trait)]})
  TT_DS_obs <- mapply(function(tt, DS){c(tt, DS)}, TT_obs, DS_obs)
  TT_NumM <- lapply(naLast_TT_data, function(x){sum(is.na(x$TT_Trait))}) # number of missing values of target trait for each variety
  DS_NumM <- lapply(naLast_DS_data, function(x){sum(is.na(x$DS_Trait))})# number of missing values of decision making trait for each variety
  # set.seed(2019)
  TT_init <- mapply(function(obs, m){
    # if there is no missing data, return NULL; otherwise generate the inital value of the missing data
    if(m==0) return(NULL) else return(rnorm(m,mean(obs)-1))
  }, TT_obs, TT_NumM, SIMPLIFY=FALSE)
  DS_init <- mapply(function(obs, m){
    # if there is no missing data, return NULL; otherwise, generate the inital value of the missing data
    if(m==0)return(NULL)else return(rnorm(m,mean(obs)-1))
  }, DS_obs, DS_NumM, SIMPLIFY=FALSE)
  TT_DS_init <- mapply(function(tt, DS){c(tt, DS)}, TT_init, DS_init, SIMPLIFY = F)
  
  DS_mu <- mean(unlist(DS_obs))##### average of all the observed value of decison-making trait, which is higher than the true value
  DS_st <- sd(unlist(DS_obs)) ## standard deviation of all the observed value of decison-making trait
  
  iter <- 1
  change <- TRUE # the change of values of parameters of main interest
  change1 <- TRUE # the change of values of nuisance parameters
  gamma_path = unlist(gamma_now)
  
  t0 = Sys.time()
  while(change & iter < max.iter){
    cat("Iteration", iter,"\n")
    
    Sigmab_all = mapply(SigmaB,naLast_TT_data,naLast_DS_data,SIMPLIFY=F)
    sigma2_err_all = mapply(function(TT_dat,DS_dat){return(diag(rep(sigm2, c(nrow(TT_dat),nrow(DS_dat)))))},naLast_TT_data, naLast_DS_data, SIMPLIFY=F)
    
    # E-step: sampling the missing value of target trait and decision score
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
    
    ## Update  beta
    Mmatt <- mapply(glmMavg,naLast_DS_data, DS_u,SIMPLIFY=FALSE)
    MmatDat <- do.call(rbind.data.frame, Mmatt)
    MmatDat <- MmatDat[MmatDat$New == 0 & MmatDat$E == 1 & MmatDat$Checks == 0,,drop=F]
    MmatDat$TraitAvg <- (MmatDat$TraitAvg - DS_mu) / DS_st # standardize
    MmatDat$TraitAvgPre <- ifelse(MmatDat$TraitAvgPre != 0, (MmatDat$TraitAvgPre - DS_mu) / DS_st, MmatDat$TraitAvgPre)
    MmatDat$Year = factor(MmatDat$Year, levels=factors_list$Year[-1])
    bet_new = glm(P_tilde ~ Year + TraitAvg, data = MmatDat, family = binomial)$coefficients
    
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
      TT_DS_Ttrait,TT_DS_X,TT_DS_Z,Sigmab_all,sigma2_err_all,SIMPLIFY = F)
    
    Sigmab_i <- mapply(b, TT_DS_Z, Sigmab_all,sigma2_err_all,b1_all, SIMPLIFY=FALSE)
    D_pre = lapply(Sigmab_i,function(x){diag(x)[match(namesAllRan,names(diag(x)))]})
    Sigmab_new <- apply(do.call(cbind,D_pre),1,mean,na.rm=T) # diagnoal elements of Sigma_b
    names(Sigmab_new) = namesAllRan
    if (sum(TT_raneff!="Variety")>0){
      for (r in TT_raneff[TT_raneff!="Variety"]){r = paste("TT_",r,sep=""); Sigmab_new[grep(r,namesAllRan)] = rep(mean(Sigmab_new[grep(r,namesAllRan)]),length(grep(r,namesAllRan)))
      }
    }
    if (sum(DS_raneff!="Variety")>0){
      for (r in DS_raneff[DS_raneff!="Variety"]){
        r = paste("DS_",r,sep="")
        Sigmab_new[grep(r,namesAllRan)] = rep(mean(Sigmab_new[grep(r,namesAllRan)]),length(grep(r,namesAllRan)))
      }
    }
    
    rho_update <- mean(sapply(Sigmab_i, function(tt){tt['TT_Variety', 'DS_Variety']})) # covaraince of TT_Variety and DS_Variety
    rho_new = rho_update/sqrt(Sigmab_new['TT_Variety']*Sigmab_new['DS_Variety'])
    names(rho_new) = NULL
    
    ## update theta
    theta_n <- mapply(function(Ttrait,x,z,b1){
      t(x) %*% apply(Ttrait-z %*% b1, 1, mean)
    },TT_DS_Ttrait,TT_DS_X,TT_DS_Z,b1_all,SIMPLIFY = FALSE)
    theta_new = c(solve(sumXtX) %*% apply(array(unlist(theta_n),
                                                dim = c(length(theta), 1, length(theta_n))), c(1,2), sum))
    names(theta_new) <- names(theta)
    
    ## update error variance
    term1 <- mapply(error, TT_DS_Ttrait, TT_DS_X, TT_DS_Z, naLast_TT_data, Sigmab_all,sigma2_err_all,b1_all, SIMPLIFY=FALSE)
    term2 <- mapply(traceZD, TT_DS_Z, naLast_TT_data, Sigmab_all,sigma2_err_all, SIMPLIFY=FALSE)
    sigm2_new = apply(mapply(function(x,y){x+y}, term1, term2,SIMPLIFY=TRUE),1,sum)/c(TT_Num, DS_Num)
    names(sigm2_new) = names(sigm2)
    
    gamma_Sigmab_TT = sapply(paste("TT",TT_raneff,sep="_"),function(r){
      if (r=="TT_Variety") return(Sigmab_new[r]) else return(unique(Sigmab_new[grep(r,names(Sigmab_new))]))
    }); names(gamma_Sigmab_TT) = paste("TT",TT_raneff,sep="_")
    gamma_Sigmab_DS = sapply(paste("DS",DS_raneff,sep="_"),function(r){
      if (r=="DS_Variety") return(Sigmab_new[r]) else return(unique(Sigmab_new[grep(r,names(Sigmab_new))]))
    }); names(gamma_Sigmab_DS) = paste("DS",DS_raneff,sep="_")
    
    gamma_new <- list(theta=theta_new, Sigmab=c(gamma_Sigmab_TT,gamma_Sigmab_DS), sigm2=sigm2_new, bet=bet_new, rho = rho_new) # all the parameters
    gamma_path = cbind(gamma_path,unlist(gamma_new))
    
    
    # convergence
    gamma_new_v = unlist(gamma_new); gamma_now_v = unlist(gamma_now)
    
    dif0 <- max(abs(gamma_new_v[grep('TT',names(gamma_new_v))]-gamma_now_v[grep('TT',names(gamma_now_v))]))
    print(dif0)
    change0 <- dif0 > err
    
    dif1 <- max(abs(gamma_new_v[-grep('TT',names(gamma_new_v))]-gamma_now_v[-grep('TT',names(gamma_now_v))]))
    print(dif1)
    change1 <- dif1 > err1
    
    change = change0 | change1
    
    assign("theta", theta_new, envir = .GlobalEnv)
    assign("Sigmab", Sigmab_new, envir = .GlobalEnv)
    assign("TT_Sigmab",Sigmab_new[match(names(TT_Sigmab),names(Sigmab_new))], envir = .GlobalEnv)
    assign("DS_Sigmab",Sigmab_new[match(names(DS_Sigmab),names(Sigmab_new))], envir = .GlobalEnv)
    assign("sigm2", sigm2_new, envir = .GlobalEnv)
    assign("rho", rho_new, envir = .GlobalEnv)
    assign("bet", bet_new, envir = .GlobalEnv)
    gamma_now <- gamma_new
    print(gamma_now)
    
    TT_DS_init <- lapply(u, function(x){
      if (length(x) == 0){ # missing value
        return(x)
      } else {
        return(apply(x, 2, mean))
      }
    })
    
    iter <- iter+1
  }
  run_time = Sys.time() - t0
  print(run_time)
  
  return(list(para=list(theta=theta,Sigmab=Sigmab,TT_Sigmab=TT_Sigmab,DS_Sigmab=DS_Sigmab,sigm2=sigm2,rho=rho,bet=bet),compdata=list(TT=naLast_TT_data,DS=naLast_DS_data,X=TT_DS_X,Z=TT_DS_Z,init=TT_DS_init,obs=TT_DS_obs,DS_mu=DS_mu,DS_st=DS_st,TT_Tobs=TT_Tobs,DS_Tobs=DS_Tobs),factors_list=factors_list,model=list(TT_mm=TT_mm,DS_mm=DS_mm,window=window,yes_full=yes_full,converg_control=converg_control),gamma_path=gamma_path,run_time=run_time))
}