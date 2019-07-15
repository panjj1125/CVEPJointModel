naLast <- function(x,nm){x[c(which(is.na(x[,nm]) == FALSE),which(is.na(x[,nm]) == TRUE)),]} 

design_matrix_X <- function(tmpdat, fixeffect){
  if (is.null(fixeffect)) {
    X <- model.matrix(~ 1, tmpdat)
  } else X <- model.matrix(reformulate(termlabels = fixeffect), tmpdat)  
  return(X=X)
}

design_matrix_Z <- function(tmpdat, raneffect){
  if (length(raneffect)==1 & is.element("Variety",raneffect)){
    ## random effects include only variety effect
    Z <- model.matrix(~ 1, tmpdat)
    } else{
    Z1 <- lapply(raneffect[raneffect!="Variety"],function(r1){
      r = strsplit(r1,split=":")[[1]][2]; tmpdat[,r] = as.character(tmpdat[,r])
      if(nrow(unique(tmpdat[r]))==1){
        out = matrix(rep(1,nrow(tmpdat)),ncol = 1); colnames(out) = paste(r1,unique(tmpdat[r]),sep="")
        return(out)
      } else {
        out = model.matrix(reformulate(termlabels = c(r,'-1')), tmpdat)
        colnames(out) = paste("Variety:",colnames(out),sep="")
        return(out)
      }
    }
    )
    Z <- cbind(model.matrix(~ 1, tmpdat), as.matrix(do.call(cbind.data.frame, Z1)))
  }
  return(Z=Z)
}



SigmaB <- function(TT_dat,DS_dat){
  TT_sigma_nm = unlist(lapply(TT_raneff,function(r){
    if (r=="Variety") return("TT_Variety") else {
      r1 = strsplit(r,split=":")[[1]][2]
      return(paste("TT",r,unique(as.character(TT_dat[,r1])),sep="_"))
    }
  }))
  DS_sigma_nm = unlist(lapply(DS_raneff,function(r){
    if (r=="Variety") return("DS_Variety") else {
      r1 = strsplit(r,split=":")[[1]][2]
      return(paste("DS",r,unique(as.character(DS_dat[,r1])),sep="_"))
    }
  }))
  TT_DS_sigma <- c(Sigmab[c(TT_sigma_nm,DS_sigma_nm)])
  Sigmab_TT_DS <- diag(c(TT_DS_sigma), length(TT_DS_sigma)); rownames(Sigmab_TT_DS) = colnames(Sigmab_TT_DS) = names(TT_DS_sigma)
  Sigmab_TT_DS['TT_Variety', 'DS_Variety'] <- Sigmab_TT_DS['DS_Variety', 'TT_Variety'] <- rho * sqrt(TT_DS_sigma["TT_Variety"]) * sqrt(TT_DS_sigma["DS_Variety"])
  return(Sigmab_TT_DS)
}

slice_samp <- function(init, obs, TT_dat, DS_dat, x, z, DS_mu, DS_st,Sigmab_TT_DS){
  if (length(init) == 0) {
    ## if no missing data
    return(NULL) 
  } else {
    ttl <- nrow(TT_dat)
    DSl <- nrow(DS_dat)
    ## marginal mean and covariance matrix of both observed and missing
    u_mu <- x%*%as.matrix(theta)
    u_Sigma <- z%*%Sigmab_TT_DS%*%t(z) + diag(rep(sigm2, c(ttl, DSl)))
    ## conditional mean and covariance matrix of missing given observed
    l1 <- sum(is.na(DS_dat$DS_Trait)) # number of missing DS data
    l2 <- sum(!is.na(TT_dat$TT_Trait)) # number of observed TT data
    ind1 <- (l2+1):ttl # index for the missing target trait values 
    ind2 <- (ttl + DSl - l1 + 1):(ttl+DSl)# index for the missing decison-making trait values
    na_ind <- c(ind1, ind2) # missing data index
    u_obs_mu <- u_mu[-na_ind] # mean for observed data
    u_m_obs <- u_mu[na_ind] # mean for missing data
    Sigma11 <- u_Sigma[-na_ind, -na_ind]
    Sigma22 <- u_Sigma[na_ind, na_ind]
    Sigma21 <- u_Sigma[na_ind, -na_ind]
    Sigma12 <- u_Sigma[-na_ind, na_ind]
    mu_miss <- u_m_obs + Sigma21%*%solve(Sigma11)%*%(obs-u_obs_mu)
    Sigma_miss <- Sigma22 - Sigma21%*%solve(Sigma11)%*%Sigma12
    
    DS_dat$Year = as.character(DS_dat$Year)
    DS_dat = DS_dat[order(DS_dat$Year),]
    Ptilde = tapply(DS_dat$P_tilde,DS_dat$Year,unique)
    Ptilde = Ptilde[order(names(Ptilde))]
    E = tapply(DS_dat$E,DS_dat$Year,unique)
    E = E[order(names(E))]
    ## sampling
    sss <- matrix(NA, length(init), nsamp)
    i <- 1
    while(i<=nsamp){
      ss <- rmvnorm(1, mu_miss, Sigma_miss)
      # calculate the acceptance indicators for all the missing values
      if (sum(E==1 & Ptilde==0)>0){
        DS_trait = DS_dat$DS_Trait
        DS_trait[is.na(DS_trait)] = ss[-c(1:length(ind1))]
        DS_TraitComp <- tapply(DS_trait,DS_dat$Year,mean)
        DS_TraitComp <- DS_TraitComp[order(names(DS_TraitComp))]
        DS_TraitComp <- (DS_TraitComp - DS_mu)/DS_st
        m1_dat <- data.frame(cbind(Year=names(DS_TraitComp), TraitAvg=DS_TraitComp, TraitAvgPre=c(0,DS_TraitComp[-length(DS_TraitComp)])))[E==1 & Ptilde==0,,drop=F]
        m1_dat$Year = factor(m1_dat$Year,levels=factors_list$Year[-1])
        m1_dat$TraitAvg = as.numeric(as.character(m1_dat$TraitAvg))
        m1_dat$TraitAvgPre = as.numeric(as.character(m1_dat$TraitAvgPre))
        m1 = model.matrix(reformulate(termlabels = c("TraitAvg","TraitAvgPre")), m1_dat)
        value <- c(exp(m1%*%as.matrix(bet))/(1+exp(m1%*%as.matrix(bet))))
        plant <- 1-rbinom(length(value), 1, value)
        if(!is.na(prod(plant)) & prod(plant) == 1){ ## accept the new value
          sss[,i] <- ss
          i <- i+1
        } else{i <- i}
      } else {
        sss[,i] <- ss
        i <- i+1
      }
    }
  }
  return(t(sss))
}


glmMavg <- function(DS_dat,u1){
  DS_dat$Year = as.character(DS_dat$Year)
  DS_dat = DS_dat[order(DS_dat$Year),]
  DS_dat_year = DS_dat %>% group_by(Year) %>% 
    summarise(Variety=unique(Variety), P_tilde=unique(P_tilde), New=unique(New),E=unique(E), Checks=unique(Checks),rep1=unique(rep1), TraitAvg=mean(DS_Trait))
  YearCur = as.numeric(DS_dat_year$Year); YearPre = YearCur-1;
  
  if(!is.null(u1)){# no missing value
    u1 = as.matrix(u1, length(u1))
    colnames(u1) = DS_dat$Year[DS_dat$rep1==nsamp]
    u_avg = apply(matrix(apply(u1,1,function(x){tapply(x,colnames(u1),mean)}),byrow=T,nrow=nsamp),2,mean); names(u_avg) = unique(colnames(u1))
    DS_dat_year$TraitAvg[is.na(DS_dat_year$TraitAvg)==T] = u_avg[match(DS_dat_year$Year[is.na(DS_dat_year$TraitAvg)==T],names(u_avg))]
  }
  DS_dat_year$New[!is.element(YearPre,YearCur)] = 1
  DS_dat_year$TraitAvgPre = c(0,DS_dat_year$TraitAvg[-nrow(DS_dat_year)])
  return(DS_dat_year)
}


b <- function(z, Sigmab_TT_DS,sigma2_err,b1){ 
  d <- Sigmab_TT_DS - Sigmab_TT_DS %*% t(z) %*% solve(z %*% Sigmab_TT_DS %*% t(z) + sigma2_err) %*% z %*% Sigmab_TT_DS # conditional covariance matrix of b given observed target trait  
  b2 <- as.matrix(b1 %*% t(b1)/nsamp + d)
  return(b2)
}
theta_sol <- function(Ttrait, x, z, TT_dat, Sigmab_TT_DS){ #solve theta
  # this function calculates x'avg(ti-zibi)
  sigma2_err <- diag(rep(sigm2, c(dim(TT_dat)[1], dim(x)[1]-dim(TT_dat)[1])))
  b1 <- matrix(apply(Ttrait, 2, function(t1){Sigmab_TT_DS %*% t(z) %*% solve(z %*% Sigmab_TT_DS %*% t(z) + sigma2_err) %*% (t1-x %*% theta)}),nrow=dim(z)[2]); names(b1) = colnames(Sigmab_TT_DS) # bi: conditional mean of b given observed target trait 
  the <- t(x) %*% apply(Ttrait-z %*% b1, 1, mean) 
  return(the)
}

error <- function(Ttrait, x, z, TT_dat, Sigmab_TT_DS,sigma2_err,b1){ #solve sigma
  TT_DS_err <- (Ttrait-c(x %*% theta) - z %*% b1)
  TT_err <- TT_DS_err[1:dim(TT_dat)[1],]
  DS_err <- TT_DS_err[-(1:dim(TT_dat)[1]),]
  TT_errm <- mean(diag(t(TT_err) %*% (TT_err)))
  DS_errm <- mean(diag(t(DS_err) %*% (DS_err)))
  err <- c(TT_errm, DS_errm)
  return(err)
}
traceZD <- function(z,TT_dat, Sigmab_TT_DS,sigma2_err){ #solve sigma
  d <- Sigmab_TT_DS - Sigmab_TT_DS %*% t(z) %*% solve(z %*% Sigmab_TT_DS %*% t(z) + sigma2_err) %*% z %*% Sigmab_TT_DS # conditional covariance matrix of b given observed target trait  
  TT2 <- sum(diag(t(z) %*% z %*% d)[1:(dim(z)[2]-1)])
  DS2 <- sum(diag(t(z) %*% z %*% d)[-c(1:(dim(z)[2]-1))])
  return(c(TT2,DS2))
}


completeData <- function(dataV, window=3,factors_list,AllLoc_Year){
  dataV$Year <- as.character(dataV$Year)
  dataV$Variety <- as.character(dataV$Variety)
  YearAll <- as.numeric(as.character(sort(unique(dataV$Year))))
  if (is.infinite(window)){
    # Yearfull = factors_list$Year
    YearFull = as.numeric(factors_list$Year)[as.numeric(factors_list$Year)>=min(YearAll)]
  } else {
    if(length(YearAll) == 1){
      YearL <- seq(YearAll, YearAll+window-1, 1)
    } else{
      YearL <- sapply(YearAll, function(x){seq(x,x+window-1,1)}, simplify = FALSE)
      i <- 2
      while (i <= length(YearL)) {
        if ((YearL[[i]][1] - YearL[[i-1]][1]) < window) {
          YearL <- YearL[-i]
        } else {i <- i+1}
      }
    }
    YearFull <- unlist(YearL)
  }
  YearFull = YearFull[which(YearFull %in% factors_list$Year)] # new: some year is after cycle 17
  
  New <- rep(NA, nrow(dataV))
  New[!duplicated(dataV$Variety)] <- 1
  P_tilde <- rep(1,nrow(dataV))
  E <- rep(1,nrow(dataV))
  data1 <- cbind(dataV, P_tilde=P_tilde, New=New, E=E)
  
  # creat complete data
  out = lapply(YearFull,function(yr){
    loc_list = AllLoc_Year[[which(names(AllLoc_Year)==as.character(yr))]]
    rep_list = factors_list$Rep
    return(cbind(Year=yr,Loc=rep(loc_list,each=length(rep_list)),Rep=rep(rep_list,length(loc_list))))
  })
  
  out = do.call(rbind,out)
  compdata <- data.frame(Variety=rep(unique(data1$Variety),nrow(out)),out)
  complete_data <- left_join(compdata, data1)
  
  if (!is.element("VarietybyYear",colnames(complete_data))) complete_data = complete_data %>% mutate(VarietybyYear=paste(Variety,Year,sep="-"))
  
  complete_data$New[which(complete_data$VarietybyYear %in% complete_data$VarietybyYear[which(!is.na(complete_data$New))])] <- 1
  complete_data$New[is.na(complete_data$New)] <- 0
  
  complete_data$E[which(complete_data$Year %in% YearFull[c(which(YearFull==YearAll[1]):length(YearFull))])] = 1
  complete_data$E[is.na(complete_data$E)] = 0
  complete_data$P_tilde[which(complete_data$Year %in% YearAll)] = 1
  complete_data$P_tilde[is.na(complete_data$P_tilde)] = 0
  complete_data$rep1 <- ifelse(complete_data$P_tilde == 1, 1, nsamp)
  complete_data$weights1 <- ifelse(complete_data$P_tilde == 1, 1, 1/nsamp)
  complete_data$Checks <- 0
  return(complete_data)
}