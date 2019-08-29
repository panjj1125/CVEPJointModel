example <- function(dat0){
  
  dat0_TT_mu = mean(dat0$TT_Trait,na.rm=T)
  dat0_TT_sd = sd(dat0$TT_Trait,na.rm=T)
  dat0$TT_Trait = scale(dat0$TT_Trait)
  dat0$DS_Trait = scale(log(dat0$DS_Trait))
  
  source("FUNS.R")
  source("CVEP_JM.R")
  res_PC_post = CVEP_JM(MET_dat=dat0,factors= c("Variety","Year","Loc","Rep"),TT_mm=list(fixed=c("Year","Loc"),random=c("Variety","Variety:Loc")),DS_mm=list(fixed=c("Year","Loc"),random=c("Variety")),window=Inf,yes_full=F,converg_control=list(nsamp=500,max.iter=1000,err=5*10^(-5),err1=5*10^(-4),seed=20190421))
  res_PC_full = CVEP_JM(MET_dat=dat0,factors= c("Variety","Year","Loc","Rep"),TT_mm=list(fixed=c("Year","Loc"),random=c("Variety","Variety:Loc")),DS_mm=list(fixed=c("Year","Loc"),random=c("Variety")),window=Inf,yes_full=T,converg_control=list(nsamp=500,max.iter=1000,err=5*10^(-5),err1=5*10^(-4),seed=20190421))
  source("CVEP_JM_Year.R")
  res_YC_post = CVEP_JM_Year(MET_dat=dat0,factors= c("Variety","Year","Loc","Rep"),TT_mm=list(fixed=c("Year","Loc"),random=c("Variety","Variety:Loc")),DS_mm=list(fixed=c("Year","Loc"),random=c("Variety")),window=Inf,yes_full=F,converg_control=list(nsamp=500,max.iter=1000,err=5*10^(-5),err1=5*10^(-4),seed=20190421))
  res_YC_full = CVEP_JM_Year(MET_dat=dat0,factors= c("Variety","Year","Loc","Rep"),TT_mm=list(fixed=c("Year","Loc"),random=c("Variety","Variety:Loc")),DS_mm=list(fixed=c("Year","Loc"),random=c("Variety")),window=Inf,yes_full=T,converg_control=list(nsamp=500,max.iter=1000,err=5*10^(-5),err1=5*10^(-4),seed=20190421))
  
  save(list=ls(),file="res.RData")
  # load("res_No_Year.RData")
  
  para_est_all = NULL
  ## Complete Case Analysis
  dat0$Loc = as.character(dat0$Loc)
  dat0$Year = as.character(dat0$Year)
  dat0$Variety = as.character(dat0$Variety)
  year_list = sort(unique(dat0$Year))
  loc_list = sort(unique(dat0$Loc))
  para_fixed_nm = c("int",paste("Year",year_list[-1],sep=""),paste("Loc",loc_list[-1],sep=""))
  library(lme4)
  fit_m = lmer(TT_Trait ~ Year + Loc + (1|Variety) + (1|Variety:Loc),data=dat0)
  theta_est = summary(fit_m)$coef[,1]; names(theta_est) = c("int",paste("Year",year_list[-1],sep=""),paste("Loc",loc_list[-1],sep=""))
  theta_est = theta_est*dat0_TT_sd; theta_est[1] = theta_est[1]+dat0_TT_mu
  Sigmab = as.data.frame(VarCorr(fit_m))
  Sigmab_est = Sigmab[match(c('Variety',"Variety:Loc"),Sigmab$grp),'vcov']
  names(Sigmab_est) = c("Variety","Variety:Loc")
  sigma2_est = Sigmab[match(c("Residual"),Sigmab$grp),'vcov']
  para_est_all = cbind(para_est_all,c(theta_est,sqrt(Sigmab_est)*dat0_TT_sd,sigma=sqrt(sigma2_est)*dat0_TT_sd))
  ## PC
  theta_est = res_PC_full$para$theta
  theta_est = theta_est[paste("TT",para_fixed_nm,sep="_")]
  theta_est = theta_est*dat0_TT_sd; theta_est[1] = theta_est[1]+dat0_TT_mu
  Sigmab_est = c(res_PC_full$para$Sigmab["TT_Variety"],unique(res_PC_full$para$Sigmab[grep("TT_Variety:Loc",names(res_PC_full$para$Sigmab))]))
  sigma2_est = res_PC_full$para$sigm2["TT_sigm2"]
  para_est_all = cbind(para_est_all,c(theta_est,sqrt(Sigmab_est)*dat0_TT_sd,sigma=sqrt(sigma2_est)*dat0_TT_sd))
  theta_est = res_PC_post$para$theta
  theta_est = theta_est[paste("TT",para_fixed_nm,sep="_")]
  theta_est = theta_est*dat0_TT_sd; theta_est[1] = theta_est[1]+dat0_TT_mu
  Sigmab_est = c(res_PC_post$para$Sigmab["TT_Variety"],unique(res_PC_post$para$Sigmab[grep("TT_Variety:Loc",names(res_PC_post$para$Sigmab))]))
  sigma2_est = res_PC_post$para$sigm2["TT_sigm2"]
  para_est_all = cbind(para_est_all,c(theta_est,sqrt(Sigmab_est)*dat0_TT_sd,sigma=sqrt(sigma2_est)*dat0_TT_sd))
  ## YC
  theta_est = res_YC_full$para$theta
  theta_est = theta_est[paste("TT",para_fixed_nm,sep="_")]
  theta_est = theta_est*dat0_TT_sd; theta_est[1] = theta_est[1]+dat0_TT_mu
  Sigmab_est = c(res_YC_full$para$Sigmab["TT_Variety"],unique(res_YC_full$para$Sigmab[grep("TT_Variety:Loc",names(res_YC_full$para$Sigmab))]))
  sigma2_est = res_YC_full$para$sigm2["TT_sigm2"]
  para_est_all = cbind(para_est_all,c(theta_est,sqrt(Sigmab_est)*dat0_TT_sd,sigma=sqrt(sigma2_est)*dat0_TT_sd))
  theta_est = res_YC_post$para$theta
  theta_est = theta_est[paste("TT",para_fixed_nm,sep="_")]
  theta_est = theta_est*dat0_TT_sd; theta_est[1] = theta_est[1]+dat0_TT_mu
  Sigmab_est = c(res_YC_post$para$Sigmab["TT_Variety"],unique(res_YC_post$para$Sigmab[grep("TT_Variety:Loc",names(res_YC_post$para$Sigmab))]))
  sigma2_est = res_YC_post$para$sigm2["TT_sigm2"]
  para_est_all = cbind(para_est_all,c(theta_est,sqrt(Sigmab_est)*dat0_TT_sd,sigma=sqrt(sigma2_est)*dat0_TT_sd))
  output = cbind(format(round(para_est_all[,1],2),nsmall=2,trim=T),
                 paste(format(round(para_est_all[,2],2),nsmall=2,trim=T)," (",format(round(100*(para_est_all[,2]/para_est_all[,1]-1),1),nsmall=1,trim=T),"%)",sep=""),
                 paste(format(round(para_est_all[,3],2),nsmall=2,trim=T)," (",format(round(100*(para_est_all[,3]/para_est_all[,1]-1),1),nsmall=1,trim=T),"%)",sep=""),
                 paste(format(round(para_est_all[,4],2),nsmall=2,trim=T)," (",format(round(100*(para_est_all[,4]/para_est_all[,1]-1),1),nsmall=1,trim=T),"%)",sep=""),
                 paste(format(round(para_est_all[,5],2),nsmall=2,trim=T)," (",format(round(100*(para_est_all[,5]/para_est_all[,1]-1),1),nsmall=1,trim=T),"%)",sep="")
  )
  rownames(output) = c("Intercept",paste("Y",year_list[-1],sep=""),loc_list[-1],"Variety","VxL","Error")
  colnames(output) = c("CC","JM:P+C full", "JM:P+C Post-entry","JM:Y+C full", "JM:Y+C Post-entry")
  write.csv(output,"Table1.csv",row.names = T)
  
  ## DS Trait
  para_est_all = NULL
  theta_est = res_PC_full$para$theta
  theta_est = theta_est[paste("DS",para_fixed_nm,sep="_")]
  Sigmab_est = c(res_PC_full$para$Sigmab["DS_Variety"],unique(res_PC_full$para$Sigmab[grep("DS_Variety:Loc",names(res_PC_full$para$Sigmab))]))
  sigma2_est = res_PC_full$para$sigm2["DS_sigm2"]
  rho_est = res_PC_full$para$rho
  para_est_all = cbind(para_est_all,c(theta_est,sqrt(Sigmab_est),sigma=sqrt(sigma2_est),rho=rho_est))
  
  theta_est = res_PC_post$para$theta
  theta_est = theta_est[paste("DS",para_fixed_nm,sep="_")]
  Sigmab_est = c(res_PC_post$para$Sigmab["DS_Variety"],unique(res_PC_post$para$Sigmab[grep("DS_Variety:Loc",names(res_PC_post$para$Sigmab))]))
  sigma2_est = res_PC_post$para$sigm2["DS_sigm2"]
  rho_est = res_PC_post$para$rho
  para_est_all = cbind(para_est_all,c(theta_est,sqrt(Sigmab_est),sigma=sqrt(sigma2_est),rho=rho_est))
  
  theta_est = res_YC_full$para$theta
  theta_est = theta_est[paste("DS",para_fixed_nm,sep="_")]
  Sigmab_est = c(res_YC_full$para$Sigmab["DS_Variety"],unique(res_YC_full$para$Sigmab[grep("DS_Variety:Loc",names(res_YC_full$para$Sigmab))]))
  sigma2_est = res_YC_full$para$sigm2["DS_sigm2"]
  rho_est = res_YC_full$para$rho
  para_est_all = cbind(para_est_all,c(theta_est,sqrt(Sigmab_est),sigma=sqrt(sigma2_est),rho=rho_est))
  
  theta_est = res_YC_post$para$theta
  theta_est = theta_est[paste("DS",para_fixed_nm,sep="_")]
  Sigmab_est = c(res_YC_post$para$Sigmab["DS_Variety"],unique(res_YC_post$para$Sigmab[grep("DS_Variety:Loc",names(res_YC_post$para$Sigmab))]))
  sigma2_est = res_YC_post$para$sigm2["DS_sigm2"]
  rho_est = res_YC_post$para$rho
  para_est_all = cbind(para_est_all,c(theta_est,sqrt(Sigmab_est),sigma=sqrt(sigma2_est),rho=rho_est))

  output = format(round(para_est_all,3),nsmall=3)

  rownames(output) = c("Intercept",paste("Y",year_list[-1],sep=""),loc_list[-1],"Variety","Error","Corr.")
  colnames(output) = c("JM:P+C full", "JM:P+C Post-entry","JM:Y+C full", "JM:Y+C Post-entry")
  write.csv(output,"Supp_TableB1.csv",row.names = T)
  
  
  ## GLM
  rownm = c("(Intercept)","Year1966","TraitAvgPre","TraitAvg")
  tab_PC = cbind(res_PC_full$para$bet,res_PC_post$para$bet); tab_PC = tab_PC[match(rownm,rownames(tab_PC)),]; tab_PC = format(round(tab_PC,3),nsmall=3); tab_PC[is.na(rownames(tab_PC)),] = ""
  tab_YC = cbind(res_YC_full$para$bet,res_YC_post$para$bet); tab_YC = tab_YC[match(rownm,rownames(tab_YC)),]; tab_YC = format(round(tab_YC,3),nsmall=3); tab_YC[is.na(rownames(tab_YC)),] = ""
  output = cbind(tab_PC,tab_YC)
  rownames(output) = c("Intercept","Y1966","Previous","Current")
  colnames(output) = c("JM:P+C full", "JM:P+C Post-entry","JM:Y+C full", "JM:Y+C Post-entry")
  write.csv(output,"Supp_TableB2.csv",row.names = T)
}