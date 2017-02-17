require(plyr)
#the model: Pass it a file of subject biomarkers (including age), 
#it will make a datafile with subject IDs as rownames, Age, first bioage, and FinalBA 
#as the final version using chronological age adjustment. 
#Note the function does its own inverting of the appropriate variables, 
#so it is important to put in unalterred data

# set sex to "M" or "F" to generate appropriate graph
# specify whether function is being run on a "training" or "longitudinal" set to dictate output directory

BAmodel<-function(datBM,pars, sex, set){
  if(lsinvert[1] != "z"){datBM<-data.frame(cbind(datBM$Age,as.matrix(scale(datBM[-1]))%*%rot))
  colnames(datBM)[1]<-"Age"}
  datBM[,lsinvert[-1]]<- datBM[,lsinvert[-1]]*(-1)
  BioAge<-data.frame(cbind(datBM$Age,adply(datBM[-1], 1, function(x){
    su<- sum(((x[1,]-pars[2,])*pars[1,])/(pars[3,]^2))/sum((pars[1,]^2)/(pars[3,]^2))
    out<- data.frame(su)
    return(out)
  })$su), row.names=row.names(datBM))
  
  colnames(BioAge)<-c('Age','bioage')
  #BioAge now contains rows by subject ID, with age and bioage
  
  #characteristic r of the model (0 to 1, bigger is better)
  rchar =sum(pars[5,]*(pars[1,])/(pars[3,]))/sum((pars[1,])/(pars[3,]))
  
  pdf(paste0("output/", set, "Set/PrelimBA_", sex, ".pdf"))
  plot(BioAge$Age, BioAge$bioage)
  abline(lm(BioAge$Age ~ BioAge$bioage))
  dev.off()
  
  Sbsqn<-(sum(((BioAge[,"bioage"]-BioAge[,"Age"])-sum(BioAge[,"bioage"]-BioAge[,"Age"])/length(BioAge$Age))^2)/length(BioAge$Age))-(((max(BioAge$Age)-min(BioAge$Age))^2/(12*length(datBM[-1])))*(1-(rchar^2))/(rchar^2))
  
  Sbeb<-((1-rchar^2)^.5)*(max(BioAge$Age)-min(BioAge$Age))/(rchar*((12*length(datBM[-1]))^.5))
  
  Sbecb<-Sbeb/(1+((Sbeb/Sbsqn)^2))
  assign("stErrorinYears",Sbecb,.GlobalEnv)
  
  BioAge$FinalBA<-NA
  #Final model
  for(i in 1:length(datBM$Age)){
    BioAge[i,"FinalBA"]<-(sum((datBM[-1][i,]-pars[2,])*pars[1,]/(pars[3,]^2))+datBM[i,"Age"]/Sbsqn)/(sum((pars[1,]^2)/(pars[3,]^2))+1/Sbsqn)
  }
  
  pdf(paste0("output/", set, "Set/FinalBA_", sex, ".pdf"))
  plot(BioAge$Age, BioAge$FinalBA)
  abline(lm(BioAge$Age ~ BioAge$FinalBA))
  dev.off()
  
  return(BioAge)}