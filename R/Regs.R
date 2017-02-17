require(plyr)
#regressions
#Will handle any data file with a column called "Age" and the rest biomarkers. Regs calculates the parameters for the model. Also assigns a list of variable to be inverted called lsinvert
Regs<-function(datBM,PCAsc=NULL){
  PC<-"z"
  if(!is.null(PCAsc)){
    PC<-"t"
    pca1= prcomp(subset(datBM,select=-Age), scale. = TRUE)
    Agee<-subset(datBM,select=Age)
    datBM<-data.frame(cbind(Agee,pca1$x))
    rot<-as.matrix(pca1$rotation)
    assign("rot",rot,.GlobalEnv)}
  
  runone<-adply(subset(datBM,select=-Age), 2, function(x){
    test<-lm(x[,1] ~ datBM$Age)
    cors<-cor(x[,1], datBM$Age, use = "complete")
    #tes<-ncvTest(test) #will give a p value on heteroscedasticity
    #ad<-shapiro.test(resid(test)) #will give a p value on non-normality of residuals
    out <-data.frame( "Coefficient"=summary(test)$coefficients[2,1], 
                      "Intercept"=summary(test)$coefficients[1,1], 
                      "RMSE"=mean(test$residuals^2)^.5, 
                      "t value"=summary(test)$coefficients[2,3], 
                      "r"=cors,
                      "r.squared"=summary(test)$r.squared, 
                      "p.value"=summary(test)$coefficients[2,4]
                      #,"heteroscedasticity"=tes$p,     #uncomment for diagnostics on regressions
                      #"Shapiro_Wilk_p"=ad$p            
    )
    return(out)
  })
  lsinvert<-runone$X1[runone$Coefficient < 0]
  datBM[,lsinvert]<- datBM[,lsinvert]*(-1)
  assign("lsinvert",c(PC,lsinvert),.GlobalEnv)
  runtwo<-adply(subset(datBM,select=-Age), 2, function(x){
    test<-lm(x[,1] ~ datBM$Age)
    cors<-cor(x[,1], datBM$Age, use = "complete")
    #tes<-ncvTest(test) #will give a p value on heteroscedasticity
    #ad<-shapiro.test(resid(test)) #will give a p value on non-normality of residuals
    out <-data.frame( "Coefficient"=summary(test)$coefficients[2,1], 
                      "Intercept"=summary(test)$coefficients[1,1], 
                      "RMSE"=mean(test$residuals^2)^.5, 
                      "t value"=summary(test)$coefficients[2,3], 
                      "r"=cors,
                      "r.squared"=summary(test)$r.squared, 
                      "p.value"=summary(test)$coefficients[2,4]
                      #,"heteroscedasticity"=tes$p,     #uncomment for diagnostics on regressions
                      #"Shapiro_Wilk_p"=ad$p            
    )
    return(out)
  })
  pars<-t(runtwo)
  colnames(pars)<-pars[1,]
  pars<-pars[-1,]
  pars<-data.frame(pars, row.names = row.names(pars))
  pars <- as.data.frame(sapply(pars, as.numeric), row.names = row.names(pars))
  return(pars)
}