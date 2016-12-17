
options(stringsAsFactors = FALSE) #do this before importing data

#imported csv BA.csv
#GROOMING, these steps are as needed

BA$Gender[BA$Gender == "M"]<- 0
BA$Gender[BA$Gender == "F"]<- 1
BA$Gender[BA$Gender == BA$Gender[1]]<- 0
BA$Gender[BA$Gender == BA$Gender[9]]<- 1

row.names(BA)<-BA$Subject.ID
BA<-BA[-1]
BA$Gender<-as.numeric(as.character(BA$Gender))
# at this point, BA is a data frame with rows named by subject ID,
#so BA only contains numeric data

########
#optional: split across gender, note men=0, fem=1
BAmen<-BA[BA$Gender==0,]
BAfem<-BA[BA$Gender==1,]

BAmen<-BAmen[-2] #2 is column number of Gender, change this as needed
BAfem<-BAfem[-2]

###############################
###############################
###############################
###############################
###############################
#the findOutliers function finds outliers in each column
# set "cutoff" to desired number of sd
findOutlier <- function(x, cutoff=3) {
  sds <- apply(x, 2, sd, na.rm=TRUE)
  meeen <- apply(x, 2, mean, na.rm=TRUE)
  result <- mapply(function(d, s, m) {
    which(d > m + cutoff*s | d < m - cutoff*s)
  },
  x, sds, meeen
  )
  result
}

#actually generating a list of ouliers
outliers <- findOutlier(BAmen) 

##########################
##########################
#removeOutlier is a function that deletes values from a
#list "outliers"
removeOutlier <- function(x, outliers) {
  result <- mapply(function(d,o) {
    res <- d
    res[o] <- NA
    return(res)
  }, x, outliers)
  return(as.data.frame(result, row.names = row.names(x)))
}

#actually removing outliers
BAmen <- removeOutlier(BAmen, outliers)

###############################
###############################
###############################
###############################
###############################

#Check regression assumptions
#generates a pdf where each page prints the graphical depictions
#of regression assumptions for each variable (in order, not named)

somePDFPath = "~/Desktop/workspaceNovember/Massumptions.pdf"#change as needed
pdf(file=somePDFPath)  
for (i in seq(2,length(BAmen)))   
{   ex<-lm(BAmen[,i] ~ BAmen$Age)
    par(mfrow = c(2,2)) 
    plot(ex)   
} 
dev.off() 

###############################
###############################
###############################
###############################
###############################

#regressions. This makes the file needed for the model, and does a few tests (p values reported)
#for the assumptions on regression
#the BA[-1] removes the age column, avoiding a "perfect fit" error

AgeRegressionsMen<-adply(BAmen[-1], 2, function(x){
  test<-lm(x[,1] ~ BAmen$Age)
  cors<-cor(x[,1], BAmen$Age, use = "complete")
  tes<-ncvTest(test) #will give a p value on heteroscedasticity
  ad<-shapiro.test(resid(test)) #will give a p value on non-normality of residuals
  out <-data.frame( "Coefficient"=summary(test)$coefficients[2,1], "Intercept"=summary(test)$coefficients[1,1], "RMSE"=mean(test$residuals^2)^.5, "t value"=summary(test)$coefficients[2,3], "r"=cors, "r.squared"=summary(test)$r.squared, "p.value"=summary(test)$coefficients[2,4], "heteroscedasticity"=tes$p, "Shapiro_Wilk_p"=ad$p)
  return(out)
})

#find the biomarkers that give a negative regression coefficient, then replace that marker with its negative
x<-AgeRegressionsMen$X1[AgeRegressionsMen$Coefficient < 0]
BAmen[,x]<- BAmen[,x]*(-1)
#note: re-run the regressions file after this step


###############################
###############################
###############################

#write a csv file with regressions output, puts in working directory
write.csv(AgeRegressionsMen,file="AgeRegressionsMen.csv")

###############################
###############################
#transpose and reformat for use
vars<-t(AgeRegressionsMen)
colnames(vars)<-vars[1,]
vars<-vars[-1,]
vars<-data.frame(vars, row.names = row.names(vars))
vars <- as.data.frame(sapply(vars, as.numeric), row.names = row.names(vars))

###############################
#The Model#
###############################
BioAgeMen<-data.frame(cbind(BAmen$Age,adply(BAmen[-1], 1, function(x){
  su<- sum(((x[1,]-vars[2,])*vars[1,])/(vars[3,]^2))/sum((vars[1,]^2)/(vars[3,]^2))
  out<- data.frame(su)
  return(out)
})$su), row.names=row.names(BAmen))

colnames(BioAgeMen)<-c('Age','bioage')
#BioAgeMen now contains rows by subject ID, with age and bioage

#characteristic r of the model (0 to 1, bigger is better)
rchar =sum(vars[5,]*(vars[1,])/(vars[3,]))/sum((vars[1,])/(vars[3,]))
