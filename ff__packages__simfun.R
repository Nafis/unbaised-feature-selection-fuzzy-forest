rm(list=ls())

#list of packages

library(pdp)
library(randomForest)
library(party)
library(gtools)
library(fuzzyforest)
library(MASS)
library(WGCNA)
library(factoextra)
library(NbClust)
library(corrplot)
library(ClustVarLV)
library(ClustOfVar)
library(varclust)
library(dplyr)
library(ggplot2)
library(clusterCrit)
library(aricode)
library(stringr)




#####################################################################################
#  simulation 250 data 4th block not correlated
######################################################################################"


sim_100_<-function(m,alpha,alpha12,alpha23,alpha34,gamma12,gamma23,gamma34){
  #m: the number of corelated features in two different modules  
  Sigma1 <- matrix(0,100,100)
  #correlation between the features inside module 1 is 0.8 
  A1=1:25#sample(1:25,sample(15:25,1))
  
  #correlation between the features inside module 2 is 0.8
  A2=26:50#sample(26:50,sample(15:25,1))
  
  #correlation between the features inside module 3 is 0.8
  A3=51:75#sample(51:75,sample(15:25,1))
  
  #the features in module 4 are independent from each other
  A4=76:100#sample(76:100,sample(15:25,1))
  
  
  
  if(m!=0){
    ind1=1:(2*m)
    ind2=1:m
    ind3=25:(25-m)
    #correlation within groups
    Sigma1[A1,A2]<-alpha12
    Sigma1[A2,A1]<-alpha12
    Sigma1[A2,A3]<-alpha23
    Sigma1[A3,A2]<-alpha23
    Sigma1[A3,A4]<-alpha34
    Sigma1[A4,A3]<-alpha34
    #correlation with groups
    Sigma1[A1,A2[ind1]]<-gamma12/4
    Sigma1[A2[ind1],A1]<-gamma12/4
    Sigma1[A2,A3[ind1]]<-gamma23/4
    Sigma1[A3[ind1],A2]<-gamma23/4
    
    Sigma1[A1,A2[ind2]]<-gamma12
    Sigma1[A2[ind2],A1]<-gamma12
    
    Sigma1[A2,A3[ind2]]<-gamma23
    Sigma1[A3[ind2],A2]<-gamma23
    
    Sigma1[A3[ind3],A1]<-gamma34
    Sigma1[A1,A3[ind3]]<-gamma34
    
    Sigma1[A1[ind2],A2]<-gamma12
    Sigma1[A2, A1[ind2]]<-gamma12
    #Sigma1[A3[ind2],A1]<-gamma12
    #Sigma1[A1, A3[ind2]]<-gamma12
    
    Sigma1[A1,A1]<-alpha
    Sigma1[A2,A2]<-alpha
    Sigma1[A3,A3]<-alpha
    Sigma1[A4,A4]<-0
    
    
    
    
  }
  
  diag(Sigma1)<- 1
  
  Mu1=rep(0,100)
  lambda1=rep(0,100)
  lambda1[c(1,2)]=5
  lambda1[c(3)]=2
  lambda1[c(98,99,100)]=10
  
  epsilon<-rnorm(250,0,0.5)
  X1=mvrnorm(n = 250, Mu1, Sigma1%*%t(Sigma1),empirical = FALSE)
  X1=scale(X1)
  y1=X1%*%lambda1+epsilon
  
  data_sim1=as.data.frame(cbind(X1,y1))
  colnames(data_sim1)=c(paste0("X",1:100),"y")
  
  return(data_sim1)
}



###################################################"""