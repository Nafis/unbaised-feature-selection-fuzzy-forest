rm(list=ls())
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
library(cluster)


#####################################################################################
#  simulation 4th block not correlated to 
######################################################################################"

##################################################################################


sim_100_<-function(m,alpha,alpha4,alpha12,alpha23,alpha34,gamma12,gamma23,gamma34){
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
    Sigma1[A4,A4]<-alpha4
    
    
    
    
  }
  
  diag(Sigma1)<- 1
  
  Mu1=rep(0,100)
  lambda1=rep(0,100)
  lambda1[c(1,2)]=5
  lambda1[3]=2
  lambda1[c(99,100)]=10
  lambda1[98]=5

  epsilon<-rnorm(250,0,0.5)
  X1=mvrnorm(n = 250, Mu1, Sigma1%*%t(Sigma1),empirical = FALSE)
  X1=scale(X1)
  y1=X1%*%lambda1+epsilon
  
  data_sim1=as.data.frame(cbind(X1,y1))
  colnames(data_sim1)=c(paste0("X",1:100),"y")
  
  return(data_sim1)
}






###################################################################################"
# Function that find different clustering partitions 
##################################################################################


ClustVarSim<-function(fun,k){
  
  data=fun[,-ncol(fun)]
  CorData=cor(data)
  cdist <- dist(1 - CorData**2)
  #cdist<-dist(cor(data))
  Methods=c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")
  groups=NULL
  grp_local=c()
  grp_directional=c()

  screen_c <- screen_control(keep_fraction = .25,
                             ntree_factor = 1,
                             min_ntree = 10)
  
  select_c <- select_control(number_selected =100,
                             ntree_factor = 1,
                             min_ntree = 10)
  
  
  WGCNA_c <- WGCNA_control(power = 6, 
                           TOMType = "unsigned", 
                           minModuleSize = 5,
                           numericLabels = TRUE, 
                           pamRespectsDendro = FALSE)

  for (i in 1:length(Methods)) {
    
    hclust_fit <- hclust(cdist, method=Methods[i])
    #k=fviz_nbclust(as.matrix(cdist), FUN=hcut,hc_method=Methods[i], method = c("silhouette"))
    cah_groups=cutree(hclust_fit, k=k)
    groups=cbind(groups,cah_groups )
  }
  #CLV:
  clv_local=CLV(data, method="local", sX=FALSE)
  grp_local=get_partition(clv_local,K=k,type="vector")
  clv_directional=CLV(data, method="directional",sX=FALSE)
  grp_directional=get_partition(clv_directional,K=k,type="vector")
  
  #Clustofvar:
  Clust_Ofvar=hclustvar(X.quanti=data, X.quali=NULL)
  grp_clust_Ofvar=cutreevar(Clust_Ofvar, k=k, matsim=TRUE)
  
  #varclust:
 
  #if number of cluster is know we should use this function in package varclust:
  var_clust=mlcc.reps(as.matrix(data), numb.clusters = k, numb.runs = 50, numb.cores = 1)
  
  #WGCNA
  wf_fit <- wff(y~., data=fun,
                WGCNA_params = WGCNA_c,
                screen_params = screen_c,
                select_params = select_c,
                final_ntree = 1,
                num_processors = 8)
  
  
  
  
  groups<- as.data.frame(cbind(groups,grp_local,grp_directional,grp_clust_Ofvar$cluster,var_clust$segmentation,wf_fit$module_membership[,2] ))
  
  
  
  colnames(groups)<-c(Methods,"clv_local","clv_dir","clustofvar","varclust","WGCNA")

  
  return(list("Modules"=groups))
}
####################################################################"
# It does FF on the differents clustering methods
######################################################################"


fuzzy_clustvar<-function(data, X){
  #X is the matrix of NxP: each column is a partition of N variable
  screen_c <- screen_control(keep_fraction = .25,
                             ntree_factor = 1,
                             min_ntree = 250)
  select_c <- select_control(number_selected = 100,
                             ntree_factor = 1,
                             min_ntree = 250)
  
  
  WGCNA_c <- WGCNA_control(power = 6, 
                           TOMType = "unsigned", 
                           minModuleSize = 5,
                           numericLabels = TRUE, pamRespectsDendro = FALSE)
  
  ff=list()
  for(j in 1:ncol(X)){
    
    ff_fit <- ff(as.data.frame(data[-ncol(data)]), as.vector(data[,ncol(data)]), module_membership = X[,j],screen_params = screen_c,select_params = select_c,final_ntree = 250)
    ff[[j]]<-ff_fit
  }
  
  names(ff)<-colnames(X)
  return(ff)
}



#########################################################
# randomForest on data
###################################################""
varimp=c()

set.seed(123)
sim_fun=sim_100_(5,0.8,0,0,0,0,0.7,0.7,0.7)
corrplot(cor(sim_fun),method="square")

for(i in 1:100){
#sim_fun=randomForest(Species~., data=iris)
  
#sim_fun=sim_100(5,0.8,0,0,0,0.7,0.7,0.7)
rf=randomForest(as.data.frame(sim_fun[,-ncol(sim_fun)]),sim_fun$y)
#rf=randomForest(Species~., data=iris)
varimp=c(varimp,rf$importance[order(rf$importance,decreasing = T),][1:100])
}
#frequency:
nn=names(varimp)
varimp_ord=matrix(nn,100,100)
tt=table(varimp_ord[1:10,])
mm=str_sort(names(tt),numeric=T)
plot(tt[mm],ylab="frequency",xlab="Variables",main=" 100 RF repetition over one simulated data \nfrequency" )
#mean
varimp_imp=matrix(varimp,100,100)


#mean:
dd=as.data.frame(cbind("var"=names(varimp),"importance"=varimp))


rf_varmean=as.data.frame(aggregate(as.numeric(dd$importance),list(dd$var),"mean"))
mm=str_sort(rf_varmean$Group.1,numeric=T )
rownames(rf_varmean)=rf_varmean$Group.1
rf_varmean=rf_varmean[mm,]
plot(rf_varmean[,2],type="b",xlab="variables",ylab="RF-Vimp mean ", main="100 RF repetition over one simulated data \nVariable importance mean")
axis(1,at=1:nrow(rf_varmean),labels=rf_varmean[,1], pch=5)


##########################################################"
# 10 fuzzy forest on sim_fun
#########################################################"


screen_c <- screen_control(keep_fraction = .25,
                           ntree_factor = 1,
                           min_ntree = 10)

select_c <- select_control(number_selected = 100,
                           ntree_factor = 1,
                           min_ntree = 10)


WGCNA_c <- WGCNA_control(power = 6, 
                         TOMType = "unsigned", 
                         minModuleSize = 5,
                         numericLabels = TRUE, 
                         pamRespectsDendro = FALSE)


wff_fit=list()
for(i in 1:100){
#sim_fun=sim_100(5,0.8,0,0,0,0.7,0.7,0.7)
ff_fit <- wff(as.data.frame(sim_fun[-ncol(sim_fun)]), as.vector(sim_fun[,ncol(sim_fun)]),WGCNA_params=WGCNA_c ,screen_params = screen_c,select_params = select_c,final_ntree = 250)
wff_fit[[i]]=ff_fit
}



#frequency:

varnames= unlist(lapply(wff_fit, function(x) x$feature_list$feature_name[1:10] ))
tt=table(varnames)
mm=str_sort(names(tt),numeric=T)
plot(tt[mm],ylab="frequency",xlab="Variables",main="100 fuzzy forest VIM for one simulated data \nfrequency" )

#mean:
varnames=unlist(lapply(wff_fit, function(x) x$feature_list$feature_name ))
varimp_=unlist(lapply(wff_fit, function(x) x$feature_list$variable_importance ))
dd=as.data.frame(cbind("var"=varnames,"importance"=varimp_))
wff_varmean=as.data.frame(aggregate(as.numeric(dd$importance),list(dd$var),"mean"))
mm=str_sort(wff_varmean$Group.1,numeric=T )
rownames(wff_varmean)=wff_varmean$Group.1
wff_varmean=wff_varmean[mm,]
plot(wff_varmean[,2],type="b",xlab="variables",ylab="wff-Vimp mean ", main="100 Fuzzy forest VIMP for one simulated data \nVariable importance mean")
axis(1,at=1:nrow(wff_varmean),labels=wff_varmean[,1], pch=5)


#######################################################################"
# check the linear relation of X and y in function sim_100_
#######################################################################
set.seed(100)
sim_fun=sim_100_(5,0.8,0,0,0,0,0.7,0.7,0.7)
par(mfrow=c(3,3))
plot(sim_fun[,1],sim_fun[,ncol(sim_fun)],xlab="X1",ylab="y")
plot(sim_fun[,2],sim_fun[,ncol(sim_fun)],xlab="X2",ylab="y")
plot(sim_fun[,3],sim_fun[,ncol(sim_fun)],xlab="X3",ylab="y")
plot(sim_fun[,15],sim_fun[,ncol(sim_fun)],xlab="X15",ylab="y")
plot(sim_fun[,52],sim_fun[,ncol(sim_fun)],xlab="X52",ylab="y")
plot(sim_fun[,65],sim_fun[,ncol(sim_fun)],xlab="X65",ylab="y")
plot(sim_fun[,75],sim_fun[,ncol(sim_fun)],xlab="X75",ylab="y")
plot(sim_fun[,99],sim_fun[,ncol(sim_fun)],xlab="X99",ylab="y")
plot(sim_fun[,100],sim_fun[,ncol(sim_fun)],xlab="X100",ylab="y")

plot(sim_fun[,100],sim_fun[,101],xlab="X100",ylab="y")






