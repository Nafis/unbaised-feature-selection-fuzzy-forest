# simulation data 28/12/2020

#simulation

#par(mfrow=c(1,4))
library(factoextra)
library(NbClust)

#set.seed(1)
#fun=fun_1=sim_100(5,0.7,0,0,0,0.6,0.6,0.6) #==> single 
set.seed(123)
fun=sim_100(5,0.8,0,0,0,0,0,0)
#set.seed(123)
#fun=fun_1=sim_100(5,0.8,0,0,0,0.6,0.6,0.6) #==> single

data=fun[,-ncol(fun)]
corrplot(cor(fun),method="square")

###################################"
# different clustering méthod with optimal number of clusters
#######################################


###hierarchical clustering on different methods

cdist<-dist(1-CorData**2)
Methods=c("ward.D", "ward.D2", "complete", "average", "mcquitty", "centroid","single","median")
groups=NULL
ff=list()

par(mfrow=c(2,4))
for (i in 1:length(Methods)) {
    hclust_fit <- hclust(cdist, method=Methods[i])
    plot( hclust_fit)
    cah_groups=cutree(hclust_fit, k=4)
    groups=cbind(groups,cah_groups )
}
  

  
  for (i in 1:length(Methods)) {
  fviz_nbclust(as.matrix(cdist), FUN=hcut,hc_method=Methods[i] ,method = "silhouette")
    print(i)
  }
  colnames(groups)=Methods
  
  groups=cbind(groups,"original"=rep(1:4,each=25))
  
  p=ncol(groups)
  l=NA
  ll_rand=ll_ARI=ll_NMI=matrix(0,p,p)
  
  #calculate indice of comprison
  
  
  for( i in 1:p){
    for(j in 1:p){
      part1=as.integer( groups[,i])
      part2=as.integer( groups[,j])
      l_rand=extCriteria(part1,part2,crit="Rand")# Rand Index (ARI)
      l_ARI=ARI(part1,part2)
      l_NMI=NMI(part1, part2)#varient=max. we can take also the mean of all varient of NMI
      ll_rand[i,j]=as.numeric(l_rand$rand)
      ll_ARI[i,j]=as.numeric(l_ARI)
      ll_NMI[i,j]=as.numeric(l_NMI)
      l=NA
    }
  }
  
  
  colnames(ll_rand)<-rownames(ll_rand)<-colnames(ll_ARI)<-rownames(ll_ARI)<-colnames(ll_NMI)<-rownames(ll_NMI)<-colnames(groups )
  
  # plot Indice of comparison the paritions
  par(mfrow=c(1,4))
  corrplot(cor(data),method = "square" )
  corrplot(ll_rand,method="square", title="Rand")
  corrplot(ll_ARI,method="square", title="ARI")
  corrplot(ll_NMI,method="square", title="NMI")
  

#Clustofvar:
Clust_Ofvar=hclustvar(X.quanti=data, X.quali=NULL)
grp_clust_Ofvar=cutreevar(Clust_Ofvar, k=k, matsim=TRUE)

#varclust:

#if number of cluster is know we should use this function in package varclust:
var_clust=mlcc.reps(as.matrix(data), numb.clusters = k, numb.runs = 50, numb.cores = 1)

#WGCNA
wf_fit <- wff(y~., data=fun,WGCNA_params = WGCNA_c,screen_params = screen_c,select_params = select_c,final_ntree = 1,num_processors = 8)

groups<- as.data.frame(cbind(groups,grp_local,grp_directional,grp_clust_Ofvar$cluster,var_clust$segmentation,wf_fit$module_membership[,2] ))
colnames(groups)<-c(Methods,"clv_local","clv_dir","clustofvar","varclust","WGCNA")
