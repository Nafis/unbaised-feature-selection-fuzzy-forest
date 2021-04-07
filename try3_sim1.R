#loop simulation for 



k=4
par(mfrow=c(1,4))


ll_rand_mean=matrix(0,13,13)
ll_ARI_mean=matrix(0,13,13)
ll_NMI_mean=matrix(0,13,13)

kk= seq(100,200)
for(r in kk){
  
  
  
  
  set.seed(r)
  print(r)
  fun=sim_100_(5,0.8,0,0,0.,0.7,0.7,0.7)
  
  data=fun[,-ncol(fun)]
 
  
  ###################################"
  # different clustering méthod with optimal number of clusters
  #######################################
  
  
  ###hierarchical clustering on different methods
  CorData=cor(data)
  cdist<-dist(1-CorData**2)
  Methods=c("ward.D", "ward.D2", "complete", "average", "mcquitty", "centroid","single","median")
  groups=NULL
  ff=list()


  for (i in 1:length(Methods)) {
    hclust_fit <- hclust(cdist, method=Methods[i])
    cah_groups=cutree(hclust_fit, k=k)
    groups=cbind(groups,cah_groups )
  }
  colnames(groups)=Methods
  groups=cbind(groups,"original"=rep(1:4,each=25))
  groups=as.data.frame(groups)
  
  #clusofvar
  
  Clust_Ofvar=hclustvar(X.quanti=data, X.quali=NULL)
  grp_clust_Ofvar=cutreevar(Clust_Ofvar, k=k, matsim=TRUE)
  
  groups$clustofvar=  grp_clust_Ofvar$cluster
  
  
  #varclust:
  
  #if number of cluster is know we should use this function in package varclust:
  var_clust=mlcc.reps(as.matrix(data), numb.clusters = k, numb.runs = 50, numb.cores = 1)
  groups$varclust<-var_clust$segmentation
  
  
  #CLV 
  clv_local=CLV(data, method="local", sX=FALSE)
  grp_local=get_partition(clv_local,K=k,type="vector")
  clv_directional=CLV(data, method="directional",sX=FALSE)
  grp_directional=get_partition(clv_directional,K=k,type="vector")
  groups$clv_local=grp_local
  groups$clv_dir=grp_directional
  
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


ll_rand_mean=ll_rand_mean+ll_rand
ll_ARI_mean=ll_ARI_mean+ll_ARI
ll_NMI_mean=ll_NMI_mean+ll_NMI
}



ran=round(ll_rand_mean/length(kk),3)
ari=round(ll_ARI_mean/length(kk),3)
nmi=round(ll_NMI_mean/length(kk),3)




par(mfrow=c(1,4))
corrplot(cor(data),method = "square" )
corrplot(ran,method="square", title="Rand")
corrplot(ari,method="square", title="ARI")
corrplot(nmi,method="square", title="NMI")






