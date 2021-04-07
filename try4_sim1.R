# fuzzy forest on different clustering méthodes results:
start_time <- Sys.time()
mean_ff=list()

fff=list()
#simulate data
zz=seq(100,150)
for(r in zz){
set.seed(r)
fun=sim_100_(5,0.8,0,0,0.,0.7,0.7,0.7)
#clustering the variables
cluster_res=ClustVarSim(fun, k=4)
#do fuzzy forest of different method of clustering
groups=cluster_res$Modules

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

#par(mfrow=c(1,4))
#corrplot(cor(fun),method = "square" )
#corrplot(ll_rand,method="square", title="rand")
#corrplot(ll_ARI,method="square", title="ARI")
#corrplot(ll_NMI,method="square", title="NMI")
#mtext(paste("seed=",r), side = 3, line = -5, outer = TRUE)

#do fuzzy forest on partitions
ff_fit=fuzzy_clustvar(fun,groups)

fff[[r-99]]=list("ff"=ff_fit, "index_rand"=ll_rand,"index_ARI"=ll_ARI,"index_NMI"=ll_NMI)

}


end_time <- Sys.time()

diff_time<-end_time-start_time
##############################################################################################


#random forests


#results of randm forest
varimp=c()
zz=seq(100,150)
for(r in zz){
  set.seed(r)
  fun=sim_100_(5,0.8,0,0,0.,0.7,0.7,0.7)
  #clustering the variables

  rf=randomForest(as.data.frame(fun[,-ncol(fun)]),fun$y)
  #rf=randomForest(Species~., data=iris)
  varimp=c(varimp,rf$importance[order(rf$importance,decreasing = T),][1:10])
  #do fuzzy forest on partitions
 print(r)
  
}


freq_rf=round(table(names(varimp))/51,3)
