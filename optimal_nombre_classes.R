#find the best number of k for different methods

#Data: 
set.seed(123)
funA=sim_100_(5,0.8,0,0,0,0,0,0,0)
set.seed(123)
funB=sim_100_(5,0.8,0,0,0,0,0.7,0.7,0.7)
fun=funA
data=fun[,-ncol(fun)]


##
# hclust different methods:

CorData=cor(data)
cdist<-dist(1-CorData**2)
Methods=c("ward.D", "ward.D2", "complete", "average", "mcquitty", "centroid","single","median")
# git and plot the optimal number of clusters by "silhouette", "wss" and "gap_stat" method:
for (i in 1:length(Methods)) {
  hclust_fit <- hclust(cdist, method=Methods[i])
  plot(hclust_fit)
  nb_class_sil=fviz_nbclust(as.matrix(cdist), FUN=hcut,hc_method=Methods[i] ,method = "silhouette", k.max=25)
  plot(nb_class_sil,title=method[i])
  nb_class_wss=fviz_nbclust(as.matrix(cdist), FUN=hcut,hc_method=Methods[i] ,method = "wss", k.max=25)
  plot(nb_class_wss)
  nb_class_gap=fviz_nbclust(as.matrix(cdist), FUN=hcut,hc_method=Methods[i] ,method = "gap_stat", k.max=25)
  plot(nb_class_gap)
}
  
  #cah_groups=cutree(hclust_fit, k=k)
  #groups=cbind(groups,cah_groups )



  
  #clusofvar
  
  Clust_Ofvar=hclustvar(X.quanti=data, X.quali=NULL)
  plot(Clust_Ofvar)
  grp_clust_Ofvar=cutreevar(Clust_Ofvar, k=k, matsim=TRUE)
  
  groups$clustofvar=  grp_clust_Ofvar$cluster
  
  
  #varclust(kmeans):
  
  #if number of cluster is know we should use this function in package varclust:
  var_clust=mlcc.reps(as.matrix(data), numb.clusters = k, numb.runs = 50, numb.cores = 1)
  groups$varclust<-var_clust$segmentation
  
  
  #CLV 
  clv_local=CLV(data, method="local", sX=FALSE)
  plot(clv_local)
  grp_local=get_partition(clv_local,K=k,type="vector")
  clv_directional=CLV(data, method="directional",sX=FALSE)
  grp_directional=get_partition(clv_directional,K=k,type="vector")
  plot(clv_directional)
  groups$clv_local=grp_local
  groups$clv_dir=grp_directional
  groups=cbind(groups,"original"=rep(1:4,each=25))
  groups=as.data.frame(groups)
######################################################################""
  
  #chose the k manually
  
  #calculate ll_rand
#####################################################################



  
  par(mfrow=c(1,4))
  
  
  ll_rand_mean=matrix(0,14,14)
  ll_ARI_mean=matrix(0,14,14)
  ll_NMI_mean=matrix(0,14,14)
  
  kk= seq(100,150)
  for(r in kk){
    
    set.seed(r)
    print(r)
    fun=sim_100_(5,0.8,0,0,0,0,0.7,0.7,0.7)
    
    data=fun[,-ncol(fun)]
    
    
    ###hierarchical clustering on different methods
    CorData=cor(data)
    cdist<-dist(1-CorData**2)
    Methods=c("ward.D", "ward.D2", "complete", "average", "mcquitty", "centroid","single","median")
    groups=NULL
    ff=list()
    
    
    for (i in 1:length(Methods)) {
      hclust_fit <- hclust(cdist, method=Methods[i])
      cah_groups=cutree(hclust_fit, k=10)
      groups=cbind(groups,cah_groups )
    }
    colnames(groups)=Methods
    
    groups=as.data.frame(groups)
    
    #clusofvar
    
    Clust_Ofvar=hclustvar(X.quanti=data, X.quali=NULL)
    grp_clust_Ofvar=cutreevar(Clust_Ofvar, k=3, matsim=TRUE)
    
    groups$clustofvar=  grp_clust_Ofvar$cluster
    
    
    #varclust:
    
    #if number of cluster is know we should use this function in package varclust:
    var_clust=mlcc.reps(as.matrix(data), numb.clusters = 4, numb.runs = 50, numb.cores = 1)
    groups$varclust<-var_clust$segmentation
    

    #CLV 
    clv_local=CLV(data, method="local", sX=FALSE)
    grp_local=get_partition(clv_local,K=4,type="vector")
    clv_directional=CLV(data, method="directional",sX=FALSE)
    grp_directional=get_partition(clv_directional,K=3,type="vector")
    groups$clv_local=grp_local
    groups$clv_dir=grp_directional
    
    #WGCNA
    
    wf_fit <- wff(y~., data=fun,WGCNA_params = WGCNA_c,screen_params = screen_c,select_params = select_c,final_ntree = 1,num_processors = 8)
    groups$WGCNA=wf_fit$module_membership[,2]
    
  
    
    
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
  
  
  
#############################################"
  
  #K manually
  ClustVarSim_manu<-function(fun){
    
    data=fun[,-ncol(fun)]
    CorData=cor(data)
    cdist <- dist(1 - CorData**2)
    #cdist<-dist(cor(data))
    Methods=c("ward.D", "complete")
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
      cah_groups=cutree(hclust_fit, k=10)
      groups=cbind(groups,cah_groups )
    }
    #CLV:
    clv_local=CLV(data, method="local", sX=FALSE)
    grp_local=get_partition(clv_local,K=4,type="vector")
    clv_directional=CLV(data, method="directional",sX=FALSE)
    grp_directional=get_partition(clv_directional,K=3,type="vector")
    
    #Clustofvar:
    Clust_Ofvar=hclustvar(X.quanti=data, X.quali=NULL)
    grp_clust_Ofvar=cutreevar(Clust_Ofvar, k=3, matsim=TRUE)
    
    #varclust:
    
    #if number of cluster is know we should use this function in package varclust:
    var_clust=mlcc.reps(as.matrix(data), numb.clusters = 4, numb.runs = 50, numb.cores = 1)
    
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
  

#############################################"
  
  
  
  start_time <- Sys.time()
  mean_ff=list()
  
  fff=list()
  #simulate data
  zz=seq(100,150)
  for(r in zz){
    set.seed(r)
    fun=sim_100_(5,0.8,0,0,0,0.,0.7,0.7,0.7)
    #clustering the variables
    cluster_res=ClustVarSim_manu(fun)
    #do fuzzy forest of different method of clustering
    groups=cluster_res$Modules
    
    p=ncol(groups)
    #l=NA
    #ll_rand=ll_ARI=ll_NMI=matrix(0,p,p)
    
    #calculate indice of comprison
    
    #for( i in 1:p){
     # for(j in 1:p){
    #    part1=as.integer( groups[,i])
    #    part2=as.integer( groups[,j])
    #    l_rand=extCriteria(part1,part2,crit="Rand")# Rand Index (ARI)
    #    l_ARI=ARI(part1,part2)
     #   l_NMI=NMI(part1, part2)#varient=max. we can take also the mean of all varient of NMI
      #  ll_rand[i,j]=as.numeric(l_rand$rand)
       # ll_ARI[i,j]=as.numeric(l_ARI)
        #ll_NMI[i,j]=as.numeric(l_NMI)
        #l=NA
      #}
    #}
    #colnames(ll_rand)<-rownames(ll_rand)<-colnames(ll_ARI)<-rownames(ll_ARI)<-colnames(ll_NMI)<-rownames(ll_NMI)<-colnames(groups )
    
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
    fun=sim_100_(5,0.8,0,0,0,0.,0.7,0.7,0.7)
    #clustering the variables
    
    rf=randomForest(as.data.frame(fun[,-ncol(fun)]),fun$y)
    #rf=randomForest(Species~., data=iris)
    varimp=c(varimp,rf$importance[order(rf$importance,decreasing = T),][1:10])
    #do fuzzy forest on partitions
    print(r)
    
  }
  
  
  freq_rf=round(table(names(varimp))/51,3)
  
##################################################
#plot the results
#################################################
  
  
  #summary the results in many matrix
  #fff_manuel=fff

  fff1=fff
  
  take_mean_freq<-function(x){
    
    m=NULL
    f=NULL
    for(i in 1:length(x)){
      xx=x[[i]]
      m=rbind(m,xx[1:10,])
      f=rbind(f,xx[1:10,])
    }
    mm= m%>%group_by(feature_name )%>%summarise(variable_importance_mean=mean(variable_importance) )
    ff=table(f[,1])/length(x)
    return(list("mean_imp_per_method"=mm,"freq_imp_per_method"=round(ff,3)))
  }
  
  #####################################################################
  # average and frequecny of variable important measure per method over 51 simulations
  ###############################################################"
  
  nn=list()
  mean=list()
  freq=list()
  imp=list()
  for(j in 1:13){
    nn[[j]]=unlist(lapply(fff,function(x) nrow(x$ff[[j]]$feature_list)))
    
    imp[[j]]=lapply(fff,function(x) x$ff[[j]]$feature_list)
  }
  
  for(j in 1:13){
    mean_fre_=lapply(imp , function(x) take_mean_freq(x))
  }
  
  names(mean_fre_)=names(fff[[1]]$ff)
  
  ###################################################""""
  #plot the average of Vimp per method:
  
  r=paste0("X",1:100)
  methods=names(imp)=names(mean_fre_)
  
  kk=matrix(0,100,length(methods))
  TT1=as.data.frame(kk)
  
  colnames(TT1)=c(methods)
  rownames(TT1)=r
  TT2=TT1
  for(i in 1:7){
    TT1[mean_fre_[[i]]$mean_imp_per_method$feature_name ,i]=mean_fre_[[i]]$mean_imp_per_method$variable_importance_mean
    TT2[names(mean_fre_[[i]]$freq_imp_per_method) ,i]=mean_fre_[[i]]$freq_imp_per_method
  }
  
  
  
  
  
  
  TT3=TT1[which(rowMeans(TT1)>10),]
  TT1=TT3
  dev.off()
  plot(as.numeric(TT1[,2]), type = "b", frame = FALSE, pch = 10, xlab = "x", ylab = "y", 
       lty = 1, lwd = 2,xaxt='n',ann=FALSE, col=rainbow(8, start = 0, end = 0.90)[1],ylim=c(0,130))
  
  axis(1,at=1:nrow(TT1),labels=rownames(TT1), pch=5)
  # 3. Add a second line
  
  for(j in 2:7){
    
    lines(as.numeric(TT1[,j]), pch = 10, type = "b", lty = j, lwd = 2, col=rainbow(8, start = 0, end = 0.90)[j])
    
  }
  
  legend("topleft", legend = methods, lty = 1, cex = 0.75,lwd =3,col=rainbow(8, start = 0, end = 0.90))
  title(main="VIM mean for 51 simutaions \nper different classification methods")
  
  
  
  #######################################################
  # plot the frequency + rf(from try4)
  varimp=c()
  zz=seq(100,150)
  for(r in zz){
    set.seed(r)
    fun=sim_100_(5,0.8,0,0,0,0.,0.7,0.7,0.7)
    #clustering the variables
    
    rf=randomForest(as.data.frame(fun[,-ncol(fun)]),fun$y)
    #rf=randomForest(Species~., data=iris)
    varimp=c(varimp,rf$importance[order(rf$importance,decreasing = T),][1:10])
    #do fuzzy forest on partitions
    print(r)
    
  }
  
  
  freq_rf=round(table(names(varimp))/51,3)
  TT2$RF=0
  TT2$RF[which(rownames(TT2)%in%names(freq_rf))]=freq_rf
  
  #
  
  TT4=TT2[which(rowSums(TT2)>1.2),]
  TT1=TT2
  dev.off()
  plot(as.numeric(TT1[,2]), type = "b", frame = FALSE, pch = 10, xlab = "x", ylab = "y", 
       lty = 1, lwd = 2,xaxt='n',ann=FALSE, col=rainbow(8, start = 0, end = 0.90)[1],ylim=c(0,1))
  
  axis(1,at=1:nrow(TT1),labels=rownames(TT1), pch=5)
  # 3. Add a second line
  
  for(j in 2:8){
    
    lines(as.numeric(TT1[,j]), pch = 10, type = "b", lty = j, lwd = 2, col=rainbow(8, start = 0, end = 0.90)[j])
    
  }
  
  legend("top", legend = c(methods,"RF"), lty = 1, cex = 0.75,lwd =3,col=rainbow(8, start = 0, end = 0.90))
  title(main=" feature selection performance for 51 simutaions \nper n=250, p=100, M4 indeendent")
  
  
  
  
  
  
  
  
  
  #save(mean_fre_,file="mean_frequency_per_method_ff_50_sim_")
  #################################################################"
  #construct a table
  
  rr=paste0("X",1:100)
  methods=names(imp)=names(mean_fre_)
  
  TT1=as.data.frame(rr)
  kk=matrix(0,100,length(methods))
  TT1=cbind(TT1,kk)
  colnames(TT1)=c("name",methods)
  rownames(TT1)=TT1$name
  
  ALL=c()
  for(i in 1:51){
    TT=TT1 
    
    for(j in 2:14){
      TT[which(rownames(TT)%in%imp[[j-1]][[i]]$feature_name),j]=round(imp[[j-1]][[i]]$variable_importance,2)
    }
    ALL[[i]]=TT
  }
  #save(ALL, file="ff_vimp_per_method_50")
  ###########################################""
  
  
  
  
  
  
  
  #summary the results in many matrix
  
  #fff=ff_50_sim
  fff1=fff
  
  take_mean_freq<-function(x){
    
    m=NULL
    f=NULL
    for(i in 1:length(x)){
      xx=x[[i]]
      m=rbind(m,xx[1:10,])
      f=rbind(f,xx[1:10,])
    }
    mm= m%>%group_by(feature_name )%>%summarise(variable_importance_mean=mean(variable_importance) )
    ff=table(f[,1])/length(x)
    return(list("mean_imp_per_method"=mm,"freq_imp_per_method"=round(ff,3)))
  }
  
  #####################################################################
  # average and frequecny of variable important measure per method over 51 simulations
  ###############################################################"
  
  nn=list()
  mean=list()
  freq=list()
  imp=list()
  for(j in 1:7){
    nn[[j]]=unlist(lapply(fff,function(x) nrow(x$ff[[j]]$feature_list)))
    
    imp[[j]]=lapply(fff,function(x) x$ff[[j]]$feature_list)
  }
  
  for(j in 1:7){
    mean_fre_=lapply(imp , function(x) take_mean_freq(x))
  }
  
  names(mean_fre_)=names(fff[[1]]$ff)
  
  ###################################################""""
  #plot the average of Vimp per method:
  
  r=paste0("X",1:100)
  methods=names(imp)=names(mean_fre_)
  
  kk=matrix(0,100,length(methods))
  TT1=as.data.frame(kk)
  
  colnames(TT1)=c(methods)
  rownames(TT1)=r
  TT2=TT1
  for(i in 1:7){
    TT1[mean_fre_[[i]]$mean_imp_per_method$feature_name ,i]=mean_fre_[[i]]$mean_imp_per_method$variable_importance_mean
    TT2[names(mean_fre_[[i]]$freq_imp_per_method) ,i]=mean_fre_[[i]]$freq_imp_per_method
  }
  
  
  
  
  
  
  TT3=TT1[which(rowMeans(TT1)>10),]
  TT1=TT3
  dev.off()
  plot(as.numeric(TT1[,2]), type = "b", frame = FALSE, pch = 10, xlab = "x", ylab = "y", 
       lty = 1, lwd = 2,xaxt='n',ann=FALSE, col=rainbow(7, start = 0, end = 0.85)[1],ylim=c(0,52))
  
  axis(1,at=1:nrow(TT1),labels=rownames(TT1), pch=5)
  # 3. Add a second line
  
  for(j in 2:13){
    
    lines(as.numeric(TT1[,j]), pch = 10, type = "b", lty = j, lwd = 2, col=rainbow(7, start = 0, end = 0.85)[j])
    
  }
  
  legend("bottomleft", legend = methods, lty = 1, cex = 0.75,lwd =3,col=rainbow(7, start = 0, end = 0.85))
  title(main="VIM mean for 51 simutaions \nper different classification methods")
  
  
  
  #######################################################
  # plot the frequency + rf(from try4)
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
  TT2$RF=0
  TT2$RF[which(rownames(TT2)%in%names(freq_rf))]=freq_rf
  
  #
  
  TT4=TT2[which(rowSums(TT2)>1.2),]
  TT1=TT4
  dev.off()
  par(mfrow=c(2,4))
  for(jj in 1:8){
  plot(as.numeric(TT1[,jj]), type = "b", frame = FALSE, pch = 10, xlab = "x", ylab = "y", 
       lty = 1, lwd = 2,xaxt='n',ann=FALSE, col='red',ylim=c(0,1))
  
  axis(1,at=1:nrow(TT1),labels=rownames(TT1), pch=5)
  # 3. Add a second line
  
  for(j in 2:8){
    
    lines(as.numeric(TT1[,j]), pch = 10, type = "b", lty = j, lwd = 2, col='grey')
    lines(as.numeric(TT1[,jj]), pch = 10, type = "b", lty = j, lwd = 2, col='red')
  }
  
  
  title(main=colnames(TT1)[jj])
  }
  
  
  
  
  
  
  
  
  #save(mean_fre_,file="mean_frequency_per_method_ff_50_sim_")
  #################################################################"
  #construct a table
  
  rr=paste0("X",1:100)
  methods=names(imp)=names(mean_fre_)
  
  TT1=as.data.frame(rr)
  kk=matrix(0,100,length(methods))
  TT1=cbind(TT1,kk)
  colnames(TT1)=c("name",methods)
  rownames(TT1)=TT1$name
  
  ALL=c()
  for(i in 1:51){
    TT=TT1 
    
    for(j in 2:14){
      TT[which(rownames(TT)%in%imp[[j-1]][[i]]$feature_name),j]=round(imp[[j-1]][[i]]$variable_importance,2)
    }
    ALL[[i]]=TT
  }
  #save(ALL, file="ff_vimp_per_method_50")
  
  
  
  
  
  
  
