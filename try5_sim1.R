# extract the information from fuzzy forest loop object over 51 sim



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
for(i in 1:13){
TT1[mean_fre_[[i]]$mean_imp_per_method$feature_name ,i]=mean_fre_[[i]]$mean_imp_per_method$variable_importance_mean
TT2[names(mean_fre_[[i]]$freq_imp_per_method) ,i]=mean_fre_[[i]]$freq_imp_per_method
}






TT3=TT1[which(rowMeans(TT1)>10),]
TT1=TT3
dev.off()
plot(as.numeric(TT1[,2]), type = "b", frame = FALSE, pch = 10, xlab = "x", ylab = "y", 
     lty = 1, lwd = 2,xaxt='n',ann=FALSE, col=rainbow(13, start = 0, end = 0.85)[1],ylim=c(0,52))

axis(1,at=1:nrow(TT1),labels=rownames(TT1), pch=5)
# 3. Add a second line

for(j in 2:13){
  
  lines(as.numeric(TT1[,j]), pch = 10, type = "b", lty = j, lwd = 2, col=rainbow(13, start = 0, end = 0.85)[j])
  
}

legend("bottomleft", legend = methods, lty = 1, cex = 0.75,lwd =3,col=rainbow(13, start = 0, end = 0.85))
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
TT1=TT2
dev.off()
plot(as.numeric(TT1[,2]), type = "b", frame = FALSE, pch = 10, xlab = "x", ylab = "y", 
     lty = 1, lwd = 2,xaxt='n',ann=FALSE, col=rainbow(14, start = 0, end = 0.85)[1],ylim=c(0,1))

axis(1,at=1:nrow(TT1),labels=rownames(TT1), pch=5)
# 3. Add a second line

for(j in 2:14){
  
  lines(as.numeric(TT1[,j]), pch = 10, type = "b", lty = j, lwd = 2, col=rainbow(14, start = 0, end = 0.85)[j])
  
}

legend("top", legend = c(methods,"RF"), lty = 1, cex = 0.75,lwd =3,col=rainbow(14, start = 0, end = 0.85))
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

#######################################"""
wb = createWorkbook()

sheet = createSheet(wb, "Sheet 1")

addDataFrame(dataframe1, sheet=sheet, startColumn=1, row.names=FALSE)
addDataFrame(dataframe2, sheet=sheet, startColumn=10, row.names=FALSE)

sheet = createSheet(wb, "Sheet 2")

addDataFrame(dataframe3, sheet=sheet, startColumn=1, row.names=FALSE)

saveWorkbook(wb, "My_File.xlsx")


library(xlsx)
write.xlsx(dataframe1, file="filename.xlsx", sheetName="sheet1", row.names=FALSE)
write.xlsx(dataframe2, file="filename.xlsx", sheetName="sheet2", append=TRUE, row.names=FALSE)










