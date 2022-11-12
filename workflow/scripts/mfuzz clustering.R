
###Removed CMRN_unex2 and MAL_Per1 and 2 cos they are outliers from the clustering algorithm###

###Remove the ribosomal reads to ensure no differences due to library prep###
rRNA = read.delim('/Users/vickyingham/Seafile/Papers/Funestus/RNAseq/rRNA.txt')
rRNA = as.vector(rRNA$Gene.ID)

data = read.delim('/Users/vickyingham/Seafile/Papers/Funestus/RNAseq/normalised_all.txt',header=T)

library(dplyr)
###Read in unique transcripts
all_data = filter(data, !(GeneID %in% rRNA))


all_data2 = all_data[,2:ncol(all_data)]
row.names(all_data2) = make.names(all_data[,1],unique=T)


all_data2 = data.frame(sapply(all_data2, function(x) as.numeric(as.character(x))))
all_data2 = as.matrix(all_data2)
rownames(all_data2)= make.names(all_data[,1],unique=T)

library(Mfuzz)

require(Biobase)
object<-new("ExpressionSet", exprs=as.matrix(all_data2))

#Since the clustering is performed in Euclidian space, the expression values of genes were standardised to have a mean value of zero and a standard deviation of one.
mestimate<- function(df){
  N <-  dim(df)[[1]]
  D <- dim(df)[[2]]
  m.sj <- 1 + (1418/N + 22.05)*D^(-2) + (12.33/N +0.243)*D^(-0.0406*log(N) - 0.1134)
  return(m.sj)
}

m <- mestimate(all_data2)
m


sumsqr <- function(x, clusters){
  sumsqr <- function(x) sum(scale(x, scale = FALSE)^2)
  wss <- sapply(split(as.data.frame(x), clusters), sumsqr)
  return(wss)
}

#get the wss for repeated clustering
iterate_fcm_WSS <- function(df,m){
  totss <- numeric()
  for (i in 2:20){
    FCMresults <- cmeans(df,centers=i,m=m)
    totss[i] <- sum(sumsqr(df,FCMresults$cluster))
  }
  return(totss)
}

all_data = as.data.frame(all_data2)
wss_2to20 <- iterate_fcm_WSS(all_data,m)
plot(1:20, wss_2to20[1:20], type="b", xlab="Number of Clusters", ylab="wss")


all_data3 = standardise(object)


cluster_object = mfuzz(all_data3,c=10,m=m)

##single allows you to draw one cluster with abline
pdf('/Users/vickyingham/Seafile/Papers/Funestus/RNAseq/Plots_a_page.pdf')
par(mar = par('mar') + c(5,0,0,0))
mfuzz.plot2(all_data3,cl=cluster_object,time.labels=colnames(all_data2),centre=T,colo='fancy',x11 = FALSE,las=2,xlab="",cex.lab=0.6,cex.axis=0.6)
dev.off()


m1 = mestimate(all_data3)

cl = acore(all_data3,cl = cluster_object)

##This writes out each cluster so you can see membership
for(i in 1:length(cl))
{
  write.table(cl[[i]],paste('/Users/vickyingham/Seafile/Papers/Funestus/RNAseq/clusters',i,'.txt',sep=''),row.names=F,sep='\t')
}

library(reshape2)
library(tidyr)
library(dplyr)

#get the centroids into a long dataframe:
fcm_centroids <- cluster_object$centers
fcm_centroids_df <- data.frame(fcm_centroids)
fcm_centroids_df$cluster <- row.names(fcm_centroids_df)
centroids_long <- tidyr::gather(fcm_centroids_df,"sample",'value',-cluster)


###Can only use this palette with n colours <=10
cbp1 <- c("#000000","#999999", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7","#293352","#F0E442")


library(ggplot2)
ggplot(centroids_long, aes(x=sample,y=value, group=cluster, colour=as.factor(cluster))) + 
  geom_point() +
  geom_line()+
  xlab("Data Set") +
  ylab("Expression") +
  labs(title= "Cluster Expression by Species",color = "Cluster") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_discrete(limits=centroids_long$sample) +
  geom_vline(xintercept = c(506),linetype='dashed',colour='black') +
  scale_color_manual(values = cbp1)+
  geom_hline(yintercept = 0,linetype='dashed') 

###Plot only centroind of interest (different in funestus)###
fcm_short = fcm_centroids_df[c(6,7,10),]
centroids_short <- tidyr::gather(fcm_short,"sample",'value',-cluster)

ggplot(centroids_short, aes(x=sample,y=value, group=cluster, colour=as.factor(cluster))) + 
  geom_point() +
  geom_line()+
  xlab("Data Set") +
  ylab("Expression") +
  labs(title= "Cluster Expression by Species",color = "Cluster") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_discrete(limits=centroids_long$sample) +
  geom_vline(xintercept = c(506),linetype='dashed',colour='black') +
  scale_color_manual(values = cbp1)+
  geom_hline(yintercept = 0,linetype='dashed') 

###Don't really want clusters with a correlation of  > 0.9
correlations = cor(t(fcm_centroids))
correlations = cor(all_data2)

library(corrplot)
corrplot(correlations, method="color",tl.cex = 0.6, tl.col = "black")


###Hyper geometric tests for clusters
assigned_names = read.delim('/Users/vickyingham/Seafile/Papers/Funestus/RNAseq/AllGenes_with_added.txt',header=T)
data = read.delim('/Users/vickyingham/Seafile/Papers/Funestus/RNAseq/normalised_all.txt',header=T)

library(dplyr)
###Read in unique transcripts
all_data = filter(data, !(GeneID %in% rRNA))


##Change this to gambiae column
hyper_data = merge(all_data,assigned_names,by.x = 'GeneID',by.y = 'AnGamGeneStableID',all=T)
hyper_data = hyper_data[,-c((ncol(hyper_data)-4):(ncol(hyper_data)-1))]
hyper_data = as.matrix(hyper_data)

hyper_data[is.na(hyper_data[,ncol(hyper_data)])==T,ncol(hyper_data)] = 'Unknown'

assignments = unique(hyper_data[,ncol(hyper_data)])  
assignments = assignments[-which(assignments== 'Unknown')]

library(stringr)
hyper_tests = c()
total_genes = nrow(all_data)
for(j in 1:length(cl))
{
  p_val = c()
  cluster_members = as.data.frame(str_replace(as.character(cl[[j]]$NAME), '\\.', '-'))
  colnames(cluster_members) = 'gene'
  ###Change to gambiae column
  cluster_assignment = merge(cluster_members,hyper_data,by.x='gene',by.y='GeneID',all=F)
  for(i in 1:length(assignments))
  {
    total_on_array = length(which(hyper_data[,ncol(hyper_data)]==assignments[i]))
    total_in_cluster = length(which(cluster_assignment[,ncol(cluster_assignment)]==assignments[i]))
    p_val = c(p_val,p.adjust(1-phyper(total_in_cluster-1,total_on_array,(total_genes-total_on_array),nrow(cluster_members)),method='bonferroni'))
  }
  hyper_tests = rbind(hyper_tests,c(p_val,j))
}
colnames(hyper_tests) = c(assignments,'cluster')

write.table(hyper_tests,'/Users/vickyingham/Seafile/Papers/Funestus/RNAseq/Clusters/hyper_tests.txt',sep='\t',row.names=F)

pca_data <- prcomp(t(all_data2), center = TRUE, scale = TRUE)



##need to add info
library("factoextra")

##Variance Explained
fviz_eig(pca_data)

##cos2 is their quality on the factor map
fviz_pca_ind(pca_data,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,     # Avoid text overlapping
             axes=c(1,2)
)








