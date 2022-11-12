counts_file = read.delim('/Users/vickyingham/Seafile/Papers/Funestus/RNAseq/gambiae_complex_counts.txt',header=T)
sampleData = read.delim('/Users/vickyingham/Seafile/Papers/Funestus/RNAseq/ALLcoldata_nofun.txt',header=T)

library(DESeq2)
library(tidyverse)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)
library(ggplot2)
library(limma)

counts = counts_file[,9:ncol(counts_file)]
counts2=counts
row.names(counts) = make.names(counts_file[,1],unique=T)


counts2$Geneid = rownames(counts)

counts = counts2 %>% column_to_rownames('Geneid')
total_samples = colSums(counts)
tot_counts = cbind(colnames(counts),total_samples)
tot_counts = as.data.frame(tot_counts)
tot_counts$total_samples = as.numeric(as.character(tot_counts$total_samples))


g = ggplot(tot_counts, aes(x=V1, y=total_samples)) + 
  geom_bar(stat='identity') + 
  theme_light() +
  ggtitle("Total reads counted (mapped to Ag (PEST))") +
  theme(axis.text.x = element_text(angle=45,size=5)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab('Number of reads') +
  xlab('Sample')

ggsave('/Users/vickyingham/Seafile/Papers/Funestus/RNAseq/counts.pdf',g)

########### make DESeq dataset for all counts ##############

###When doing all R vs all S can correct for batch effects####

dds = DESeqDataSetFromMatrix(countData = counts, 
                             colData = sampleData, 
                             design = ~resistance)

dds = estimateSizeFactors(dds)
dds = estimateDispersions(dds)
vsd = varianceStabilizingTransformation(dds)

mat <- assay(vsd)
mat <- limma::removeBatchEffect(mat, vsd$Batch)
assay(vsd) <- mat
vstcounts <- assay(vsd)

vstcounts = vstcounts[order(apply(vstcounts,1,sum),decreasing =TRUE),]

correlations = cor(vstcounts)
pheatmap(correlations,fontsize_row=5)

conds=as.factor(as.character(sampleData$condition))

library(randomcoloR)
cond_colours = distinctColorPalette(length(levels(conds)))[conds]
names(cond_colours)=sampleData$condition

pca2=prcomp(t(vstcounts),center=TRUE)
plot(pca2$x[,1],pca2$x[,2], col=cond_colours,  pch=19, cex=2, main="Sample to Sample PCA (VST)",xlab='PC1',ylab='PC2')
text(pca2$x[,1],pca2$x[,2], as.vector(colnames(counts)), pos=3, cex=0.4)

pca2=prcomp(t(vstcounts),center=TRUE)
plot(pca2$x[,2],pca2$x[,3], col=cond_colours,  pch=19, cex=2, main="Sample to Sample PCA (VST)",xlab='PC2',ylab='PC3')
text(pca2$x[,2],pca2$x[,3], as.vector(colnames(counts)), pos=3, cex=0.4)


library("factoextra")
pdf("/Users/vickyingham/Seafile/Papers/Funestus/RNAseq/DEseq/pca_variance_all.pdf")
fviz_eig(pca2)
dev.off()

sampleData$batch = factor(sampleData$batch)
dds = DESeqDataSetFromMatrix(countData = counts, 
                             colData = sampleData, 
                             design = ~resistance + batch)

dds$resistance<- relevel(dds$resistance, ref = "susceptible")
dds$species <- relevel(dds$species, ref = "arabiensis")
dds <- DESeq(dds)
resultsNames(dds)

plotDispEsts( dds, ylim = c(1e-6, 1e1) )

library(apeglm)
library(ashr)

type_A_result = results(dds, contrast=c("resistance","resistant","susceptible"))

type_A_dds <- lfcShrink(dds, coef=c("resistance_resistant_vs_susceptible"), type="apeglm")

sum(type_A_dds$pvalue < 0.05, na.rm=TRUE )

plotMA(type_A_dds)

write.table(type_A_dds,paste('/Users/vickyingham/Seafile/Papers/Funestus/RNAseq/','gambiae_vs_arabiensis_resistance_corrected.txt',sep=''),row.names=T,sep='\t')

library(viridis)
library(hrbrthemes)

library(pals)
col = glasbey(n=29)

goi = 'AGAP008212'
d <- plotCounts(dds, gene=goi,  intgroup="resistance",
                returnData=TRUE)
d$condition = sampleData$condition
d$species = sampleData$species


ggplot(d, aes(x=resistance, y=count,fill=resistance)) + 
  geom_boxplot(outlier.color=NA) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6,option="A") +
  geom_point(position=position_jitter(w=0.1,h=0),aes(color=condition,shape=species),size=3) + 
  scale_y_log10() + ggtitle(goi) + xlab('Resistance Status') + ylab('Normalised Counts') +
  theme_minimal() +
  scale_fill_manual(values=alpha(c("#56B4E9","#D55E00"),0.4)) +
  scale_color_manual(values=alpha(col,0.9))



###Plotting for multiple genes!####

gsoi = read.delim('/Users/vickyingham/Seafile/Papers/Funestus/RNAseq/Top_20.txt',header=F)
gsoi = as.matrix(gsoi)
plot_data=c()
for(i in 1:length(gsoi))
{
  skip=FALSE
  
  tryCatch(plotCounts(dds, gene=gsoi[i],  intgroup="resistance",returnData=TRUE), error = function(e) { skip <<- TRUE})

  if(skip == FALSE)
  {
    to_plot = plotCounts(dds, gene=gsoi[i],  intgroup="resistance",
                         returnData=TRUE)
    to_plot$gene = rep(gsoi[i],nrow(to_plot))
    to_plot$condition = sampleData$condition
    to_plot$species = sampleData$species
    plot_data = rbind(plot_data,to_plot)
  }
  else
  {
    next
  }
}

#plot_data$condition = rep(sampleData$condition,length(unique(plot_data$gene)))
#plot_data$species = rep(sampleData$species,length(unique(plot_data$gene)))

ggplot(plot_data, aes(x=resistance, y=count,fill=resistance)) + 
  geom_boxplot(outlier.color=NA) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6,option="A") +
  geom_point(position=position_jitter(w=0.1,h=0),aes(color=condition,shape=species),size=3) + 
  scale_y_log10() + ggtitle('Genes of interest') + xlab('Resistance Status') + ylab('Normalised Counts') +
  theme_minimal() +
  scale_fill_manual(values=alpha(c("#56B4E9","#D55E00"),0.4)) +
  scale_color_manual(values=alpha(col,0.9)) +
  theme(axis.text.x = element_text(angle=45)) +
  facet_wrap(~gene,ncol=10,scales='free')




#####Species differences#####
col = glasbey(n=29)

goi = 'AGAP012572.1'
d <- plotCounts(dds, gene=goi,  intgroup="species",
                returnData=TRUE)
d$condition = sampleData$condition
d$species = sampleData$species


ggplot(d, aes(x=species, y=count,fill=species)) + 
  geom_boxplot(outlier.color=NA) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6,option="A") +
  geom_point(position=position_jitter(w=0.1,h=0),aes(color=condition,shape=species),size=3) + 
  scale_y_log10() + ggtitle(goi) + xlab('Resistance Status') + ylab('Normalised Counts') +
  theme_minimal() +
  scale_fill_manual(values=alpha(c("#56B4E9","#D55E00",'#000000'),0.4)) +
  scale_color_manual(values=alpha(col,0.9))



####Plot for res vs sus gambiae complex

col = glasbey(n=29)

goi = 'AGAP006222'
d <- plotCounts(dds, gene=goi,  intgroup="condition3",
                returnData=TRUE)
d$resistance = sampleData$resistance
d$species = sampleData$species
d$resistance <- factor(d$resistance, levels = c('susceptible','resistant'))


ggplot(d, aes(x=resistance, y=count,fill=resistance)) + 
  geom_boxplot(outlier.color=NA) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6,option="A") +
  geom_point(position=position_jitter(w=0.1,h=0),aes(color=condition3,shape=species),size=3) + 
  scale_y_log10() + ggtitle(goi) + xlab('Resistance Status') + ylab('Normalised Counts') +
  theme_minimal() +
  scale_fill_manual(values=alpha(c("#56B4E9","#D55E00"),0.4)) +
  scale_color_manual(values=alpha(col,0.9))

