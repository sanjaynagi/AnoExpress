## Differential expression
library(DESeq2)
library(pheatmap)
library(data.table)
library(ggrepel)
library(openxlsx)
library(glue)
library(RColorBrewer)
library(tidyverse)
library(plotly)
######

name = 'meta'


metadata = fread("config/sample_metadata.tsv") %>% rename("sampleID"="colData")
metadata[(metadata$batch == 1) &(metadata$species == 'arabiensis'), 'batch'] = 0
metadata$batch = factor(metadata$batch)
counts_meta = fread("results/final_raw_counts.tsv") %>% as.data.frame()
counts = counts_meta[,10:length(counts_meta)] # remove info columns
rownames(counts) = make.unique(counts_meta$GeneID)

gene_names = fread("../rna-seq-busia/resources/exampleGene2TranscriptMap.tsv", sep="\t") %>% distinct()
names_df = data.frame("GeneID2" = make.unique(counts_meta$GeneID))
names_df$GeneID = gsub("\\..*","",names_df$GeneID2)
names_df = names_df %>% left_join(gene_names)


### PCA 
# make DESeq dataset
dds = DESeqDataSetFromMatrix(countData = counts, 
                             colData = metadata, 
                             design = ~ batch)
###### estimate paramters and normalise 
dds = estimateSizeFactors(dds)
dds = estimateDispersions(dds)
vsd = varianceStabilizingTransformation(dds)
vstcounts = assay(vsd)
vstcounts = vstcounts[order(apply(vstcounts,1,sum),decreasing=TRUE),]

#### write pca of samples to pdf
pca2 = prcomp(t(vstcounts),center=TRUE)
pc = data.frame(pca2$x) %>% rownames_to_column("sampleID")
pc = left_join(pc, metadata)
### plot interactive PCA data 
plot_ly(data=pc, x=~PC1, y=~PC2, color=~species,   
        text = ~paste("dataset: ", condition , '<br>resistance:', resistance))

vstcounts_df = data.frame(vstcounts)
vstcounts_df %>% round_df(2) %>% fwrite("results/vst_counts.tsv.gz", sep="\t")

#### heatmap ####
correlations = cor(vstcounts)
pdf("results/plots/heatmap_correlations.pdf")
pheatmap(correlations)
garbage = dev.off()


##### meta analysis #####
table(metadata$batch, metadata$species)
table(metadata$batch, metadata$resistance)





volcano = function(data, title){
  
  data = data %>% mutate(col = case_when(padj < 0.05 & abs(log2FoldChange) > 1 ~ "#ff6666",
                                         padj > 0.05 & abs(log2FoldChange) > 1 ~ "#39e600",
                                         padj < 0.05 & abs(log2FoldChange) < 1 ~ "#3385ff",
                                         padj > 0.05 & abs(log2FoldChange) < 1 ~ "#b3b3b3"))
  
  data$labels = data %>% dplyr::mutate("Gene_name" = case_when(GeneName == "" ~ GeneID,
                                                               is.na(GeneName) ~ GeneID,
                                                               TRUE ~ GeneName)) %>% select(Gene_name) %>% deframe()
  
  plot = ggplot(data=data, aes(x=log2FoldChange, y=-log10(padj), color=col, alpha=0.4), size=2) + 
    geom_point() + 
    scale_color_identity() +
    xlim(-10, 10) +
    ylab("-log10(P)") +
    xlab("log2 Fold Change") +
    ggtitle(title) + 
    theme_light() + 
    theme(legend.position = "none")  +
    geom_vline(xintercept = log2(2), linetype='dashed', colour='grey') + 
    geom_vline(xintercept = log2(0.5), linetype='dashed', colour='grey') + 
    geom_hline(yintercept = -log10(0.05), linetype='dashed', colour='grey')# + 
  #geom_text_repel(data = subset(data, data[['padj']] < 0.01 & abs(data[['log2FoldChange']]) > 2), aes(label = subset(data, data[['padj']] < 0.01 & abs(data[['log2FoldChange']]) > 2)[["labels"]], colour='black'))
  
  print(plot)
}

#### need to loop through soecies as well as batches 
diff_exp = function(metadata, counts_meta, counts){
  results_list = list()
  nsig_list = list()
  for (dataset in unique(metadata$batch)){
    if (dataset == 5){
      next
    }
    meta = metadata %>% filter(batch == dataset)
    counts2 = counts[, metadata$batch == dataset]
  
    resistants = unique(meta[meta$resistance == 'resistant']$condition)
    susceptibles = unique(meta[meta$resistance == 'susceptible']$condition)
    comparisons = crossing(resistants, susceptibles)
    print(dataset)
    print(as.data.frame(comparisons))

    for (i in 1:nrow(comparisons)){
      res = comparisons[i, 'resistants']
      sus = comparisons[i, 'susceptibles']
      comp = glue("{res}_v_{sus}")
      print(comp)
    
      controls = which(meta$condition %in% sus)
      cases = which(meta$condition %in% res)
    
      idxs = c(controls, cases)
      subcounts = counts2[, idxs]
      subsamples = meta[idxs,]
      
      # make treatment a factor with the 'susceptible' as reference
      subsamples$treatment = as.factor(subsamples$resistance)
      subsamples$treatment = relevel(subsamples$treatment, "susceptible")
      # make DESeq dataset
      print("subcounts shape")
      print(dim(subcounts))
      print("subsamples shape")
      print(dim(subsamples))
      dds = DESeqDataSetFromMatrix(countData = subcounts, 
                                   colData = subsamples, 
                                   design = ~ treatment)
      ###### estimate paramters and normalise 
      dds = estimateSizeFactors(dds)
      dds = estimateDispersions(dds)
      dds = estimateDispersions(dds)
      cds = nbinomWaldTest(dds)
      results = results(cds, contrast = c("treatment", "susceptible", "resistant")) %>% as.data.frame() 
      results = results[order(results$padj),] #order by pvalue 
      results = results %>% rownames_to_column("GeneID") %>% dplyr::mutate("FC" = (2^log2FoldChange))
      
      ### absolute difference
      #### Get rowsums of counts, grouping by case/control. Then get difference of counts and join with DE results
      readdiff = data.frame(t(rowsum(t(subcounts), group = subsamples$treatment, na.rm = T))) #transpose and get rowsums for each group
      readdiff$absolute_diff = readdiff[,"resistant"] - readdiff[,"susceptible"] #get difference
      readdiff = data.frame(readdiff) %>% rownames_to_column('GeneID')
      results = unique(left_join(results, readdiff[,c('GeneID','absolute_diff')]))
      
      # join DE results with normal gene names
      results = unique(left_join(results, gene_names))
      results_list[[comp]] = results

      fwrite(results, glue("results/genediff/{comp}_diffexp.csv")) #write to csv 
      # volcano plot for each comparison, First make vector of labels which is AGAPs unless a gene name exists
      results$labels = results %>% dplyr::mutate("Gene_name" = case_when(GeneName == "" ~ GeneID,
                                                                         is.na(GeneName) ~ GeneID,
                                                                         TRUE ~ GeneName)) %>% select(Gene_name) %>% deframe()
      #get number of sig genes 
      res1 = results %>% filter(padj < 0.05) %>% 
        count("direction" = FC > 1) %>% 
        dplyr::mutate("direction" = case_when(direction == FALSE ~ "Downregulated, padj = 0.05",
                                              direction == TRUE ~ "Upregulated, padj = 0.05")) %>%
        dplyr::rename(!!glue("{comp}_ngenes") := "n")
      
      res2 = results %>% filter(padj < 0.001) %>% 
        count("direction" = FC > 1) %>% 
        dplyr::mutate("direction" = case_when(direction == FALSE ~ "Downregulated, padj = 0.001",
                                              direction == TRUE ~ "Upregulated, padj = 0.001")) %>%
        dplyr::rename(!!glue("{comp}_ngenes") := "n")
      
      nsig_list[[comp]] = bind_rows(res1, res2)

      pdf(glue("results/genediff/{comp}_Volcano_plot.pdf"))
      volcano(data = results, title = comp)
      null = dev.off()
      cat("\n", glue("{comp} complete!"), "\n")
      }
  }
  return(list(results_list, nsig_list))
}



res_list = diff_exp(metadata, counts_meta, counts)
results_list = res_list[[1]]
nsig_list = res_list[[2]]

#### write to excel file on diff sheets #### 
sheets = names(results_list)
wb <- createWorkbook("Workbook")
for (i in 1:length(sheets)){
  addWorksheet(wb, glue("{sheets[[i]]}"))
  writeData(wb, sheets[i], results_list[[i]], rowNames = FALSE, colNames = TRUE)
}
#### save workbook to disk once all worksheets and data have been added ####
saveWorkbook(wb,file="results/genediff/meta_genediff.xlsx", overwrite = TRUE)

# Join different comparisons together and write out number of sig genes 
purrr::reduce(nsig_list, inner_join) %>% fwrite("results/genediff/nsig_genes.tsv", sep="\t", col.names = TRUE)



fc_data = data.frame("GeneID" = results_list$Tiefora_v_Ngousso$GeneID)
pval_data = data.frame("GeneID" = results_list$Tiefora_v_Ngousso$GeneID)
for (i in 1:length(results_list)){
  name = sheets[i]
  name_var = glue("{name}_log2FoldChange")
  name_pval = glue("{name}_padj")
  df = results_list[[i]] %>% 
    select(c("GeneID", "log2FoldChange")) %>% 
    rename({{ name_var }} := log2FoldChange)
  
  pval_df = results_list[[i]] %>% 
    select(c("GeneID", "padj")) %>% 
    rename({{ name_var }} := padj)
  
  
  fc_data = fc_data %>% inner_join(df) %>% distinct()
  pval_data = pval_data %>% inner_join(pval_df) %>% distinct()
}

round_df = function(df, digits) {
  
  #' This function rounds all the numeric columns of a data.frame
  nums = vapply(df, is.numeric, FUN.VALUE = logical(1))
  df[,nums] = round(df[,nums], digits = digits)
  (df)
}
#fc_data['mean'] = (apply(fc_data, 1, mean))
#fc_data['median'] = (apply(fc_data, 1, median))

fc_data = fc_data %>%  inner_join(names_df)
pval_data = pval_data  %>%  inner_join(names_df)

pval_data %>% select(-c(GeneID2, TranscriptID)) %>% round_df(3) %>% distinct() %>% fwrite(., file="results/pval_data.tsv", sep="\t")
fc_data %>% select(-c(GeneID2, TranscriptID)) %>% round_df(2) %>% distinct() %>% fwrite(., file="results/fc_data.tsv", sep="\t")


