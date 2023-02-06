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


f = glue
datasets = c("gamb_colu", "gamb_colu_arab", "gamb_colu_arab_fun")

counts = fread("results/counts.gamb_colu.tsv") %>% as.data.frame()
AGAMnames_df = fread("../rna-seq-busia/resources/exampleGene2TranscriptMap.tsv", sep="\t") %>% distinct()
#names_df = data.frame("GeneID2" = make.unique(counts_meta$GeneID))
#names_df$GeneID = gsub("\\..*","",names_df$GeneID2)
#names_df = names_df %>% left_join(gene_names)

AFUNnames_df = fread("resources/AfunGenes2TranscriptMap.tsv", sep="\t") %>% distinct()


### log2 counts 
for (analysis in c("gamb_colu", "gamb_colu_arab", "gamb_colu_arab_fun", "fun")){
  counts = fread(f("results/counts.{analysis}.tsv")) %>% 
  as.data.frame() %>% 
  column_to_rownames("GeneID") %>% 
  mutate_if(is.numeric, as.integer)

  metadata = fread("config/sample_metadata.tsv")
  metadata$batch = as.factor(metadata$batch)
  
  if (analysis == 'gamb_colu'){
    sp_bool = metadata$species %in% c("gambiae", "coluzzii")
  } else if (analysis == 'gamb_colu_arab'){
    sp_bool = metadata$species %in% c("gambiae", "coluzzii", "arabiensis")
  } else if (analysis == 'gamb_colu_arab_fun'){
    sp_bool = metadata$species %in% c("gambiae", "coluzzii", "arabiensis", "funestus")
  } else if (analysis == 'fun'){
    sp_bool = metadata$species == "funestus"
  }
    # subset to analysis
  meta = metadata[sp_bool, ]
  
  dds = DESeqDataSetFromMatrix(countData = counts, 
                               colData = meta, 
                               design = ~ resistance)
  dds = estimateSizeFactors(dds)
  dds = estimateDispersions(dds)
  norm_counts = round_df(log2(as.data.frame(counts(dds, normalized=TRUE)) + 1), 2)  %>% 
    rownames_to_column("GeneID")
  
  print(dim(norm_counts))
  fwrite(norm_counts, f("results/log2counts.{analysis}.tsv"), sep="\t")
}


### PCA ###
# make DESeq dataset

counts = fread("results/counts.gamb_colu_arab_fun.tsv") %>% 
  as.data.frame() %>% 
  column_to_rownames("GeneID") %>% 
  mutate_if(is.numeric, as.integer)


metadata = fread("config/sample_metadata.tsv")
metadata$batch = as.factor(metadata$batch)

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

round_df = function(df, digits) {
  #' This function rounds all the numeric columns of a data.frame
  nums = vapply(df, is.numeric, FUN.VALUE = logical(1))
  df[,nums] = round(df[,nums], digits = digits)
  (df)
}

#### need to loop through species as well as batches 
diff_exp = function(dataset, names_df){
  results_list = list()
  nsig_list = list()
  
  metadata = fread("config/sample_metadata.tsv") %>% rename("sampleID"="colData") %>% as.data.frame()
  metadata[(metadata$batch == 1) &(metadata$species == 'arabiensis'), 'batch'] = 0
  metadata$batch = factor(metadata$batch)

  # load dataset 
  counts = fread(f("results/counts.{dataset}.tsv"), sep="\t") %>% 
    as.data.frame() %>%
    column_to_rownames("GeneID") %>% 
    mutate_if(is.numeric, as.integer)
  
  if (dataset == 'gamb_colu'){
    sp_bool = metadata$species %in% c("gambiae", "coluzzii")
  } else if (dataset == 'gamb_colu_arab'){
    sp_bool = metadata$species %in% c("gambiae", "coluzzii", "arabiensis")
  } else if (dataset == 'gamb_colu_arab_fun'){
    sp_bool = metadata$species %in% c("gambiae", "coluzzii", "arabiensis", "funestus")
  } else if (dataset == 'fun'){
    sp_bool = metadata$species == "funestus"
  }
  
  print(dataset)
  # subset to dataset
  meta = metadata[sp_bool, ]
  print(dim(meta))
  print(dim(counts))
  
  # analyse each experiment separately 
  for (experiment in unique(meta$batch)){
    if (experiment == 5){
      next
    }
    
    stopifnot(nrow(meta) == length(counts))
    
    # subset to batch 
    meta2 = meta %>% filter(batch == experiment)
    counts2 = counts[, meta2$sampleID]
    
    stopifnot(all(meta2$sampleID == colnames(counts2)))
    
    resistants = unique(meta2[meta2$resistance == 'resistant',]$condition)
    susceptibles = unique(meta2[meta2$resistance == 'susceptible',]$condition)
    comparisons = crossing(resistants, susceptibles)
    print(experiment)
    print(as.data.frame(comparisons))

    for (i in 1:nrow(comparisons)){
      res = comparisons[i, 'resistants']
      sus = comparisons[i, 'susceptibles']
      comp = glue("{res}_v_{sus}")
      print(comp)
    
      controls = which(meta2$condition %in% sus)
      cases = which(meta2$condition %in% res)
    
      idxs = c(controls, cases)
      subcounts = counts2[, idxs]
      subsamples = meta2[idxs,]
      
      # make treatment a factor with the 'susceptible' as reference
      subsamples$treatment = as.factor(subsamples$resistance)
      subsamples$treatment = relevel(subsamples$treatment, "susceptible")
      # make DESeq dataset
      print("subcounts shape")
      print(dim(subcounts))
      print(head(subcounts))
      print("subsamples shape")
      print(dim(subsamples))
      print(head(subsamples))
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
      results = results %>% mutate(log2FoldChange=log2FoldChange*-1)
      results = results %>% rownames_to_column("GeneID") %>% dplyr::mutate("FC" = (2^log2FoldChange))
      
      ### absolute difference
      #### Get rowsums of counts, grouping by case/control. Then get difference of counts and join with DE results
      readdiff = data.frame(t(rowsum(t(subcounts), group = subsamples$treatment, na.rm = T))) #transpose and get rowsums for each group
      readdiff$absolute_diff = readdiff[,"resistant"] - readdiff[,"susceptible"] #get difference
      readdiff = data.frame(readdiff) %>% rownames_to_column('GeneID')
      results = unique(left_join(results, readdiff[,c('GeneID','absolute_diff')]))
      
      # join DE results with normal gene names
      results = unique(left_join(results, names_df))
      results_list[[comp]] = results

      fwrite(results, glue("results/genediff/{comp}_diffexp.csv")) #write to csv 
      # volcano plot for each comparison, First make vector of labels which is AGAPs unless a gene name exists
      #results$labels = results %>% dplyr::mutate("Gene_name" = case_when(GeneName == "" ~ GeneID,
      #                                                                   is.na(GeneName) ~ GeneID,
      #                                                                   TRUE ~ GeneName)) %>% select(Gene_name) %>% deframe()
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

      #pdf(glue("results/genediff/{comp}_Volcano_plot.pdf"))
      #volcano(data = results, title = comp)
      #null = dev.off()
      cat("\n", glue("{comp} complete!"), "\n")
      }
  }
  #### write to excel file on diff sheets #### 
  sheets = names(results_list)
  wb <- createWorkbook("Workbook")
  for (i in 1:length(sheets)){
    addWorksheet(wb, glue("{sheets[[i]]}"))
    writeData(wb, sheets[i], results_list[[i]], rowNames = FALSE, colNames = TRUE)
  }
  #### save workbook to disk once all worksheets and data have been added ####
  saveWorkbook(wb,file=f("results/genediff/{dataset}_genediff.xlsx"), overwrite = TRUE)
  # Join different comparisons together and write out number of sig genes 
  purrr::reduce(nsig_list, inner_join) %>% fwrite(f("results/genediff/{dataset}_nsig_genes.tsv"), sep="\t", col.names = TRUE)
  
  fc_data = data.frame("GeneID" = results_list[[1]]$GeneID)
  pval_data = data.frame("GeneID" = results_list[[1]]$GeneID)
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
  
  fc_data = fc_data %>% inner_join(names_df)
  pval_data = pval_data %>% inner_join(names_df)
  
  pval_data %>% 
    select(-TranscriptID) %>% 
    round_df(3) %>% 
    distinct() %>% 
    fwrite(., file=f("results/pvals.{dataset}.tsv"), sep="\t")
  
  fc_data %>% 
    select(-TranscriptID) %>%
    round_df(2) %>% 
    distinct() %>% 
    fwrite(., file=f("results/fcs.{dataset}.tsv"), sep="\t")

  return(list(results_list, nsig_list))
}



for (dataset in datasets){
  res_list = diff_exp(dataset, names_df = AGAMnames_df)
}


res_list = diff_exp("gamb_colu", names_df = AGAMnames_df)


res = diff_exp(dataset = "fun", names_df = AFUNnames_df)


all(metadata$sampleID == colnames(counts))

colnames(counts)
