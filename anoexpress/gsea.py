import pandas as pd
import numpy as np

from .candidates import load_candidates
from .utils import resolve_gene_id


def load_genes_for_enrichment(analysis, func, gene_ids, percentile, microarray, low_count_filter=None, gff_method='malariagen_data'):
   
    assert func is not None or gene_ids is not None, "either a ranking function (func) or gene_ids must be provided"
    assert func is None or gene_ids is None, "Only a ranking function (func) or gene_ids must be provided, not both"

    fc_data = pd.read_csv(f"https://raw.githubusercontent.com/sanjaynagi/AnoExpress/main/results/fcs.{analysis}.tsv", sep="\t")
    fc_genes = fc_data.reset_index()['GeneID'].to_list()
    name = 'enrich'

    if func:
      # get top % percentile genes ranked by func
      fc_ranked = load_candidates(analysis=analysis, name='enrich', func=func, microarray=microarray, low_count_filter=low_count_filter)
      percentile_idx = fc_ranked.reset_index()['GeneID'].unique().shape[0] * percentile
      top_geneIDs = fc_ranked.reset_index().loc[:, 'GeneID'][:int(percentile_idx)] 
    elif gene_ids:
      top_geneIDs = resolve_gene_id(gene_id=gene_ids, analysis=analysis, gff_method=gff_method)

    return top_geneIDs, fc_genes

def go_hypergeometric(analysis, func=None, gene_ids=None, percentile=0.05, microarray=False, low_count_filter=None):
    """
    Perform a hypergeometric test on GO terms of the the top % percentile genes ranked by user input function, or on 
    a user inputted gene_id list

    Parameters
    ----------
    analysis: {"gamb_colu", "gamb_colu_arab", "gamb_colu_arab_fun", "fun"}
      which analysis to load gene expression data for. analyses with more species will have less genes
      present, due to the process of finding orthologs.
    func: function
      function to rank genes by (such as np.nanmedian, np.nanmean)
    gene_ids: list, optional
      list of gene ids to perform hypergeometric test on. Defaults to None
    percentile: float, optional
      percentile of genes to use for the enriched set in hypergeometric test. Defaults to 0.05
    microarray: bool, optional
      whether to include the IR-Tex microarray data in the gene ranking. Default is False.
    low_count_filter: int, optional
      if provided, genes with a median count below the threshold will be removed before gene ranking. Default is None.

    Returns
    -------
    go_hypergeo_results: pd.DataFrame
    """

    top_geneIDs, fc_genes = load_genes_for_enrichment(analysis=analysis, func=func, gene_ids=gene_ids, percentile=percentile, microarray=microarray, low_count_filter=low_count_filter)

    # load gene annotation file 
    gaf_df = pd.read_csv("https://raw.githubusercontent.com/sanjaynagi/AnoExpress/main/resources/AgamP4.gaf", sep="\t")
    go_annotations = gaf_df[['go_term', 'descriptions']].rename(columns={'go_term':'annotation'}).drop_duplicates()
    gaf_df = gaf_df[['GeneID', 'go_term']].drop_duplicates()
    gaf_df = gaf_df.query("GeneID in @fc_genes")
    N = gaf_df.GeneID.unique().shape[0] #Total number of genes with some annotation 
    k = np.isin(gaf_df.loc[:, 'GeneID'].unique(), top_geneIDs).sum() 

    hyper_geo = _hypergeometric(
        annotation_df=gaf_df, 
        column_name='go_term', 
        target_gene_list=top_geneIDs,
        N=N,
        k=k)    
    hyper_geo = hyper_geo.merge(go_annotations, how='left')
    return(hyper_geo)


def pfam_hypergeometric(analysis, func=None, gene_ids=None, percentile=0.05, microarray=False, low_count_filter=None):
    """
    Perform a hypergeometric test on PFAM domains of the the top % percentile genes ranked by user input function,
    or on a user inputted gene_id list

    Parameters
    ----------
    analysis: {"gamb_colu", "gamb_colu_arab", "gamb_colu_arab_fun", "fun"}
      which analysis to load gene expression data for. analyses with more species will have less genes
      present, due to the process of finding orthologs.
    name: str
      name of the function to rank genes by
    func: function
      function to rank genes by (such as np.nanmedian, np.nanmean)
    gene_ids: list, optional
      list of gene ids to perform hypergeometric test on. Defaults to None
    percentile: float, optional
      percentile of genes to use for the enriched set in hypergeometric test. Defaults to 0.05
    microarray: bool, optional
      whether to include the IR-Tex microarray data in the gene ranking. Default is False.
    low_count_filter: int, optional
      if provided, genes with a median count below the threshold will be removed before gene ranking. Default is None.
    
    Returns
    -------
    pfam_hypergeo_results: pd.DataFrame
    """

    top_geneIDs, fc_genes = load_genes_for_enrichment(analysis=analysis, func=func, gene_ids=gene_ids, percentile=percentile, microarray=microarray, low_count_filter=low_count_filter)

    # load gene annotation file 
    pfam_df = pd.read_csv("https://github.com/sanjaynagi/AnoExpress/blob/main/resources/Anogam_long.pep_Pfamscan.seqs.gz?raw=true", sep="\s+", header=None, compression='gzip').iloc[:, [0,4]]
    pfam_df.loc[:, 0] = pfam_df.loc[:, 0].str.replace("Anogam_", "").str.replace("-R[A-Z]", "", regex=True)
    pfam_df.columns = ['GeneID', 'pfam']
    pfam_df = pfam_df.query("GeneID in @fc_genes")
    N = pfam_df.GeneID.unique().shape[0] #Total number of genes with some annotation 
    k = np.isin(pfam_df.loc[:, 'GeneID'].unique(), top_geneIDs).sum()  

    # run hypergeometric test
    hyper_geo = _hypergeometric(
        annotation_df=pfam_df, 
        column_name='pfam', 
        target_gene_list=top_geneIDs,
        N=N,
        k=k)
        
    return(hyper_geo)

def kegg_hypergeometric(analysis, func=None, gene_ids=None, percentile=0.05, microarray=False, low_count_filter=None):
    """
    Perform a hypergeometric test on GO terms of the the top % percentile genes ranked by user input function.

    Parameters
    ----------
    analysis: {"gamb_colu", "gamb_colu_arab", "gamb_colu_arab_fun", "fun"}
      which analysis to load gene expression data for. analyses with more species will have less genes
      present, due to the process of finding orthologs.
    name: str
      name of the function to rank genes by
    func: function
      function to rank genes by (such as np.nanmedian, np.nanmean)
    gene_ids: list, optional
      list of gene ids to perform hypergeometric test on. Defaults to None
    percentile: float, optional
      percentile of genes to use for the enriched set in hypergeometric test. Defaults to 0.05
    microarray: bool, optional
      whether to include the IR-Tex microarray data in the gene ranking. Default is False.
    low_count_filter: int, optional
      if provided, genes with a median count below the threshold will be removed before gene ranking. Default is None.

    Returns
    -------
    go_hypergeo_results: pd.DataFrame
    """

    top_geneIDs, fc_genes = load_genes_for_enrichment(analysis=analysis, func=func, gene_ids=gene_ids, percentile=percentile, microarray=microarray, low_count_filter=low_count_filter)

    # load gene annotation file 
    kegg_df = pd.read_csv("https://raw.githubusercontent.com/sanjaynagi/AnoExpress/main/resources/AgamP4.kegg", sep="\t")
    kegg_annotations = kegg_df[['kegg_pathway', 'description']].rename(columns={'kegg_pathway':'annotation'}).drop_duplicates()
    kegg_df = kegg_df[['GeneID', 'kegg_pathway']].drop_duplicates()
    kegg_df = kegg_df.query("GeneID in @fc_genes")
    N = kegg_df.GeneID.unique().shape[0] #Total number of genes with some annotation 
    k = np.isin(kegg_df.loc[:, 'GeneID'].unique(), top_geneIDs).sum() 

    hyper_geo = _hypergeometric(
        annotation_df=kegg_df, 
        column_name='kegg_pathway', 
        target_gene_list=top_geneIDs,
        N=N,
        k=k)    
    hyper_geo = hyper_geo.merge(kegg_annotations, how='left')
    return(hyper_geo)

def _hypergeometric(annotation_df, column_name, target_gene_list, N, k):
    """
    This function performs a hypergeometric test on a given annotation column
    """
    from scipy.stats import hypergeom
    from tqdm import tqdm
    from statsmodels.stats.multitest import fdrcorrection

    # get unique annotations
    unique_annots = annotation_df.loc[:, column_name].unique()

    sig_list = []
    res_list = []
    for annot in tqdm(unique_annots):

        annot_genes = annotation_df.query("{col} == @annot".format(col=column_name))['GeneID']
        m = len(annot_genes)

        x = annot_genes.isin(target_gene_list).sum()
        res = hypergeom(M=N, 
                        n=m, 
                        N=k).sf(x-1)
        sig_list.append(annot)
        res_list.append(res)    

    hyper_geo =  pd.DataFrame({'annotation': sig_list, 'pval':res_list})
    hypo, hyper_geo.loc[:, 'padj'] = fdrcorrection(hyper_geo['pval'])           
    return(hyper_geo.sort_values(by='pval'))

