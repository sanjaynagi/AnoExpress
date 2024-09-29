import numpy as np
import pandas as pd

from .data import data, load_annotations
from .utils import _gene_ids_from_annotation


def load_candidates(analysis, name='median', func=np.nanmedian, query_annotation=None, query_fc=None, microarray=False, low_count_filter=None, fraction_na_allowed=None):
    """
    Load the candidate genes for a given analysis. Optionally, filter by annotation or fold change data.
    
    Parameters
    ----------
    analysis: {"gamb_colu", "gamb_colu_arab", "gamb_colu_arab_fun", "fun"}
      which analysis to load gene expression data for. analyses with more species will have less genes
      present, due to the process of finding orthologs.
    name: str, optional
      name of the function to rank genes by. Defaults to 'median'
    func: function, optional
      function to rank genes by. Defaults to np.nanmedian
    query_annotation: str or list, optional
      filter genes by GO or PFAM annotation. Defaults to None
    query_fc: float, optional
      filter genes by fold change. Defaults to None
    microarray: bool, optional
      whether to include the IR-Tex microarray data in the requested data. Default is False.
    low_count_filter: int, optional
      if provided, genes with a median count below the threshold will be removed before gene ranking. Default is None.
    fraction_nas_allowed: float, optional
      fraction of missing values allowed in the data. Defaults to 0.5
    
    Returns
    -------
    fc_ranked: pd.DataFrame
    """
    
    fc_data = data(data_type='fcs', analysis=analysis, microarray=microarray, annotations=True, sort_by=None, low_count_filter=low_count_filter, fraction_na_allowed=fraction_na_allowed)

    if query_annotation is not None:
      gene_annot_df = load_annotations()
      gene_ids = _gene_ids_from_annotation(gene_annot_df=gene_annot_df, annotation=query_annotation)
      fc_data = fc_data.query("GeneID in @gene_ids")
      assert not fc_data.empty, "No genes were found for the selection. It is possible these genes were removed by the ortholog finding process"
    
    fc_ranked = fc_data.apply(func, axis=1).to_frame().rename(columns={0:f'{name} log2 Fold Change'}).copy()
    fc_ranked = fc_ranked.sort_values(f'{name} log2 Fold Change', ascending=False)
    fc_ranked = fc_ranked.reset_index()
    fc_ranked.loc[:, f'{name} Fold Change'] = np.round(2**fc_ranked.loc[:, f'{name} log2 Fold Change'], 2)

    if query_fc is not None:
       fc_ranked = fc_ranked.query(f'`{name} Fold Change` > {query_fc}')

    return(fc_ranked)





def consistent_genes(analysis, direction, n_experiments, low_count_filter=None):
    
    fc_data = data(data_type="fcs", analysis=analysis, microarray=False, annotations=True, low_count_filter=low_count_filter)
    pval_data = data(data_type='pvals', analysis=analysis, microarray=False, annotations=True, low_count_filter=low_count_filter)

    print(f"There are {fc_data.shape[0]} genes and {fc_data.shape[1]} differential expression comparisons in {analysis}")
    if direction == 'up':
        res_df = fc_data[fc_data.apply(lambda x: (x > 0).sum() >= n_experiments , axis=1)]
        res_pval = pval_data[pval_data.apply(lambda x: x < 0.05, axis=1).sum(axis=1) > n_experiments].reset_index()['GeneID'].to_list()
        res_df = res_df.query("GeneID in @res_pval")

        if res_df.empty: 
            print(f"There are no genes expressed direction={direction} in {n_experiments} experiments")
            return
        else:
            return(res_df)
    else: 
        res_df = fc_data[fc_data.apply(lambda x: (x < 0).sum() >= n_experiments, axis=1)]
        res_pval = pval_data[pval_data.apply(lambda x: x < 0.05, axis=1).sum(axis=1) > n_experiments].reset_index()['GeneID'].to_list()
        res_df = res_df.query("GeneID in @res_pval")
    if res_df.empty:
        print(f"There are no genes expressed {direction} in {n_experiments} experiments")
        return
    return(res_df)



def contig_expression(contig, analysis, data_type='fcs', microarray=False, pvalue_filter=None, size=10, step=5, fraction_na_allowed=None):
    """
    Calculate fold change data for a given contig. Returns both raw fold change data for each experiment and a moving average
    
    
    Parameters
    ----------
    contig: str
      contig, one of 2L, 2R, 3L, 3R, X or 2RL, 3RL
    analysis: {"gamb_colu", "gamb_colu_arab", "gamb_colu_arab_fun", "fun"}
      which analysis to load gene expression data for. analyses with more species will have less genes
      present, due to the process of finding orthologs.
    data_type: {"fcs", "log2counts"}, optional
      whether to load fold change data or log2 counts data. Defaults to 'fcs'
    microarray: bool, optional
      whether to include the IR-Tex microarray data in the plot
    pvalue_filter: float, optional
      filter genes by pvalue. Defaults to None
    size: int, optional
      size of window in genes for moving average. Defaults to 10
    step: int, optional
      step size in genes for moving average. Defaults to 5
    fraction_na_allowed: float, optional
      fraction of missing values allowed for each gene in the data. Defaults to no filter.
    """
    import malariagen_data
    import allel
    
    fc_data = data(
                    data_type=data_type, 
                    analysis=analysis, 
                    microarray=microarray, 
                    annotations=True, 
                    pvalue_filter=pvalue_filter,
                    fraction_na_allowed=fraction_na_allowed
                    ).reset_index()
    
    ag3 = malariagen_data.Ag3()
    gff = ag3.genome_features(contig).query("type == 'gene'").assign(midpoint = lambda x: (x.start + x.end)/2)
    genes_df = gff[['ID', 'midpoint']].rename(columns={'ID':'GeneID'})
    genes_df = genes_df.merge(fc_data, how='left').set_index(['GeneID', 'GeneName', 'GeneDescription', 'midpoint'])
    # calculate moving average of fold changes
    median_fc = genes_df.mean(axis=1).to_frame().rename(columns={0:'median_fc'}).reset_index()
    moving_mid = allel.moving_statistic(median_fc.midpoint, np.nanmedian, size=size, step=step)
    moving_med = allel.moving_statistic(median_fc['median_fc'], np.nanmedian, size=size, step=step)
    windowed_fold_change_df = pd.DataFrame({'midpoint':moving_mid, 'median_fc':moving_med})

    fold_change_df = genes_df.reset_index().melt(id_vars=['midpoint', 'GeneID', 'GeneName', 'GeneDescription'], var_name='comparison', value_name='fold_change')  
    return fold_change_df, windowed_fold_change_df

