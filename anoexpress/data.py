import pandas as pd
from tqdm import tqdm
import numpy as np

from .utils import (
    resolve_gene_id, 
    _gene_ids_from_annotation,
    filter_low_counts,
    load_gff
)

index_col = {'fcs':'comparison',
            'pvals': 'comparison',
            'log2counts': 'sampleID'}


taxon_query_dict = { 'gamb_colu':"species in ['gambiae','coluzzii']", 
                  'gamb_colu_arab':"species in ['gambiae','coluzzii','arabiensis']", 
                  'gamb_colu_arab_fun':"species in ['gambiae','coluzzii','arabiensis','funestus']", 
                  'fun':"species in ['funestus']"
                    }


def load_results_arrays(data_type, analysis):
    """
    Load the counts data for a given analysis and sample query
    """
    assert data_type in ['log2counts', 'fcs', 'pvals'], "data_type must be either 'log2counts', 'fcs' or 'pvals'"
    assert analysis in ['gamb_colu', 'gamb_colu_arab', 'gamb_colu_arab_fun', 'fun', 'irtex'], "analysis must be either 'gamb_colu', 'gamb_colu_arab', 'gamb_colu_arab_fun', 'fun' or 'irtex'"
    df_data = pd.read_csv(f"https://raw.githubusercontent.com/sanjaynagi/AnoExpress/main/results/{data_type}.{analysis}.tsv", sep="\t")
    df_data = df_data.set_index('GeneID')
    return(df_data)


def sample_metadata(analysis):
    """
    Load the sample metadata in a pandas dataframe
    """
    sample_metadata = pd.read_csv("https://raw.githubusercontent.com/sanjaynagi/AnoExpress/main/config/sample_metadata.tsv", sep="\t")
    sample_metadata = sample_metadata.query(taxon_query_dict[analysis])
    return sample_metadata


def xpress_metadata():
    comparison_metadata = pd.read_csv("https://raw.githubusercontent.com/sanjaynagi/AnoExpress/main/config/comparison_metadata.tsv", sep="\t")
    return(comparison_metadata)


def irtex_metadata():
    """
    Load the ir-tex fold change data for a given analysis and sample query
    """   
    comparison_metadata = pd.read_csv("https://raw.githubusercontent.com/sanjaynagi/AnoExpress/main/config/irtex_metadata.tsv", sep="\t")
    return(comparison_metadata)


def metadata(analysis, microarray=False):
    """
    Load the comparisons metadata from both AnoExpress and IR-Tex in a pandas dataframe

    Parameters
    ----------
    analysis: {"gamb_colu", "gamb_colu_arab", "gamb_colu_arab_fun", "fun"}
      which analysis to load gene expression data for. analyses with more species will have less genes
      present, due to the process of finding orthologs.
    microarray: bool, optional
      whether to include the IR-Tex microarray metadata in the metadata. Default is False.
    sample_metadata: bool, optional
      whether to return the sample metadata in addition to the comparisons metadata. Default is False.

    Returns
    -------
    comparisons_df: pandas dataframe
    samples_df: pandas dataframe
    """
    # load metadata from AnoExpress
    metadata = xpress_metadata()   
    metadata = metadata.assign(technology='rnaseq')

    # load metadata from IR-Tex
    if microarray == True:
        irtex_meta = irtex_metadata()
        irtex_meta = irtex_meta.assign(technology='microarray')
        metadata = pd.concat([metadata, irtex_meta]).reset_index(drop=True)

    # subset to the species of interest
    metadata = metadata.query(taxon_query_dict[analysis])

    return metadata



def data(data_type, analysis, microarray=False, gene_id=None, sample_query=None, sort_by=None, annotations=False, pvalue_filter=None, low_count_filter=None, fraction_na_allowed=None):
    """
    Load the combined data for a given analysis and sample query

    Parameters
    ----------
    data_type: {"fcs", "pvals", "log2counts"}
      which data type to load gene expression data for. log2counts are the log2 counts, fcs are the fold changes, and pvals are the adjusted p-values.
    analysis: {"gamb_colu", "gamb_colu_arab", "gamb_colu_arab_fun", "fun"}
      which analysis to load gene expression data for. analyses with more species will have less genes
      present, due to the process of finding orthologs.
    microarray: bool, optional
      whether to include the IR-Tex microarray data in the requested data. Default is False.
    gene_id: str or list, optional
      A string (AGAP/AFUN identifier or genomic span in the format 2L:500-10000), or list of strings, or path to a file containing a list of gene ids in the first column. 
      Input file can be .tsv, .txt, or .csv, or .xlsx.
    sample_query: str, optional
      A string containing a pandas query to subset the samples of interest from the comparison metadata file. For example, to plot only the
      samples from Burkina Faso, use "country == 'Burkina Faso'". Defaults to None.
    sort_by: {"median", "mean", "agap", "position", None}, optional
      sort by median/mean of fold changes (descending), or by AGAP, or by position in the genome, or dont sort input gene ids.
    annotations: bool, optional
      whether to add gene name and description to the dataframe as index. Default is False.
    pvalue_filter: float, optional
      if provided, fold-change entries with an adjusted p-value below the threshold will be set to NaN. Default is None.
      ignored if the data_type is not 'fcs'. 
    low_count_filter: int, optional
      if provided, genes with a median count below the threshold will be removed from the dataframe. Default is None.
    fraction_na_allowed: float, optional
      fraction of missing values allowed in the data. Defaults to 0.5
    
    Returns
    -------
    results_data: pandas dataframe
    """
    assert analysis in ['gamb_colu', 'gamb_colu_arab', 'gamb_colu_arab_fun', 'fun'], "analysis must be one of 'gamb_colu', 'gamb_colu_arab', 'gamb_colu_arab_fun', 'fun'"
    assert data_type in ['log2counts', 'fcs', 'pvals'], "data_type must be one of 'log2counts', 'fcs', 'pvals'"
    # load the metadata and subset to the species of interest
    if data_type in ['fcs', 'pvals']:
      df_metadata = metadata(analysis=analysis, microarray=microarray)
    else:
      # if count data load sample level metadata 
      df_metadata = sample_metadata(analysis=analysis)

    # load the data and merge wit microarray data if requested
    df = load_results_arrays(data_type=data_type, analysis=analysis)
    # load ir-tex data and merge
    if microarray == True and data_type != 'log2counts':
        irtex_df = load_results_arrays(data_type=data_type, analysis='irtex')
        df = df.reset_index().merge(irtex_df, on='GeneID', how='left').set_index('GeneID')

    # get the sample or comparison ids 
    metadata_ids = df_metadata.loc[:, index_col[data_type]].to_list()
    # subset to the species comparisons of interest
    df = df.loc[:, metadata_ids]

    if sample_query:
      # subset to the sample ids of interest
      mask = df_metadata.eval(sample_query).to_list()
      df = df.loc[:, mask]

    # subset to the gene ids of interest including reading file 
    if gene_id is not None:
      gene_id = resolve_gene_id(gene_id=gene_id, analysis=analysis)
      df = df.query("GeneID in @gene_id")

    if annotations: # add gene name and description to the dataframe as index 
      df = add_annotations_to_array(df)
    
    if data_type == 'fcs' and pvalue_filter is not None:
      pval_df = data(data_type="pvals", analysis=analysis, microarray=microarray, gene_id=gene_id, sort_by=sort_by, annotations=False)
      df = null_fold_changes(pval_df=pval_df, fc_df=df, threshold=pvalue_filter)

    # sort genes 
    df = _sort_genes(df=df, analysis=analysis, sort_by=sort_by)

    # remove low count genes
    if low_count_filter is not None:
      df = filter_low_counts(data_df=df, data_type=data_type, analysis=analysis, gene_id=gene_id, count_threshold=low_count_filter, func=np.nanmedian)
      
    # remove genes with lots of NA
    if fraction_na_allowed:
      df = filter_nas(df=df, fraction_na_allowed=fraction_na_allowed)
    
    return df

def null_fold_changes(pval_df, fc_df, threshold=0.05):
    fold_changes_null = fc_df.copy() # make a copy of fold_changes to modify
    fold_changes_null[pval_df > threshold] = np.nan # set values to NaN where pvalue > 0.05
    return fold_changes_null
 
def add_annotations_to_array(df):
    df_annots = pd.read_csv("https://raw.githubusercontent.com/sanjaynagi/AnoExpress/main/resources/AgamP4.annots.tsv", sep="\t")
    df = df.reset_index().merge(df_annots, on="GeneID", how="left").set_index(["GeneID", "GeneName", "GeneDescription"])   
    return df 

def _sort_genes(df, analysis, sort_by=None):
  if sort_by is None:
     return df.copy()
  if sort_by == 'median':
    sort_idxs = np.argsort(df.apply(np.nanmedian, axis=1)).values
  elif sort_by == 'mean':
    sort_idxs = np.argsort(df.apply(np.nanmean, axis=1)).values
  elif sort_by == 'agap':
    sort_idxs = np.argsort(df.reset_index()['GeneID'].values)
  elif sort_by == 'position':
    assert analysis != 'fun', "funestus cannot be sorted by position yet"
    
    gff = load_gff()
    gene_ids = gff.reset_index()['GeneID'].to_list()
    ordered_genes = gff.query(f"GeneID in {gene_ids}")['GeneID'].to_list()
    sort_idxs = [np.where(df.reset_index()['GeneID'] == gene)[0][0] for gene in ordered_genes if gene in df.reset_index()['GeneID'].to_list()]

  return df.iloc[sort_idxs, :].copy()


def query_fc_count_data(fc_data, count_data, comparison_metadata, sample_metadata, query):
    mask = comparison_metadata.eval(query).to_list()
    comparison_metadata = comparison_metadata[mask] 
    fc_data = fc_data.loc[:, mask]
    
    resistant_strains = comparison_metadata['resistant'].to_list()
    sample_mask = sample_metadata.eval("condition in @resistant_strains").to_list()
    sample_metadata = sample_metadata[sample_mask]
    count_data = count_data.loc[:, sample_mask]

    return fc_data, count_data, comparison_metadata, sample_metadata


def load_annotations():
    """
    Load pfam or go annotations for Anopheles gambiae 
    """
    pfam_df = pd.read_csv("https://github.com/sanjaynagi/AnoExpress/blob/main/resources/Anogam_long.pep_Pfamscan.seqs.gz?raw=true", sep="\s+", header=None, compression='gzip')
    go_df = pd.read_csv("https://github.com/sanjaynagi/AnoExpress/blob/main/resources/Anogam_long.pep_eggnog_diamond.emapper.annotations.GO.gz?raw=true", sep="\t", header=None, compression='gzip')
    pfam_df.columns = ["transcript", "pstart", "pend", "pfamid", "domain", "domseq"]
    go_df.columns = ['transcript', 'GO_terms']

    gene_annot_df = pfam_df.merge(go_df)
    gene_annot_df.loc[:, 'gene_id'] = gene_annot_df.loc[:, 'transcript'].str.replace("Anogam_", "").str.replace("-R[A-Z]", "", regex=True)
    return(gene_annot_df)

def filter_nas(df, fraction_na_allowed):
    """
    Filter genes with more than fraction_na_allowed of missing values
    """
    n_cols = df.shape[1]
    na_mask = df.apply(lambda x: x.isna().sum() / n_cols > fraction_na_allowed, axis=1)
    print(f"Removing {na_mask.sum()} genes with higher proportion of NAs than the threshold ({fraction_na_allowed})")
    return df.loc[~na_mask, :]




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

