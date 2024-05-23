import pandas as pd
import numpy as np
from tqdm.notebook import tqdm

taxon_query_dict = { 'gamb_colu':"species in ['gambiae','coluzzii']", 
                  'gamb_colu_arab':"species in ['gambiae','coluzzii','arabiensis']", 
                  'gamb_colu_arab_fun':"species in ['gambiae','coluzzii','arabiensis','funestus']", 
                  'fun':"species in ['funestus']"
                    }

index_col = {'fcs':'comparison',
            'pvals': 'comparison',
            'log2counts': 'sampleID'}


gff_url =  'https://vectorbase.org/common/downloads/release-68/AgambiaePEST/gff/data/VectorBase-68_AgambiaePEST.gff'


# AnoExpress
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


def load_gff(type='protein_coding_gene', query=None):

    df = pd.concat([chunk for chunk in tqdm(pd.read_csv(gff_url, sep="\t", comment="#", chunksize=10000), desc='Loading gff data from VectorBase')])
    df.columns = ['contig', 'source', 'type', 'start', 'end', 'na', 'strand', 'na2', 'attributes']
    df = df.assign(contig=lambda x: x.contig.str.split("_").str.get(1))
    
    if type:
     df = df.query(f"type == '{type}'")

    # may only work for protein_coding_genes 
    df = df.assign(GeneID=df.attributes.str.split(";", expand=True).iloc[:, 0].str.split("=").str.get(1))

    # combine 2R and 2L, 3R and 3L
    offset_2R = 61545105
    offset_3R = 53200684

    gffs = []
    for contig in tqdm(['2R', '2L', '3R', '3L']):
        df_contig = df.query("contig == @contig").copy()
        if contig == '2L':
            df_contig = df_contig.assign(contig='2RL', start=lambda x: x.start + offset_2R, end=lambda x: x.end + offset_2R)
        if contig == '3L':
            df_contig = df_contig.assign(contig='3RL', start=lambda x: x.start + offset_3R, end=lambda x: x.end + offset_3R)
        elif contig in ['3R', '2R']:
            df_contig = df_contig.assign(contig=lambda x: x.contig + 'L')
        gffs.append(df_contig)

    gff = pd.concat(gffs)
    gff = pd.concat([gff, df]).sort_values(['contig', 'start', 'end'])

    if query:
        gff = gff.query(query)
        
    return gff

def resolve_gene_id(gene_id, analysis):
    
    if isinstance(gene_id, str):
      if gene_id.startswith(('2L', '2R', '3L', '3R', 'X', '2RL', '3RL')):
        if analysis == 'fun':
          assert "Unfortunately the genome feature file does not contain AFUN identifiers, so we cannot subset by genomic span for An. funestus."
        else:

          contig, start_end = gene_id.split(':')
          start, end = start_end.replace(",", "").split('-')
          start, end = int(start), int(end)

          gff = load_gff(query=f"contig == '{contig}' and start <= {end} and end >= {start}")
          gene_id = gff.GeneID.to_list()

      elif gene_id.endswith(('.tsv', '.txt')):
          gene_id = pd.read_csv(gene_id, sep="\t", header=None).iloc[:, 0].to_list()
      elif gene_id.endswith('.csv'):
          gene_id = pd.read_csv(gene_id, header=None).iloc[:, 0].to_list()
      elif gene_id.endswith('.xlsx'):
          gene_id = pd.read_excel(gene_id, header=None).iloc[:, 0].to_list()
      
    return gene_id

def filter_low_counts(data_df, data_type, analysis, gene_id, count_threshold=5, func=np.nanmedian):
    if data_type != 'log2counts':
        count_data = data(data_type='log2counts', analysis=analysis, gene_id=gene_id)
        mask = 2**count_data.apply(func=func, axis=1) > count_threshold
    else:
        mask = 2**data_df.apply(func=func, axis=1) > count_threshold
        
    print(f"Removing {(mask == 0).sum()} genes with median counts below the threshold ({count_threshold})")
    mask = mask[mask].index.to_list()
    
    return data_df.query("GeneID in @mask")

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


def plot_gene_expression(gene_id, analysis="gamb_colu_arab_fun", microarray=False, sample_query=None, title=None, plot_type='strip', sort_by='agap', pvalue_filter=None, width=1600, height=None, save_html=None):
    """Plot fold changes of provided AGAP gene IDs from RNA-Seq 
    meta-analysis dataset

    Parameters
    ----------

    gene_id : str or list
      An AGAP identifier or list of AGAP identifiers, or AFUN if the analysis == 'fun'. Can also be a path to a file 
      containing a list of gene ids in the first column. Input file can be .tsv, .txt, or .csv, or .xlsx.
    analysis: {"gamb_colu", "gamb_colu_arab", "gamb_colu_arab_fun", "fun"}
      which analysis to load gene expression data for. analyses with more species will have less genes
      present, due to the process of finding orthologs.
    microarray: bool, optional
      whether to include the IR-Tex microarray data in the plot
    sample_query: str, optional
      A string containing a pandas query to subset the samples of interest from the comparison metadata file. For example, to plot only the
      samples from Burkina Faso, use "country == 'Burkina Faso'". Defaults to None.
    title : str
      Plot title
    plot_type : {"strip", "boxplot"}, optional
      valid options are 'strip' or 'boxplot' 
    sort_by : {"median", "mean", "agap", None}, optional
      sort by median/mean of fold changes (descending), or by AGAP, or dont sort input gene ids. 
      identifier
    pvalue_filter: float, optional
      if provided, fold-change entries with an adjusted p-value below the threshold will be removed from the plot. Default is None.
    width : int
      Width in pixels of the plotly figure
    height: int, optional
      Height in pixels of the plotly figure. Defaults to automatic sizing
    save_html : str, optional
      Path to save the plotly figure as an html file
    """
    import plotly.express as px
    import plotly.subplots as sp
    
    # load the metadata and subset to the species of interest
    df_metadata = metadata(analysis=analysis, microarray=microarray)
    df_samples = sample_metadata(analysis=analysis)

    # load fold change data, make long format and merge with metadata for hovertext
    fc_data = data(data_type="fcs", analysis=analysis, microarray=microarray, sample_query=sample_query, gene_id=gene_id, sort_by=sort_by, annotations=True, pvalue_filter=pvalue_filter).reset_index()
    # load count data, make long format and merge with metadata for hovertext
    count_data = data(data_type="log2counts", analysis=analysis, microarray=microarray, gene_id=gene_id, sample_query=sample_query, sort_by=None)
    count_data = count_data.loc[fc_data['GeneID']].reset_index()

    if sample_query:
      fc_data, count_data, df_metadata, df_samples = query_fc_count_data(fc_data=fc_data, count_data=count_data, df_metadata=df_metadata, sample_metadata=df_samples, query=sample_query)

    count_data = count_data.melt(id_vars='GeneID', var_name='sampleID', value_name='log2_counts')
    count_data = count_data.merge(df_samples, how='left').assign(counts = lambda x: np.round(2**x.log2_counts, 0))

    fc_data.loc[:, 'Label'] = [id_ + " | " + name if name != "" else id_ for id_, name in zip(fc_data['GeneID'].fillna(""), fc_data['GeneName'].fillna(""))]
    fc_data = fc_data.drop(columns=['GeneName', 'GeneID', 'GeneDescription']).melt(id_vars='Label', var_name='comparison', value_name='log2FC')
    fc_data = fc_data.merge(df_metadata, how='left')

    #reorder fc data to match count_data ordering (important for legends/colours matching)
    fc_data = fc_data.set_index('species').loc[count_data.species.unique()].reset_index()

    if not height:
      height = np.min([fc_data.shape[0]*12, 2500])
    
    if not title:
      title = ""

    myplot = px.box if plot_type == 'boxplot' else px.strip
    figure1 = myplot(
          fc_data, 
          y='Label', 
          x='log2FC', 
          color='species',
          title="title", 
          hover_data=['resistant', 'susceptible', 'species', 'country', 'technology'],
          template='simple_white'
        )

    figure2 = myplot(
        count_data, 
        x='counts', 
        y='GeneID', 
        color='species', 
        orientation='h', 
        hover_data=['sampleID', 'species', 'country'],
        template='simple_white',
        )

    for i in range(len(figure2.data)):
      figure2.data[i]['showlegend'] = False
      
    figure1_traces = []
    figure2_traces = []
    for trace in range(len(figure1["data"])):
        figure1_traces.append(figure1["data"][trace])
    for trace in range(len(figure2["data"])):
        figure2_traces.append(figure2["data"][trace])

    # Create a 1x2 subplot
    final_figure = sp.make_subplots(rows=1, cols=2, column_widths=[0.7, 0.3]) 
    # Get the Express fig broken down as traces and add the traces to the proper plot within in the subplot
    for traces in figure1_traces:
        final_figure.append_trace(traces, row=1, col=1)
    for traces in figure2_traces:
        final_figure.append_trace(traces, row=1, col=2)

    # reset boxmode to group so species are separate, resize plots, set axes titles, labels and vlines
    final_figure.layout['boxmode'] = 'group'
    final_figure.layout['xaxis2']['domain'] = (0.65, 1.0)
    final_figure.update_layout(title_text=title, title_x=0.5, width=width, height=height, template='simple_white')
    final_figure.update_yaxes(title_text="Gene", row=1, col=1, title_font = {"size": 18}, tickfont={"size":14}),
    final_figure.update_yaxes(showticklabels=False, row=1, col=2)
    final_figure.update_xaxes(title_text="log2 fold change", row=1, col=1, title_font = {"size": 18})
    final_figure.update_xaxes(title_text="counts", row=1, col=2, title_font = {"size": 18})
    for i in [1,2]: final_figure.add_vline(0,  line_width=1, line_dash="dash", line_color="grey", row=1, col=i)
    
    if save_html:
       final_figure.write_html(save_html)

    return(final_figure)


def _gene_ids_from_annotation(gene_annot_df, annotation):
    """
    Extract gene ids from gene_annot_df based on annotation
    """      
    if isinstance(annotation, str):
      annotation = [annotation]

    gene_list = np.array([]) 
    for annot in annotation:
      if annot.startswith("GO"):
        ids = gene_annot_df.query(f"GO_terms.str.contains('{annot}', na=False)", engine='python')['gene_id'].to_numpy()
      else:
        ids = gene_annot_df.query("domain == @annot")['gene_id'].to_numpy()
      gene_list = np.hstack([gene_list, ids])
    
    return np.unique(gene_list)
          
   
    


def plot_gene_family_expression(gene_identifier, analysis, title, microarray=False, plot_type='strip', sort_by='median', width=1600, height=None):
  """Plot gene expression of gene families belonging to GO terms or PFAM domains

  Parameters
  ----------

  gene_identifier : str or list
    An AGAP identifier or list of AGAP identifiers
  analysis: {"gamb_colu", "gamb_colu_arab", "gamb_colu_arab_fun"}
    which analysis to load gene expression data for. analyses with more species will have less genes
    present, due to the process of finding orthologs.
  title : str
    Plot title
  microarray: bool, optional,
    if True, load microarray data from Ir-tex. If False, load RNAseq data only. Defaults to False
  plot_type : {"strip", "boxplot"}, optional
    valid options are 'strip' or 'boxplot' 
  sort_by : {"median", "mean", "agap"}, optional
    sort by median/mean of fold changes (descending), or by AGAP
    identifier
  width : int
    Width in pixels of the plotly figure
  height: int, optional
    Height in pixels of the plotly figure. Defaults to automatic sizing
  """
  assert analysis != 'fun', "GO terms and PFAM domains cannot be searched against An. funestus. Please use a different analysis."
  
  # Read in .csv file containing pfam and go terms
  gene_annot_df = load_annotations()
  gene_ids = _gene_ids_from_annotation(gene_annot_df, gene_identifier)
  fig = plot_gene_expression(gene_id=gene_ids, microarray=microarray, title=title, analysis=analysis, plot_type=plot_type, sort_by=sort_by, width=width, height=height)

  return(fig)

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


def load_genes_for_enrichment(analysis, func, gene_ids, percentile, microarray, low_count_filter=None):
   
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
      top_geneIDs = resolve_gene_id(gene_id=gene_ids, analysis=analysis)

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



def plot_heatmap(analysis, gene_id=None, query_annotation=None, query_func=np.nanmedian, query_fc=None, query_name='median', cmap=None, cbar_pos=None, figsize=None):
    """
    Plot a heatmap of the top 100 genes ranked by user input function.
    
    PARAMETERS
    ----------
    analysis: {"gamb_colu", "gamb_colu_arab", "gamb_colu_arab_fun", "fun"}
      which analysis to load gene expression data for. analyses with more species will have less genes
      present, due to the process of finding orthologs.
    query_annotation: str, optional
      Either a GO term or a PFAM domain to select genes to plot heatmap of. Defaults to None.
    query_func: function, optional
      function to rank genes by (such as np.nanmedian, np.nanmean). Defaults to np.nanmedian
    query_fc: float, optional
      fold change threshold to select genes to plot heatmap of. Defaults to None.
    query_name: str, optional
      name of the function to rank genes by. Defaults to 'median'.
    cmap: str, optional
      colormap to use for heatmap.
    cbar_pos: {"left", "right", "top", "bottom"}, optional
      position of colorbar. Defaults to None.
    figsize: tuple, optional
      size of figure.
    """
    import seaborn as sns
    # load metadata

    if gene_id:      
      fc_data = data(data_type='fcs', gene_id=gene_id, analysis=analysis, microarray=False, annotations=True, sort_by=None).reset_index()
    elif not gene_id:
      fc_data = data(data_type='fcs', analysis=analysis, microarray=False, annotations=True, sort_by=None).reset_index()
      fc_ranked = load_candidates(analysis=analysis, name=query_name, func=query_func, query_annotation=query_annotation, query_fc=query_fc)
      fc_genes = fc_ranked.loc[:, 'GeneID'].to_list()
      fc_data = fc_data.query(f"GeneID in {fc_genes}").copy()
    
    fc_data.loc[:, 'Label'] = [id_ + " | " + name if name != "" else id_ for id_, name in zip(fc_data['GeneID'].fillna(""), fc_data['GeneName'].fillna(""))]
    fc_data = fc_data.set_index("Label").drop(columns=['GeneName', 'GeneID', 'GeneDescription'])
    fc_data.columns = [c.replace("_log2FoldChange", "").replace("_", " ") for c in fc_data.columns]
    mask = fc_data.isnull()

    if fc_data.empty or fc_data.shape[0] == 1:
      print(f"Too few observations for gene selection or {query_annotation} and FC of greater than {query_fc}")
      return
    
    if not figsize:
        height = np.max([fc_data.shape[0]/2.5, 4])
        figsize = [10, height]

    cg = sns.clustermap(
        fc_data.fillna(0), 
        mask=mask, 
        cbar_pos=cbar_pos, 
        cmap=cmap,
        figsize=figsize, 
        tree_kws={'linewidths':3, 'colors':'darkgrey'},
        linewidths=2,
        yticklabels=True,
    )

    cg.ax_col_dendrogram.set_visible(False)




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




def plot_contig_expression_track(
    contig,
    analysis="gamb_colu_arab_fun",
    data_type='fcs',
    microarray=False,
    pvalue_filter=None,
    size=10,
    step=5,
    title=None, 
    width=800, 
    height=600, 
    palette=None,
    sizing_mode='stretch_width',
    x_range=None,
    y_range=None,
    show=False
):
    """
    Plot fold change data for a given contig. Plots both raw fold change data for each experiment and a moving average

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
    title: str, optional
      plot title. Defaults to None
    width: int, optional
      width of plot. Defaults to 800
    height: int, optional
      height of plot. Defaults to 600
    color: bool, optional
      Colour points by taxon of the experiment
    sizing_mode: {"stretch_width", "scale_width", "scale_both", "scale_height", "stretch_both"}, optional
      sizing mode for plot. Defaults to 'stretch_width'
    x_range: tuple, optional
      x axis range. Defaults to None
    y_range: tuple, optional
      y axis range. Defaults to None
    show: bool, optional
      whether to show the plot. Defaults to False
    """
    import bokeh
    import bokeh.plotting as bkplt

    fold_change_df, windowed_fold_change_df = contig_expression(
       contig=contig, 
       analysis=analysis, 
       data_type=data_type, 
       microarray=microarray,
       pvalue_filter=pvalue_filter, 
       size=size, 
       step=step
       )
    
    df_metadata = metadata(analysis=analysis, microarray=microarray)
    fold_change_df = fold_change_df.merge(df_metadata)
    
    if palette:
        # add color column 
        fold_change_df['color'] = np.select(
            [
                fold_change_df['species'] == "coluzzii", 
                fold_change_df['species'] == "gambiae",
                fold_change_df['species'] == "arabiensis",
                fold_change_df['species'] == "funestus"
            ], 
              list(palette), 
            
            default='grey'
            )
        color = 'color'
    else:
        color = 'grey'
    
    # determine X axis range
    x_min = fold_change_df.midpoint.to_numpy()[0]
    x_max = fold_change_df.midpoint.to_numpy()[-1]
    if x_range is None:
        x_range = bokeh.models.Range1d(x_min, x_max, bounds="auto")

    # create a figure
    xwheel_zoom = bokeh.models.WheelZoomTool(
        dimensions="width", maintain_focus=False
    )
    
    tooltips = [
            ("Gene ID", "@GeneID"),
            ("Gene Name","@GeneName"),
            ("Gene Description", "@GeneDescription"),
            ("Experiment", "@comparison"),
            ("Taxon", "@species"),
            ("Technology", "@technology"),
            ("Country", "@country"),
            ("Log2 FC", "@fold_change"),
        ]

    fig = bkplt.figure(
        title=title,
        tools=["xpan", "xzoom_in", "xzoom_out", xwheel_zoom, "reset", "hover"],
        active_scroll=xwheel_zoom,
        active_drag="xpan",
        sizing_mode=sizing_mode,
        width=width,
        height=height,
        toolbar_location="above",
        tooltips=tooltips,
        x_range=x_range,
        y_range=y_range,
        output_backend='webgl'
    )

    # plot 
    fig.circle(
        x="midpoint",
        y="fold_change",
        size=4,
        line_width=1,
        line_color=color,
        fill_color=None,
        source=fold_change_df,
    )
    
    fig.line(
        x="midpoint",
        y="median_fc",
        line_width=2,
        line_color="black",
        source=windowed_fold_change_df,
    )

    # tidy up the plot
    fig.yaxis.axis_label = "Log2 Fold Change"

    if show:
        bkplt.show(fig)
    return fig 






def plot_contig_expression(contig, analysis, data_type='fcs', microarray=False, size=10, step=5, pvalue_filter=None, palette=None,y_range=(-10,15), height=400, width=600, title=None, show=False):
    """
    Plot fold change data for a given contig with a gene track. 

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
    title: str, optional
      plot title. Defaults to None
    width: int, optional
      width of plot. Defaults to 800
    height: int, optional
      height of plot. Defaults to 600
    y_range: tuple, optional
      y axis range. Defaults to None
    show: bool, optional
      whether to show the plot. Defaults to False
    """

    import bokeh
    import bokeh.plotting as bkplt
    import malariagen_data
    ag3 = malariagen_data.Ag3()

    fig1 = plot_contig_expression_track(
                                    contig=contig,
                                    analysis=analysis,
                                    data_type=data_type,
                                    microarray=microarray,
                                    pvalue_filter=pvalue_filter,
                                    palette=palette,
                                    size=size,
                                    step=step,
                                    y_range=y_range,
                                    height=height,
                                    width=width,
                                    title=title,
                                    show=False,
                                    )
    fig1.xaxis.visible = False
    # plot genes
    fig2 = ag3.plot_genes(
        region=contig,
        sizing_mode="stretch_width",
        width=width,
        height=100,
        x_range=fig1.x_range,
        show=False,
    )
    # combine plots into a single figure
    fig = bokeh.layouts.gridplot(
        [fig1, fig2],
        ncols=1,
        toolbar_location="above",
        merge_tools=True,
        sizing_mode="stretch_width",
    )
    if show:
      bkplt.show(fig)
    return fig


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