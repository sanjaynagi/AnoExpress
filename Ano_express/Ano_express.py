import pandas as pd
import numpy as np

# Ano-express
def load_results_arrays(data_type, analysis):
    """
    Load the counts data for a given analysis and sample query
    """
    assert data_type in ['log2counts', 'fcs', 'pvals'], "data_type must be either 'log2counts', 'fcs' or 'pvals'"
    assert analysis in ['gamb_colu', 'gamb_colu_arab', 'gamb_colu_arab_fun', 'fun', 'irtex'], "analysis must be either 'gamb_colu', 'gamb_colu_arab', 'gamb_colu_arab_fun', 'fun' or 'irtex'"
    df_data = pd.read_csv(f"https://raw.githubusercontent.com/sanjaynagi/Ano-express/main/results/{data_type}.{analysis}.tsv", sep="\t")
    df_data = df_data.drop(columns=['GeneDescription']).set_index(['GeneID', 'GeneName']) if data_type == 'fcs' and analysis != 'irtex' else df_data.set_index("GeneID")
    return(df_data)

def xpress_metadata():
    """
    Load the sample and comparisons metadata in a pandas dataframe
    """
    sample_metadata = pd.read_csv("https://raw.githubusercontent.com/sanjaynagi/ano-express/main/config/sample_metadata.tsv", sep="\t").rename(columns={'colData':'sampleID'})
    comparison_metadata = pd.read_csv("https://raw.githubusercontent.com/sanjaynagi/ano-express/main/config/comparison_metadata.tsv", sep="\t")
    return(comparison_metadata, sample_metadata)

def irtex_metadata():
    """
    Load the ir-tex fold change data for a given analysis and sample query
    """   
    comparison_metadata = pd.read_csv("https://github.com/sanjaynagi/Ano-express/blob/main/config/irtex_metadata.tsv?raw=true", sep="\t")
    return(comparison_metadata)

def metadata(analysis, microarray=True, sample_metadata=False):
    """
    Load the comparisons metadata from both Ano-express and IR-Tex in a pandas dataframe
    """
    # load metadata from Ano-express
    metadata, xpress_samples = xpress_metadata()   
    metadata = metadata.assign(technology='rnaseq')

    # load metadata from IR-Tex
    if microarray == True:
        irtex_meta = irtex_metadata()
        irtex_meta = irtex_meta.assign(technology='microarray')
        metadata = pd.concat([metadata, irtex_meta]).reset_index(drop=True)

    # subset to the species of interest
    comparisons_df = species_query(metadata, analysis)
    
    if sample_metadata == True:
        return(comparisons_df, xpress_samples)
    else:
      return(comparisons_df)


def data(data_type, analysis, microarray=True):
    """
    Load the combined data for a given analysis and sample query
    """
    assert analysis in ['gamb_colu', 'gamb_colu_arab', 'gamb_colu_arab_fun', 'fun'], "analysis must be one of 'gamb_colu', 'gamb_colu_arab', 'gamb_colu_arab_fun', 'fun'"
    # load the metadata and subset to the species of interest
    comparisons_df = metadata(analysis=analysis, microarray=microarray)
    comparisons_df = species_query(comparisons_df=comparisons_df, analysis=analysis)
    comparisons_ids = comparisons_df.loc[:, 'comparison'].to_list()

    # load the data and merge
    df = load_results_arrays(data_type=data_type, analysis=analysis)
    
    # load ir-tex data and merge
    if microarray == True:
        irtex_df = load_results_arrays(data_type=data_type, analysis='irtex')
        df = df.reset_index().merge(irtex_df, on='GeneID', how='left').set_index(['GeneID', 'GeneName'])

    # subset to the species comparisons of interest
    df = df.loc[:, comparisons_ids]
    return(df)

def species_query(comparisons_df, analysis):
    """
    Subset the comparisons metadata to the species of interest
    """
    if analysis == 'gamb_colu':
        comparisons_df = comparisons_df.query('species in ["gambiae","coluzzii"]')
    elif analysis == 'gamb_colu_arab':
        comparisons_df = comparisons_df.query('species in ["gambiae","coluzzii","arabiensis"]')
    elif analysis == 'gamb_colu_arab_fun':
        comparisons_df = comparisons_df.query('species in ["gambiae","coluzzii","arabiensis","funestus"]')
    elif analysis == 'fun':
        comparisons_df = comparisons_df.query('species in ["funestus"]')
    return(comparisons_df)


def plot_gene_expression(gene_id, analysis="gamb_colu_arab_fun", microarray=False, title=None, plot_type='strip', sort_by='agap', width=1600, height=None):
    """Plot fold changes of provided AGAP gene IDs from RNA-Seq 
    meta-analysis dataset

    Parameters
    ----------
    gene_id : str or list
      An AGAP identifier or list of AGAP identifiers, or AFUN if the analysis == 'fun'.
    analysis: {"gamb_colu", "gamb_colu_arab", "gamb_colu_arab_fun", "fun"}
      which analysis to load gene expression data for. analyses with more species will have less genes
      present, due to the process of finding orthologs.
    microarray: bool, optional
      whether to include the IR-Tex microarray data in the plot
    title : str
      Plot title
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
    import plotly.express as px
    import plotly.subplots as sp
      
    fc_data = data("fcs", analysis=analysis, microarray=microarray).reset_index()
    count_data = pd.read_csv(f"https://raw.githubusercontent.com/sanjaynagi/ano-expressir/main/results/log2counts.{analysis}.tsv", sep="\t")
    comp_metadata, sample_metadata = metadata(analysis=analysis, microarray=microarray, sample_metadata=True)

    fam_fc_data = fc_data.query("GeneID in @gene_id").copy()
    fam_count_data = count_data.query("GeneID in @gene_id").copy()

    if sort_by == 'median':
        sort_idxs = np.argsort(fam_fc_data.set_index(['GeneID', 'GeneName']).apply(np.nanmedian, axis=1)).values
    elif sort_by == 'mean':
        sort_idxs = np.argsort(fam_fc_data.set_index(['GeneID', 'GeneName']).apply(np.nanmean, axis=1)).values
    elif sort_by == 'agap':
        sort_idxs = np.argsort(fam_fc_data['GeneID'].values)[::-1] 
        
    fam_fc_data = fam_fc_data.iloc[sort_idxs, :].copy()
    fam_count_data = fam_count_data.set_index("GeneID").loc[fam_fc_data['GeneID'].to_list(), :].reset_index().copy()

    fam_fc_data.loc[:, 'Label'] = [id_ + " | " + name if name != "" else id_ for id_, name in zip(fam_fc_data['GeneID'].fillna(""), fam_fc_data['GeneName'].fillna(""))]
    fam_fc_data = fam_fc_data.drop(columns=['GeneName', 'GeneID']).melt(id_vars='Label', var_name='comparison', value_name='log2FC')
    fam_count_data = fam_count_data.melt(id_vars='GeneID', var_name='sampleID', value_name='log2_counts')
    fam_fc_data.loc[:, 'comparison'] = fam_fc_data['comparison'].str.replace("_log2FoldChange", "")
    fam_fc_data = fam_fc_data.merge(comp_metadata, how='left')
    fam_count_data = fam_count_data.merge(sample_metadata, how='left').assign(counts = lambda x: np.round(2**x.log2_counts, 0))

    if not height:
      height = np.min([fam_fc_data.shape[0]*12, 2500])
    
    if not title:
      title = ""

    myplot = px.box if plot_type == 'boxplot' else px.strip
    figure1 = myplot(
          fam_fc_data, 
          y='Label', 
          x='log2FC', 
          color='species',
          title="title", 
          hover_data=['resistant', 'susceptible', 'species', 'country', 'technology'],
          template='ggplot2'
        )

    figure2 = myplot(
        fam_count_data, 
        x='counts', 
        y='GeneID', 
        color='species', 
        orientation='h', 
        hover_data=['sampleID', 'species', 'country'],
        template='ggplot2',
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
    final_figure.update_layout(title_text=title, title_x=0.5, width=width, height=height)
    final_figure.update_yaxes(title_text="Gene", row=1, col=1, title_font = {"size": 18}, tickfont={"size":14}),
    final_figure.update_yaxes(showticklabels=False, row=1, col=2)
    final_figure.update_xaxes(title_text="log2 fold change", row=1, col=1, title_font = {"size": 18})
    final_figure.update_xaxes(title_text="counts", row=1, col=2, title_font = {"size": 18})
    for i in [1,2]: final_figure.add_vline(0,  line_width=1, line_dash="dash", line_color="grey", row=1, col=i)
    
    return(final_figure)


def gene_ids_from_annotation(gene_annot_df, annotation):
    """
    Extract gene ids from gene_annot_df based on annotation
    """
    if isinstance(annotation, list):
        gene_list = np.array([])
        if annotation[0].startswith("GO"):
            for go in annotation:
              ids = gene_annot_df.query(f"GO_terms.str.contains('{go}', na=False)", engine='python')['gene_id'].to_numpy()
              gene_list = np.hstack([gene_list, ids])
            return(np.unique(gene_list))
        else:
          for dom in annotation:
              ids = gene_annot_df.query("domain == @annotation")['gene_id'].to_numpy()
              gene_list = np.hstack([gene_list, ids])
          return(np.unique(gene_list))
    else:
        if annotation.startswith("GO"): 
          return(gene_annot_df.query(f"GO_terms.str.contains('{annotation}', na=False)", engine='python')['gene_id'].to_numpy())
        else:
          return(gene_annot_df.query("domain == @annotation")['gene_id'].to_numpy())



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
  gene_ids = gene_ids_from_annotation(gene_annot_df, gene_identifier)
  fig = plot_gene_expression(gene_id=gene_ids, microarray=microarray, title=title, analysis=analysis, plot_type=plot_type, sort_by=sort_by, width=width, height=height)

  return(fig)

def load_annotations():
    """
    Load pfam or go annotations for Anopheles gambiae 
    """
    pfam_df = pd.read_csv("https://github.com/sanjaynagi/ano-expressir/blob/main/resources/Anogam_long.pep_Pfamscan.seqs.gz?raw=true", sep="\s+", header=None, compression='gzip')
    go_df = pd.read_csv("https://github.com/sanjaynagi/ano-expressir/blob/main/resources/Anogam_long.pep_eggnog_diamond.emapper.annotations.GO.gz?raw=true", sep="\t", header=None, compression='gzip')
    pfam_df.columns = ["transcript", "pstart", "pend", "pfamid", "domain", "domseq"]
    go_df.columns = ['transcript', 'GO_terms']

    gene_annot_df = pfam_df.merge(go_df)
    gene_annot_df.loc[:, 'gene_id'] = gene_annot_df.loc[:, 'transcript'].str.replace("Anogam_", "").str.replace("-R[A-Z]", "", regex=True)
    return(gene_annot_df)



def load_candidates(analysis, name='median', func=np.nanmedian, query_annotation=None, query_fc=None):
    """
    Load the candidate genes for a given analysis. Optionally, filter by annotation or fold change data.
    """
    fc_data = pd.read_csv(f"https://raw.githubusercontent.com/sanjaynagi/ano-expressir/main/results/fcs.{analysis}.tsv", sep="\t")
    fc_data = fc_data.set_index(['GeneID', 'GeneName', 'GeneDescription'])

    if query_annotation is not None:
      gene_annot_df = load_annotations()
      gene_ids = gene_ids_from_annotation(gene_annot_df=gene_annot_df, annotation=query_annotation)
      fc_data = fc_data.query("GeneID in @gene_ids")
      assert not fc_data.empty, "No genes were found for the selection. It is possible these genes were removed by the ortholog finding process"
    
    fc_ranked = fc_data.apply(func, axis=1).to_frame().rename(columns={0:f'{name} log2 Fold Change'}).copy()
    fc_ranked = fc_ranked.sort_values(f'{name} log2 Fold Change', ascending=False)
    fc_ranked = fc_ranked.reset_index()
    fc_ranked.loc[:, f'{name} Fold Change'] = np.round(2**fc_ranked.loc[:, f'{name} log2 Fold Change'], 2)

    if query_fc is not None:
       fc_ranked = fc_ranked.query(f'`{name} Fold Change` > {query_fc}')

    return(fc_ranked)

def go_hypergeometric(analysis, name, func, percentile=0.05):

    fc_data = pd.read_csv(f"https://raw.githubusercontent.com/sanjaynagi/ano-expressir/main/results/fcs.{analysis}.tsv", sep="\t")
    fc_genes = fc_data.reset_index()['GeneID'].to_list()

    # get top % percentile genes ranked by func
    fc_ranked = load_candidates(analysis=analysis, name=name, func=func)
    percentile_idx = fc_ranked.reset_index()['GeneID'].unique().shape[0] * percentile
    top_geneIDs = fc_ranked.reset_index().loc[:, 'GeneID'][:int(percentile_idx)] 

    # load gene annotation file 
    gaf_df = pd.read_csv("https://raw.githubusercontent.com/sanjaynagi/ano-expressir/main/resources/AgamP4.gaf", sep="\t")
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

def pfam_hypergeometric(analysis, name, func, percentile=0.05):
    """
    This function performs a hypergeometric test on pfam domains on the top percentile of genes ranked by func
    """

    # get all genes
    fc_data = pd.read_csv(f"https://raw.githubusercontent.com/sanjaynagi/ano-expressir/main/results/fcs.{analysis}.tsv", sep="\t")
    fc_genes = fc_data.reset_index()['GeneID'].to_list()

    # get top 5% percentile genes ranked by median
    fc_ranked = load_candidates(analysis=analysis, name=name, func=func)
    percentile_idx = fc_ranked.reset_index()['GeneID'].unique().shape[0] * percentile
    top_geneIDs = fc_ranked.reset_index().loc[:, 'GeneID'][:int(percentile_idx)] 

    # load gene annotation file 
    pfam_df = pd.read_csv("https://github.com/sanjaynagi/ano-expressir/blob/main/resources/Anogam_long.pep_Pfamscan.seqs.gz?raw=true", sep="\s+", header=None, compression='gzip').iloc[:, [0,4]]
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



def plot_heatmap(analysis, query_annotation=None, query_func=np.nanmedian, query_fc=None, query_name='median', cmap=None, cbar_pos=None, figsize=None):
    
    import seaborn as sns
    # load metadata
    fc_data = pd.read_csv(f"https://raw.githubusercontent.com/sanjaynagi/ano-expressir/main/results/fcs.{analysis}.tsv", sep="\t") 
    fc_ranked = load_candidates(analysis=analysis, name=query_name, func=query_func, query_annotation=query_annotation, query_fc=query_fc)
    fc_genes = fc_ranked.loc[:, 'GeneID']
    fam_data = fc_data.query("GeneID in @fc_genes").copy()
    
    fam_data.loc[:, 'Label'] = [id_ + " | " + name if name != "" else id_ for id_, name in zip(fam_data['GeneID'].fillna(""), fam_data['GeneName'].fillna(""))]
    fam_data = fam_data.set_index("Label").drop(columns=['GeneName', 'GeneID', 'GeneDescription'])
    fam_data.columns = [c.replace("_log2FoldChange", "").replace("_", " ") for c in fam_data.columns]
    mask = fam_data.isnull()

    if fam_data.empty or fam_data.shape[0] == 1:
      print(f"Too few observations for {query_annotation} and FC of greater than {query_fc}")
      return
    
    if not figsize:
        height = np.max([fam_data.shape[0]/2.5, 4])
        figsize = [10, height]

    cg = sns.clustermap(
        fam_data.fillna(0), 
        mask=mask, 
        cbar_pos=cbar_pos, 
        cmap=cmap,
        figsize=figsize, 
        tree_kws={'linewidths':3, 'colors':'darkgrey'},
        linewidths=2,
        yticklabels=True,
    )

    cg.ax_col_dendrogram.set_visible(False)