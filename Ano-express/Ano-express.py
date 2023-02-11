import pandas as pd
import numpy as np

# Ano-express


def load_results_arrays(data_type, analysis, sample_query=None):
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

    #Create a 1x2 subplot
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


