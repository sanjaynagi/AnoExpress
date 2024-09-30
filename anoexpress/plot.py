import numpy as np
import pandas as pd
import plotly.express as px
import plotly.subplots as sp
import seaborn as sns

from .data import data, metadata, sample_metadata, query_fc_count_data, load_annotations
from .utils import _gene_ids_from_annotation
from .candidates import load_candidates, contig_expression




def plot_gene_expression(gene_id, analysis="gamb_colu_arab_fun", microarray=False, sample_query=None, title=None, plot_type='strip', sort_by='agap', gff_method='malariagen_data', pvalue_filter=None, width=1600, height=None, save_html=None):
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
    sort_by : {"median", "mean", "agap", "position", None}, optional
      sort by median/mean of fold changes (descending), or by AGAP, or dont sort input gene ids. 
      identifier
    gff_method : {"malariagen_data", "vectorbase"}, optional
      method to use to load gff, for sorting genes by position. Defaults to 'malariagen_data
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
    fc_data = data(data_type="fcs", analysis=analysis, microarray=microarray, sample_query=sample_query, gene_id=gene_id, sort_by=sort_by, annotations=True, pvalue_filter=pvalue_filter, gff_method=gff_method).reset_index()
    # load count data, make long format and merge with metadata for hovertext
    count_data = data(data_type="log2counts", analysis=analysis, microarray=microarray, gene_id=gene_id, sample_query=sample_query, sort_by=None, gff_method=gff_method)
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







def plot_gene_family_expression(gene_identifier, analysis, title, microarray=False, plot_type='strip', sort_by='median', gff_method="malariagen_data", width=1600, height=None):
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
  fig = plot_gene_expression(gene_id=gene_ids, microarray=microarray, title=title, analysis=analysis, plot_type=plot_type, sort_by=sort_by, gff_method=gff_method, width=width, height=height)

  return(fig)




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
