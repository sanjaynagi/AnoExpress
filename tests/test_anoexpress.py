import anoexpress as xpress
import pytest
import pandas as pd
import numpy as np
from pandas.testing import assert_frame_equal
from numpy.testing import assert_allclose

gene= 'AGAP006227'
gene_ids = ['AGAP006222', 'AGAP006226', 'AGAP006227']

@pytest.mark.parametrize(    
        "data_type",
    ["fcs", "pvals", "log2counts"],
    )
@pytest.mark.parametrize(
    "analysis",
    ['gamb_colu', 'gamb_colu_arab', 'gamb_colu_arab_fun', 'fun']
)
def test_load_results_arrays(data_type, analysis):
    df = xpress.load_results_arrays(data_type=data_type, analysis=analysis)
    assert isinstance(df, pd.DataFrame)
    assert not df.empty





@pytest.mark.parametrize(    
        "data_type",
    ["fcs", "pvals"],
    )
def test_load_irtex_arrays(data_type):
    df = xpress.load_results_arrays(data_type=data_type, analysis="irtex")
    assert isinstance(df, pd.DataFrame)
    assert not df.empty

@pytest.mark.parametrize("microarray", [True, False])
def test_data_irtex(microarray):

    data_df = xpress.data(
        data_type="fcs",
        analysis="gamb_colu_arab", 
        gene_id=None, 
        microarray=microarray
    )
    assert data_df is not None
    assert not data_df.empty
    assert isinstance(data_df, pd.DataFrame)


@pytest.mark.parametrize("analysis", ['gamb_colu', 'gamb_colu_arab', 'gamb_colu_arab_fun'])
def test_data_analyses(analysis):

    data_df = xpress.data(
        data_type="pvals",
        analysis=analysis, 
        gene_id=None, 
        microarray=False
    )
    assert data_df is not None
    assert not data_df.empty
    assert isinstance(data_df, pd.DataFrame)

@pytest.mark.parametrize('data_type', ["fcs", "pvals", "log2counts"])
@pytest.mark.parametrize("gene_id",    [None, gene, gene_ids, "2RL:28,500,500-28,530,000"])
def test_data_types_genes(data_type, gene_id):

    data_df = xpress.data(
        data_type=data_type,
        analysis="gamb_colu", 
        gene_id=gene_id, 
        microarray=True
    )
    assert data_df is not None
    assert not data_df.empty
    assert isinstance(data_df, pd.DataFrame)


@pytest.mark.parametrize("gene_id",    [None, gene, gene_ids])
def test_data_sorting(gene_id):

    data_df = xpress.data(
        data_type="fcs",
        analysis="gamb_colu", 
        gene_id= gene_id, 
        microarray=False,
        sort_by="position"
    )
    assert data_df is not None
    assert not data_df.empty
    assert isinstance(data_df, pd.DataFrame)

def test_data_pvalue_filter():

    data_df = xpress.data(
        data_type="fcs",
        analysis="gamb_colu", 
        microarray=False,
        sort_by=None,
        pvalue_filter=0.05
    )
    assert data_df is not None
    assert not data_df.empty
    assert isinstance(data_df, pd.DataFrame)

def test_data_low_count_filter(low_count_filter=10):

    data_df = xpress.data(
        data_type="fcs",
        analysis="gamb_colu", 
        microarray=False,
        sort_by=None,
    )

    data_low_df = xpress.data(
        data_type="fcs",
        analysis="gamb_colu", 
        microarray=False,
        low_count_filter=low_count_filter,
        sort_by=None)

    assert data_low_df is not None
    assert not data_low_df.empty
    assert isinstance(data_low_df, pd.DataFrame)
    assert data_df.shape[0] > data_low_df.shape[0]

def test_data_sample_query(query="country == 'Burkina Faso'"):

    data_df = xpress.data(
        data_type="fcs",
        analysis="gamb_colu", 
        microarray=False,
        sort_by=None,
    )

    data_low_df = xpress.data(
        data_type="fcs",
        analysis="gamb_colu", 
        microarray=False,
        sample_query=query,
        sort_by=None)

    assert data_low_df is not None
    assert not data_low_df.empty
    assert isinstance(data_low_df, pd.DataFrame)
    assert data_df.shape[1] > data_low_df.shape[1]

@pytest.mark.parametrize(
    "gene_ids",
    [gene, gene_ids]
)
@pytest.mark.parametrize(
    "microarray",
    [True, False]
)
def test_plot_gene_expression(gene_ids, microarray):

    xpress.plot_gene_expression(
        gene_id=gene_ids, 
        analysis="gamb_colu_arab_fun", 
        microarray=microarray, 
        sort_by="agap"
        )
    

@pytest.mark.parametrize(
    "plot_type", 
    ['strip', 'boxplot']
)
def test_plot_gene_expression_type(plot_type):

    xpress.plot_gene_expression(
        gene_id="AGAP006227", 
        analysis="gamb_colu", 
        microarray=False, 
        plot_type=plot_type,
        )


@pytest.mark.parametrize(
    "gene_id", 
    ['2RL:28,480,500-28,500,000', 'X:8,500,500-8,530,000']
)
def test_plot_gene_expression_spans(gene_id):

    xpress.plot_gene_expression(
        gene_id=gene_id, 
        analysis="gamb_colu", 
        microarray=False, 
        plot_type='strip',
        )


def test_load_candidates():
    
    df = xpress.load_candidates(analysis="gamb_colu_arab_fun", name="median", func=np.nanmedian, fraction_na_allowed=0.5, low_count_filter=5)

    assert isinstance(df, pd.DataFrame)
    assert not df.empty




### hypergeometric functions ###
def test_go_hypergeometric():
    go = xpress.go_hypergeometric(analysis='gamb_colu_arab_fun', func=np.nanmedian
        )
    assert isinstance(go, pd.DataFrame)
    assert go.iloc[0,0] == 'GO:0042302' # check first value is correct

def test_pfam_hypergeometric():
    pfam = xpress.pfam_hypergeometric(analysis='gamb_colu_arab_fun', func=np.nanmedian
        )
    assert isinstance(pfam, pd.DataFrame)
    assert pfam.iloc[0,0] == 'C_tripleX' # check first value is correct


def test_pfam_hypergeometric_gene_ids():
    pfam = xpress.pfam_hypergeometric(analysis='gamb_colu_arab_fun', gene_ids=gene_ids
        )
    assert isinstance(pfam, pd.DataFrame)


def test_kegg_hypergeometric():
    kegg = xpress.kegg_hypergeometric(analysis='gamb_colu_arab_fun', func=np.nanmedian
        )
    assert isinstance(kegg, pd.DataFrame)
    assert kegg.iloc[0,0] == 'aga00982' # check first value is correct



