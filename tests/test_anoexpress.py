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
@pytest.mark.parametrize("gene_id",    [None, gene, gene_ids])
def test_data_types_genes(data_type, gene_id):

    data_df = xpress.data(
        data_type=data_type,
        analysis="gamb_colu", 
        gene_id= gene_id, 
        microarray=True
    )
    assert data_df is not None
    assert not data_df.empty
    assert isinstance(data_df, pd.DataFrame)




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