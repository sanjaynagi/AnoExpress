import Ano_express as xpress
import pytest
import pandas as pd
import numpy as np
from pandas.testing import assert_frame_equal
from numpy.testing import assert_allclose




@pytest.mark.parametrize(    
        "data_type",
    ["fcs", "pvals", "log2counts"],
    )
@pytest.mark.parametrize(
    "analysis",
    ['gamb_colu', 'gamb_colu_arab', 'gamb_colu_arab_fun', 'fun', 'irtex']
)
def test_load_results_arrays(data_type, analysis):
    df = xpress.load_results_arrays(data_type=data_type, analysis=analysis)
    assert isinstance(df, pd.DataFrame)
    assert not df.empty




@pytest.mark.parametrize('data_type', ["fcs", "pvals", "log2counts"])
@pytest.mark.parametrize("analysis", ['gamb_colu', 'gamb_colu_arab', 'gamb_colu_arab_fun', 'fun'])
@pytest.mark.parametrize("microarray", [True, False])
@pytest.mark.parametrize("gene_id",    [None, 'AGAP006227', ['AGAP006222', 'AGAP006226', 'AGAP006227']])
def test_data(data_type, analysis, gene_id, microarray):

    data_df = xpress.data(
        data_type=data_type,
        analysis=analysis, 
        gene_id= gene_id, 
        microarray=microarray
    )
    assert data_df is not None
    assert not data_df.empty
    assert isinstance(data_df, pd.DataFrame)
#    assert data_df.dtype












@pytest.mark.parametrize(
        "gene_ids",
        ['AGAP006227', ['AGAP006222', 'AGAP006226', 'AGAP006227']]
)
@pytest.mark.parametrize(
    "analysis",
    ['gamb_colu', 'gamb_colu_arab', 'gamb_colu_arab_fun']
)
@pytest.mark.parametrize(
    "microarray",
    [True, False]
)
@pytest.mark.parametrize(
    "plot_type", 
    ['strip', 'boxplot']
)
@pytest.mark.parametrize(
    'sort_by',
    ['median', 'mean', 'agap']
)
def test_plot_gene_expression(gene_ids, analysis, microarray, plot_type, sort_by):

    xpress.plot_gene_expression(
        gene_id=gene_ids, 
        analysis=analysis, 
        microarray=microarray, 
        plot_type=plot_type,
        sort_by=sort_by
        )
    

