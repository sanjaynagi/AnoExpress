import pandas as pd
import numpy as np
from tqdm.notebook import tqdm


gff_url =  'https://vectorbase.org/common/downloads/release-68/AgambiaePEST/gff/data/VectorBase-68_AgambiaePEST.gff'


def resolve_gene_id(gene_id, analysis, gff_method='malariagen_data'):
    
    if isinstance(gene_id, str):
      if gene_id.startswith(('2L', '2R', '3L', '3R', 'X', '2RL', '3RL')):
        if analysis == 'fun':
          assert "Unfortunately the genome feature file does not contain AFUN identifiers, so we cannot subset by genomic span for An. funestus."
        else:

          contig, start_end = gene_id.split(':')
          start, end = start_end.replace(",", "").split('-')
          start, end = int(start), int(end)

          gff = load_gff(query=f"contig == '{contig}' and start <= {end} and end >= {start}", method=gff_method)
          gene_id = gff.GeneID.to_list()

      elif gene_id.endswith(('.tsv', '.txt')):
          gene_id = pd.read_csv(gene_id, sep="\t", header=None).iloc[:, 0].to_list()
      elif gene_id.endswith('.csv'):
          gene_id = pd.read_csv(gene_id, header=None).iloc[:, 0].to_list()
      elif gene_id.endswith('.xlsx'):
          gene_id = pd.read_excel(gene_id, header=None).iloc[:, 0].to_list()
      
    return gene_id



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


def load_gff(method='malariagen_data', override_type=None, query=None):

    if method == 'malariagen_data':
        import malariagen_data
        record_type = override_type if override_type else 'gene'
        
        ag3 = malariagen_data.Ag3()
        gff = ag3.genome_features(['2RL', '3RL', 'X'])
        gff = gff.query(f"type == '{record_type}'").rename(columns={'ID':'GeneID'})

    elif method == 'vectorbase':
        record_type = override_type if override_type else 'protein_coding_gene'
        
        df = pd.concat([chunk for chunk in tqdm(pd.read_csv(gff_url, sep="\t", comment="#", chunksize=10000), desc='Loading gff data from VectorBase')])
        df.columns = ['contig', 'source', 'type', 'start', 'end', 'na', 'strand', 'na2', 'attributes']
        df = df.assign(contig=lambda x: x.contig.str.split("_").str.get(1))
        
        if type:
         df = df.query(f"type == '{record_type}'")
    
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