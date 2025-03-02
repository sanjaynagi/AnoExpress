import hashlib
import os
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm.notebook import tqdm

gff_url = "https://vectorbase.org/common/downloads/release-68/AgambiaePEST/gff/data/VectorBase-68_AgambiaePEST.gff"
CACHE_DIR = Path(os.path.expanduser("~/.cache/anoexpress"))


def resolve_gene_id(gene_id, analysis, gff_method="vectorbase", use_cache=True):
    """
    Resolve a gene identifier, which may be a coordinate range, a file path, or already a gene ID.

    Parameters
    ----------
    gene_id : str or list
        Gene identifier(s) to resolve
    analysis : str
        Analysis type, used to check compatibility
    gff_method : str, default='vectorbase'
        Method to use for loading GFF data if needed
    use_cache : bool, default=True
        Whether to use cached GFF data if available

    Returns
    -------
    list
        Resolved gene identifiers
    """
    if isinstance(gene_id, str):
        if gene_id.startswith(("2L", "2R", "3L", "3R", "X", "2RL", "3RL")):
            if analysis == "fun":
                assert "Unfortunately the genome feature file does not contain AFUN identifiers, so we cannot subset by genomic span for An. funestus."
            else:

                contig, start_end = gene_id.split(":")
                start, end = start_end.replace(",", "").split("-")
                start, end = int(start), int(end)

                gff = load_gff(
                    query=f"contig == '{contig}' and start <= {end} and end >= {start}",
                    method=gff_method,
                    use_cache=use_cache,
                )
                gene_id = gff.GeneID.to_list()

        elif gene_id.endswith((".tsv", ".txt")):
            gene_id = pd.read_csv(gene_id, sep="\t", header=None).iloc[:, 0].to_list()
        elif gene_id.endswith(".csv"):
            gene_id = pd.read_csv(gene_id, header=None).iloc[:, 0].to_list()
        elif gene_id.endswith(".xlsx"):
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
            ids = gene_annot_df.query(
                f"GO_terms.str.contains('{annot}', na=False)", engine="python"
            )["gene_id"].to_numpy()
        else:
            ids = gene_annot_df.query("domain == @annot")["gene_id"].to_numpy()
        gene_list = np.hstack([gene_list, ids])

    return np.unique(gene_list)


def load_gff(method="vectorbase", override_type=None, query=None, use_cache=True):
    """
    Load genome feature file (GFF) data from either VectorBase or MalariaGEN.

    Parameters
    ----------
    method : str, default='vectorbase'
        Source for the GFF data, either 'vectorbase' or 'malariagen_data'
    override_type : str, optional
        Override the default feature type filter
    query : str, optional
        Query string to filter the GFF data
    use_cache : bool, default=True
        Whether to use cached data if available

    Returns
    -------
    pandas.DataFrame
        DataFrame containing the GFF data
    """
    # Create cache directory if it doesn't exist
    if use_cache and not CACHE_DIR.exists():
        CACHE_DIR.mkdir(parents=True, exist_ok=True)

    if method == "malariagen_data":
        import malariagen_data

        record_type = override_type if override_type else "gene"

        # Generate cache key for malariagen data
        cache_key = f"malariagen_{record_type}"
        if query:
            # Hash the query to create a unique filename
            query_hash = hashlib.md5(query.encode()).hexdigest()[:10]
            cache_key = f"{cache_key}_{query_hash}"
        cache_file = CACHE_DIR / f"{cache_key}.parquet"

        # Try to load from cache
        if use_cache and cache_file.exists():
            try:
                gff = pd.read_parquet(cache_file)
                return gff
            except Exception:
                # If loading from cache fails, continue with normal loading
                pass

        # Normal loading path
        ag3 = malariagen_data.Ag3()
        gff = ag3.genome_features(["2RL", "3RL", "X"])
        gff = gff.query(f"type == '{record_type}'").rename(columns={"ID": "GeneID"})

        # Apply query if provided
        if query:
            gff = gff.query(query)

        # Cache the result
        if use_cache:
            try:
                gff.to_parquet(cache_file)
            except Exception:
                # If caching fails, just continue
                pass

    elif method == "vectorbase":
        record_type = override_type if override_type else "protein_coding_gene"

        # Generate cache key for vectorbase data
        cache_key = f"vectorbase_{record_type}"
        if query:
            # Hash the query to create a unique filename
            query_hash = hashlib.md5(query.encode()).hexdigest()[:10]
            cache_key = f"{cache_key}_{query_hash}"
        cache_file = CACHE_DIR / f"{cache_key}.parquet"

        # Try to load from cache
        if use_cache and cache_file.exists():
            try:
                gff = pd.read_parquet(cache_file)
                return gff
            except Exception:
                # If loading from cache fails, continue with normal loading
                pass

        # Raw GFF download cache
        raw_cache_file = CACHE_DIR / "vectorbase_raw.parquet"

        # Check if raw data is cached
        if use_cache and raw_cache_file.exists():
            try:
                df = pd.read_parquet(raw_cache_file)
            except Exception:
                # If loading raw cache fails, download the data
                df = pd.concat(
                    [
                        chunk
                        for chunk in tqdm(
                            pd.read_csv(
                                gff_url, sep="\t", comment="#", chunksize=10000
                            ),
                            desc="Loading gff data from VectorBase",
                        )
                    ]
                )
                # Cache the raw data
                try:
                    df.to_parquet(raw_cache_file)
                except Exception:
                    pass
        else:
            # Download the data
            df = pd.concat(
                [
                    chunk
                    for chunk in tqdm(
                        pd.read_csv(gff_url, sep="\t", comment="#", chunksize=10000),
                        desc="Loading gff data from VectorBase",
                    )
                ]
            )
            # Cache the raw data
            if use_cache:
                try:
                    df.to_parquet(raw_cache_file)
                except Exception:
                    pass

        # Process the data
        df.columns = [
            "contig",
            "source",
            "type",
            "start",
            "end",
            "na",
            "strand",
            "na2",
            "attributes",
        ]
        df = df.assign(contig=lambda x: x.contig.str.split("_").str.get(1))

        if record_type:  # Fixed bug: was checking 'type' instead of 'record_type'
            df = df.query(f"type == '{record_type}'")

        # may only work for protein_coding_genes
        df = df.assign(
            GeneID=df.attributes.str.split(";", expand=True)
            .iloc[:, 0]
            .str.split("=")
            .str.get(1)
        )

        # combine 2R and 2L, 3R and 3L
        offset_2R = 61545105
        offset_3R = 53200684

        gffs = []
        for contig in tqdm(["2R", "2L", "3R", "3L"], desc="Processing contigs"):
            df_contig = df.query("contig == @contig").copy()
            if contig == "2L":
                df_contig = df_contig.assign(
                    contig="2RL",
                    start=lambda x: x.start + offset_2R,
                    end=lambda x: x.end + offset_2R,
                )
            if contig == "3L":
                df_contig = df_contig.assign(
                    contig="3RL",
                    start=lambda x: x.start + offset_3R,
                    end=lambda x: x.end + offset_3R,
                )
            elif contig in ["3R", "2R"]:
                df_contig = df_contig.assign(contig=lambda x: x.contig + "L")
            gffs.append(df_contig)

        gff = pd.concat(gffs)
        gff = pd.concat([gff, df]).sort_values(["contig", "start", "end"])

        # Apply query if provided
        if query:
            gff = gff.query(query)

        # Cache the processed result
        if use_cache:
            try:
                gff.to_parquet(cache_file)
            except Exception:
                # If caching fails, just continue
                pass
    else:
        raise ValueError(
            f"Unknown method: {method}. Use 'vectorbase' or 'malariagen_data'."
        )

    return gff
