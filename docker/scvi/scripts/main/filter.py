# import muon.pp.filter_obs as filter_obs ???
# import muon as mu
import scanpy as sc
import argparse
import pandas as pd
import sys
from anndata import AnnData

sys.path.append("/opt/scripts/utility")
from helpers import update_validation_metrics


def filter_adata(adata_fname: str) -> AnnData:
    # TODO: make these cutoffs arguments...
    ## Fixed Parameters
    
    pct_counts_mt_max = 10
    doublet_score_max = 0.2
    total_counts_limits = [100, 100000]
    n_genes_by_counts_limits = [100, 10000]

    print(f"----load adata")

    # 0. load data
    adata_backed = sc.read_h5ad(args.adata_input, backed='r')  # type: ignore
    print(f"read {args.adata_input}")
    
    # Step 2: Choose your filter condition on the `.obs` or `.var` (must be present in backed mode)
    # For example, filter cells with a certain metadata label
    # muon api is better than the scanpyt api for this...
    # mu.pp.filter_obs(adata, "pct_counts_mt", lambda x: x <= pct_counts_mt_max)
    keep_mt = adata_backed.obs['pct_counts_mt'] <= pct_counts_mt_max

    print(f"----   filter 2")
    # mu.pp.filter_obs(adata, "doublet_score", lambda x: x < doublet_score_max)
    keep_doublet = adata_backed.obs['doublet_score'] < doublet_score_max

    print(f"----      filter 3")
    # mu.pp.filter_obs(
    #     adata,
    #     "total_counts",
    #     lambda x: (x >= total_counts_limits[0]) & (x <= total_counts_limits[1]),
    # )
    keep_total_counts = (adata_backed.obs['total_counts'] >= total_counts_limits[0]) & (adata_backed.obs['total_counts'] <= total_counts_limits[1])

    print(f"----         filter 4")
    # mu.pp.filter_obs(
    #     adata,
    #     "n_genes_by_counts",
    #     lambda x: (x >= n_genes_by_counts_limits[0])
    #     & (x <= n_genes_by_counts_limits[1]),
    # )
    keep_n_genes_by_counts = (adata_backed.obs['n_genes_by_counts'] >= n_genes_by_counts_limits[0]) & (adata_backed.obs['n_genes_by_counts'] <= n_genes_by_counts_limits[1])


    keep_cells = keep_mt & keep_doublet & keep_total_counts & keep_n_genes_by_counts

    
    # Step 3: Read the filtered subset into memory
    print(f"APPLY FILTER")
    adata_filtered = adata_backed[keep_cells].to_memory()
    
    # Step 4: Optional - Close the backed object to release file handles
    adata_backed.file.close()
    
    return adata_filtered


def main(args: argparse.Namespace):
    """
    basic logic with args as input

    """
    # Set CPUs to use for parallel computing
    sc._settings.ScanpyConfig.n_jobs = -1

    
    # 1. filter data
    adata = filter_adata(args.adata_input)
    print(f"---- filtered -----")
    
    # 2. save the filtered adata
    #    save the filtered adata
    # adata.write_h5ad(filename=args.adata_output, compression="gzip")
    adata.write_h5ad(filename=args.adata_output)
    print(f"wrote {args.adata_output}")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Filter")
    parser.add_argument(
        "--adata-input",
        dest="adata_input",
        type=str,
        help="AnnData object for a dataset",
    )
    parser.add_argument(
        "--adata-output",
        dest="adata_output",
        type=str,
        help="Output file to save AnnData object to",
    )
    # TODO: add filter parameters as arguments

    args = parser.parse_args()
    main(args)
