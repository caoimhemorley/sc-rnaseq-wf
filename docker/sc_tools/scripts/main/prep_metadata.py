import argparse
import sys
import scanpy as sc
from anndata import AnnData

sys.path.append("/opt/scripts/utility")
from helpers import anndata_from_h5


def prep_metadata(
    adata: AnnData, sample_id: str, batch_id: str, dataset_id: str, team_id: str
):
    # Set CPUs to use for parallel computing
    sc._settings.ScanpyConfig.n_jobs = -1

    adata.var_names_make_unique()
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["rb"] = adata.var_names.str.startswith(("RPL", "RPS"))
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["rb", "mt"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )

    sc.external.pp.scrublet(adata)
    # Keep the metadata clean and drop predicted_doublet; rely on the probability of doublet (i.e., doublet_score)
    adata.obs.drop("predicted_doublet", axis=1, inplace=True)

    # Add metadata
    adata.obs["sample"] = sample_id
    adata.obs["batch"] = batch_id
    adata.obs["team"] = team_id
    adata.obs["dataset"] = dataset_id
    adata.obs["batch_id"] = f"{team_id}_{dataset_id}_{batch_id}"
    return adata


def main(args: argparse.Namespace):
    # Load the data from cellbender output
    adata = anndata_from_h5(args.adata_input)
    adata = prep_metadata(adata, args.sample_id, args.batch, args.dataset, args.team)
    adata.write_h5ad(filename=args.adata_output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Preprocess by calculating doublets and QC metrics")
    parser.add_argument(
        "--adata-input",
        type=str,
        required=True,
        help="AnnData object for a dataset",
    )
    parser.add_argument(
        "--sample-id",
        type=str,
        required=True,
        help="Sample/dataset ID"
    )
    parser.add_argument(
        "--batch",
        type=str,
        required=True,
        help="Batch from which the sample/dataset originated",
    )
    parser.add_argument(
        "--dataset",
        type=str,
        required=True,
        help="Dataset ID"
    )
    parser.add_argument(
        "--team",
        type=str,
        required=True,
        help="Team ID"
    )
    parser.add_argument(
        "--adata-output",
        type=str,
        required=True,
        help="Output file to save AnnData object to",
    )

    args = parser.parse_args()
    main(args)
