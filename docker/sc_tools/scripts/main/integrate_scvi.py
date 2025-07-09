import os
import argparse
import anndata as ad
import scvi

def integrate_with_scvi(
    adata: ad.AnnData,
    batch_key: str,
) -> tuple[ad.AnnData, scvi.model.SCVI]:
    """
    fit scvi model to AnnData object
    """

    ## fixed parameters
    n_latent = 30
    n_layers = 2
    train_size = 0.85
    scvi_epochs = 300
    batch_size = 1024
    accelerator = "gpu"
    dispersion = "gene-cell"  # "gene"
    gene_likelihood = "zinb"
    latent_distribution = "normal"
    early_stopping = True
    early_stopping_patience = 20
    ## DEPRICATE these training parameters. defaults are good
    # plan_kwargs = {"lr_factor": 0.1, "lr_patience": 20, "reduce_lr_on_plateau": True}

    # integrate the data with `scVI`
    # noise = ["doublet_score", "pct_counts_mt", "pct_counts_rb"]
    categorical_covariate_keys = None
    scvi.model.SCVI.setup_anndata(
        adata,
        layer="counts",
        batch_key=batch_key,
        # continuous_covariate_keys=noise,
        # categorical_covariate_keys=categorical_covariate_keys,
    )

    model = scvi.model.SCVI(
        adata,
        n_layers=n_layers,
        n_latent=n_latent,
        dispersion=dispersion,
        gene_likelihood=gene_likelihood,
    )

    model.train(
        train_size=train_size,
        max_epochs=scvi_epochs,
        early_stopping=early_stopping,
        early_stopping_patience=early_stopping_patience,
        accelerator=accelerator,
        # early_stopping_patience=early_stopping_patience,
        # plan_kwargs=plan_kwargs,
    )

    adata.obsm[args.latent_key] = model.get_latent_representation()  # type: ignore

    return (adata, model)


def main(args: argparse.Namespace):
    # Set the number of data loader workers
    # scvi.settings.dl_num_workers = 2 # Or a value appropriate for your system

    # 0. load data
    adata = ad.read_h5ad(args.adata_input)  # type: ignore
    # 2. process data
    adata, model = integrate_with_scvi(adata, args.batch_key)
    # 3. save the integrated adata and scvi model
    model.save(args.output_scvi_dir, overwrite=True)
    adata.write_h5ad(filename=args.adata_output, compression="gzip")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run scVI integration")
    parser.add_argument(
        "--latent-key",
        type=str,
        required=True,
        help="Latent key to save the scVI latent to",
    )
    parser.add_argument(
        "--batch-key",
        dest="batch_key",
        type=str,
        required=True,
        help="Key in AnnData object for batch information",
    )
    parser.add_argument(
        "--adata-input",
        dest="adata_input",
        type=str,
        required=True,
        help="AnnData object for a dataset",
    )
    parser.add_argument(
        "--adata-output",
        dest="adata_output",
        type=str,
        required=True,
        help="Output file to save AnnData object to",
    )
    parser.add_argument(
        "--output-scvi-dir",
        dest="output_scvi_dir",
        type=str,
        required=True,
        help="Output folder to save `scvi` model",
    )

    args = parser.parse_args()
    main(args)
