import os
import argparse

import anndata as ad

# import scvi.model as scvi_model
import scvi


SCANVI_LATENT_KEY = "X_scANVI"
SCANVI_PREDICTIONS_KEY = "C_scANVI"


def label_with_scanvi(
    adata: ad.AnnData, model: scvi.model.SCVI
) -> tuple[ad.AnnData, scvi.model.SCANVI]:
    """
    fit scANVI model to AnnData object
    """
    scanvi_epochs = 300
    batch_size = 1024
    accelerator = "gpu"
    dispersion = "gene-cell"  # "gene"
    gene_likelihood = "zinb"
    latent_distribution = "normal"
    early_stopping = True
    early_stopping_patience = 20

    print("scanvi model from scvi")
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        model, adata=adata, labels_key="cell_type", unlabeled_category="Unknown"
    )

    print("train scanvi model")
    scanvi_model.train(
        accelerator=accelerator,
        max_epochs=scanvi_epochs,
        early_stopping=early_stopping,
        early_stopping_patience=early_stopping_patience,
    )

    print("scanvi latents and predictions")

    adata.obsm[SCANVI_LATENT_KEY] = scanvi_model.get_latent_representation(adata)
    adata.obs[SCANVI_PREDICTIONS_KEY] = scanvi_model.predict(adata)

    return (adata, scanvi_model)


def main(args: argparse.Namespace):
    """
    basic logic with args as input

    """

    # Set the number of data loader workers
    # scvi.settings.dl_num_workers = 2 # Or a value appropriate for your system

    # 0. load data
    adata = ad.read_h5ad(args.adata_input)  # type: ignore

    model_path = args.output_scvi_dir
    model = scvi.model.SCVI.load(model_path, adata)

    print("scanvi latents and predictions")
    # 4. get scANVI model
    adata, scanvi_model = label_with_scanvi(adata, model)
    # 5. save the integrated adata and scvi model
    scanvi_model.save(args.output_scanvi_dir, overwrite=True)

    # # 6. save the latent space
    # fname = args.adata_output
    # fname = fname.replace("integrated", "scanvi")
    # adata.write_h5ad(filename=fname, compression="gzip")

    # adata = ad.read_h5ad(args.adata_output)  # type: ignore

    # # 7. save the cell types to feather
    # adata.obs[[SCANVI_PREDICTIONS_KEY]].to_feather(args.cell_type_output)
    # 7. save the cell types to parquet
    adata.obs[[SCANVI_PREDICTIONS_KEY]].to_parquet(args.cell_type_output)
    adata.write_h5ad(filename=args.adata_output, compression="gzip")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run scVI integration")
    # parser.add_argument(
    #     "--latent-key",
    #     dest="latent_key",
    #     type=str,
    #     default="X_scvi",
    #     help="Latent key to save the scvi latent to",
    # )
    # parser.add_argument(
    #     "--scanvi-latent-key",
    #     dest="scanvi_latent_key",
    #     type=str,
    #     default="X_scanvi",
    #     help="Latent key to save the scanvi latent to",
    # )
    parser.add_argument(
        "--batch-key",
        dest="batch_key",
        type=str,
        help="Key in AnnData object for batch information",
    )
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
    parser.add_argument(
        "--output-scvi-dir",
        dest="output_scvi_dir",
        type=str,
        help="Output folder to save `scvi` model",
    )
    parser.add_argument(
        "--output-scanvi-dir",
        dest="output_scanvi_dir",
        type=str,
        help="Output folder to save `scanvi` model",
    )

    parser.add_argument(
        "--output-cell-types-file",
        dest="cell_type_output",
        type=str,
        help="Output file to write celltypes to",
    )

    # TODO: optional scvi arguments
    args = parser.parse_args()
    main(args)
