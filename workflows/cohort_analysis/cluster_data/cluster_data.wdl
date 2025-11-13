version 1.0

# Perform dataset integration, annotation, and UMAP clustering steps

workflow cluster_data {
	input {
		String cohort_id
		File mmc_adata_object

		# Sample integration
		String scvi_latent_key
		String scanvi_latent_key 
		String scanvi_predictions_key
		String batch_key

		# Clustering parameters
		Int n_neighbors
		Array[Float] leiden_res

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	call integrate_sample_data {
		input:
			cohort_id = cohort_id,
			mmc_adata_object = mmc_adata_object,
			scvi_latent_key = scvi_latent_key,
			batch_key = batch_key,
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	call assign_remaining_cells {
		input:
			cohort_id = cohort_id,
			integrated_adata_object = integrate_sample_data.integrated_adata_object,
			scvi_model_tar_gz = integrate_sample_data.scvi_model_tar_gz, #!FileCoercion
			scanvi_latent_key = scanvi_latent_key,
			scanvi_predictions_key = scanvi_predictions_key,
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	call cluster_cells {
		input:
			cohort_id = cohort_id,
			labeled_cells_adata_object = assign_remaining_cells.labeled_cells_adata_object,
			scvi_latent_key = scvi_latent_key,
			n_neighbors = n_neighbors,
			leiden_res = leiden_res,
			container_registry = container_registry,
			zones = zones
	}

	output {
		File integrated_adata_object = integrate_sample_data.integrated_adata_object
		File scvi_model_tar_gz = integrate_sample_data.scvi_model_tar_gz #!FileCoercion
		File labeled_cells_adata_object = assign_remaining_cells.labeled_cells_adata_object
		File scanvi_model_tar_gz = assign_remaining_cells.scanvi_model_tar_gz #!FileCoercion
		File scanvi_cell_types_parquet = assign_remaining_cells.scanvi_cell_types_parquet #!FileCoercion
		File umap_clustered_adata_object = cluster_cells.umap_clustered_adata_object
	}
}

task integrate_sample_data {
	input {
		String cohort_id
		File mmc_adata_object

		String scvi_latent_key
		String batch_key

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int mem_gb = ceil(size(mmc_adata_object, "GB") * 5 + 20)
	Int disk_size = ceil(size(mmc_adata_object, "GB") * 3 + 50)

	command <<<
		set -euo pipefail

		nvidia-smi

		integrate_scvi \
			--latent-key ~{scvi_latent_key} \
			--batch-key ~{batch_key} \
			--adata-input ~{mmc_adata_object} \
			--adata-output ~{cohort_id}.scvi_integrated.h5ad \
			--output-scvi-dir ~{cohort_id}_scvi_model

		# Model name cannot be changed because scvi models serialization expects a path containing a model.pt object
		tar -czvf "~{cohort_id}.scvi_model.tar.gz" "~{cohort_id}_scvi_model"

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{cohort_id}.scvi_model.tar.gz"
	>>>

	output {
		File integrated_adata_object = "~{cohort_id}.scvi_integrated.h5ad"
		String scvi_model_tar_gz = "~{raw_data_path}/~{cohort_id}.scvi_model.tar.gz"
	}

	runtime {
		docker: "~{container_registry}/sc_tools:1.0.1"
		cpu: 4
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 40
		zones: zones
		gpuType: "nvidia-tesla-t4"
		gpuCount: 1
		nvidiaDriverVersion: "545.23.08" #!UnknownRuntimeKey
	}
}

task assign_remaining_cells {
	input {
		String cohort_id
		File integrated_adata_object
		File scvi_model_tar_gz

		String scanvi_latent_key
		String scanvi_predictions_key

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int mem_gb = ceil(size(integrated_adata_object, "GB") * 20 + 50)
	Int disk_size = ceil(size(integrated_adata_object, "GB") * 3 + 50)

	command <<<
		set -euo pipefail

		nvidia-smi

		mkdir -p "~{cohort_id}_scvi_model"
		tar -xzvf ~{scvi_model_tar_gz} -C "~{cohort_id}_scvi_model" --strip-components=1

		label_scanvi \
			--latent-key ~{scanvi_latent_key} \
			--predictions-key ~{scanvi_predictions_key} \
			--adata-input ~{integrated_adata_object} \
			--scvi-outputs-dir "~{cohort_id}_scvi_model" \
			--adata-output ~{cohort_id}.scanvi_labeled.h5ad \
			--output-scanvi-dir ~{cohort_id}_scanvi_model \
			--output-cell-types-file "~{cohort_id}.scanvi_cell_types.parquet"

		# Model name cannot be changed because scvi models serialization expects a path containing a model.pt object
		tar -czvf "~{cohort_id}.scanvi_model.tar.gz" "~{cohort_id}_scanvi_model"

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{cohort_id}.scanvi_model.tar.gz" \
			-o "~{cohort_id}.scanvi_cell_types.parquet"
	>>>

	output {
		File labeled_cells_adata_object = "~{cohort_id}.scanvi_labeled.h5ad"
		String scanvi_model_tar_gz = "~{raw_data_path}/~{cohort_id}.scanvi_model.tar.gz"
		String scanvi_cell_types_parquet = "~{raw_data_path}/~{cohort_id}.scanvi_cell_types.parquet"
	}

	runtime {
		docker: "~{container_registry}/sc_tools:1.0.1"
		cpu: 4
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 40
		zones: zones
		gpuType: "nvidia-tesla-t4"
		gpuCount: 1
		nvidiaDriverVersion: "545.23.08" #!UnknownRuntimeKey
	}
}

task cluster_cells {
	input {
		String cohort_id
		File labeled_cells_adata_object

		String scvi_latent_key

		Int n_neighbors
		Array[Float] leiden_res

		String container_registry
		String zones
	}

	Int mem_gb = ceil(size(labeled_cells_adata_object, "GB") * 10 + 20)
	Int disk_size = ceil(size([labeled_cells_adata_object], "GB") * 6 + 50)

	command <<<
		set -euo pipefail

		clustering_umap \
			--latent-key ~{scvi_latent_key} \
			--n-neighbors ~{n_neighbors} \
			--leiden-res ~{sep=' ' leiden_res} \
			--adata-input ~{labeled_cells_adata_object} \
			--adata-output ~{cohort_id}.umap_clustered.h5ad
	>>>

	output {
		File umap_clustered_adata_object = "~{cohort_id}.umap_clustered.h5ad"
	}

	runtime {
		docker: "~{container_registry}/sc_tools:1.0.1"
		cpu: 16
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 40
		zones: zones
	}
}
