version 1.0

# Perform dataset integration, umap, clustering, and annotation steps

workflow cluster_data {
	input {
		String cohort_id
		File mmc_adata_object

		String scvi_latent_key
		String scanvi_latent_key 
		String scanvi_predictions_key
		String batch_key

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
			labeled_cells_adata_object = assign_remaining_cells.labeled_cells_adata_object,
			scvi_latent_key = scvi_latent_key,
			container_registry = container_registry,
			zones = zones
	}

	output {
		File integrated_adata_object = integrate_sample_data.integrated_adata_object
		File scvi_model_tar_gz = integrate_sample_data.scvi_model_tar_gz #!FileCoercion
		File labeled_cells_adata_object = assign_remaining_cells.labeled_cells_adata_object
		File scanvi_model_tar_gz = assign_remaining_cells.scanvi_model_tar_gz #!FileCoercion
		File scanvi_cell_types_parquet_gzip = assign_remaining_cells.scanvi_cell_types_parquet_gzip #!FileCoercion
		File umap_cluster_adata_object = cluster_cells.umap_cluster_adata_object
	}
}

task integrate_sample_data {
	input {
		String cohort_id
		File mmc_adata_object

		String batch_key
		String scvi_latent_key

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

		python3 /opt/scripts/main/integrate_scvi.py \
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
		docker: "~{container_registry}/sc_tools:1.0.0"
		cpu: 2
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

	Int mem_gb = ceil(size(integrated_adata_object, "GB") * 5 + 20)
	Int disk_size = ceil(size(integrated_adata_object, "GB") * 3 + 50)

	command <<<
		set -euo pipefail

		nvidia-smi

		mkdir -p "~{cohort_id}_scvi_model"
		tar -xzvf ~{scvi_model_tar_gz} -C "~{cohort_id}_scvi_model" --strip-components=1

		python3 /opt/scripts/main/label_scanvi.py \
			--latent-key ~{scanvi_latent_key} \
			--predictions-key ~{scanvi_predictions_key} \
			--adata-input ~{integrated_adata_object} \
			--scvi-outputs-dir "~{cohort_id}_scvi_model" \
			--adata-output ~{cohort_id}.scanvi_labeled.h5ad \
			--output-scanvi-dir ~{cohort_id}_scanvi_model \
			--output-cell-types-file "~{cohort_id}.scanvi_cell_types.parquet.gzip"

		# Model name cannot be changed because scvi models serialization expects a path containing a model.pt object
		tar -czvf "~{cohort_id}.scanvi_model.tar.gz" "~{cohort_id}_scanvi_model"

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{cohort_id}.scanvi_model.tar.gz" \
			-o "~{cohort_id}.scanvi_cell_types.parquet.gzip"
	>>>

	output {
		File labeled_cells_adata_object = "~{cohort_id}.scanvi_labeled.h5ad"
		String scanvi_model_tar_gz = "~{raw_data_path}/~{cohort_id}.scanvi_model.tar.gz"
		String scanvi_cell_types_parquet_gzip = "~{raw_data_path}/~{cohort_id}.scanvi_cell_types.parquet.gzip"
	}

	runtime {
		docker: "~{container_registry}/sc_tools:1.0.0"
		cpu: 2
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
		File labeled_cells_adata_object

		String scvi_latent_key

		String container_registry
		String zones
	}

	String integrated_adata_object_basename = basename(labeled_cells_adata_object, ".h5ad")
	Int mem_gb = ceil(size(labeled_cells_adata_object, "GB") * 8.7 + 20)
	Int disk_size = ceil(size([labeled_cells_adata_object], "GB") * 6 + 50)

	command <<<
		set -euo pipefail

		python3 /opt/scripts/main/clustering_umap.py \
			--adata-input ~{labeled_cells_adata_object} \
			--adata-output ~{integrated_adata_object_basename}.umap_cluster.h5ad \
			--latent-key ~{scvi_latent_key}
	>>>

	output {
		File umap_cluster_adata_object = "~{integrated_adata_object_basename}.umap_cluster.h5ad"
	}

	runtime {
		docker: "~{container_registry}/sc_tools:1.0.0"
		cpu: 16
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 40
		zones: zones
	}
}
