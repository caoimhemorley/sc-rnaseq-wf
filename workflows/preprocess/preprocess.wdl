
version 1.0
import "../structs.wdl"

workflow preprocess {
	input {
		String team_id
		String dataset_id
		String dataset_doi_url
		Array[Sample] samples

		Boolean multimodal_sc_data
		File cellranger_reference_data
		Float cellbender_fpr

		String workflow_name
		String workflow_version
		String workflow_release
		String run_timestamp
		String raw_data_path_prefix
		String billing_project
		String container_registry
		String zones
	}

	String sub_workflow_name = "preprocess"
	String cellranger_task_version = "1.1.0"
	String cellbender_task_version = "1.0.0"
	String adata_task_version = "1.1.0"

	Array[Array[String]] workflow_info = [[run_timestamp, workflow_name, workflow_version, workflow_release]]

	String workflow_raw_data_path_prefix = "~{raw_data_path_prefix}/~{sub_workflow_name}"
	String cellranger_raw_data_path = "~{workflow_raw_data_path_prefix}/cellranger/~{cellranger_task_version}"
	String cellbender_raw_data_path = "~{workflow_raw_data_path_prefix}/remove_technical_artifacts/~{cellbender_task_version}"
	String adata_raw_data_path = "~{workflow_raw_data_path_prefix}/counts_to_adata/~{adata_task_version}"

	scatter (sample_object in samples) {
		String cellranger_count_output = "~{cellranger_raw_data_path}/~{sample_object.sample_id}.raw_feature_bc_matrix.h5"
		String cellbender_count_output = "~{cellbender_raw_data_path}/~{sample_object.sample_id}.cellbender.h5"
		String initial_adata_object_output = "~{adata_raw_data_path}/~{sample_object.sample_id}.cleaned_unfiltered.h5ad"
	}

	call check_output_files_exist {
		input:
			cellranger_count_output_files = cellranger_count_output,
			remove_technical_artifacts_output_files = cellbender_count_output,
			initial_adata_object_output_files = initial_adata_object_output,
			billing_project = billing_project,
			zones = zones
	}

	scatter (index in range(length(samples))) {
		Sample sample = samples[index]

		Array[String] project_sample_id = [team_id, sample.sample_id, dataset_doi_url]

		String cellranger_count_complete = check_output_files_exist.sample_preprocessing_complete[index][0]
		String cellbender_remove_background_complete = check_output_files_exist.sample_preprocessing_complete[index][1]
		String initial_adata_object_complete = check_output_files_exist.sample_preprocessing_complete[index][2]

		String cellranger_raw_counts = "~{cellranger_raw_data_path}/~{sample.sample_id}.raw_feature_bc_matrix.h5"
		String cellranger_filtered_counts = "~{cellranger_raw_data_path}/~{sample.sample_id}.filtered_feature_bc_matrix.h5"
		String cellranger_molecule_info = "~{cellranger_raw_data_path}/~{sample.sample_id}.molecule_info.h5"
		String cellranger_metrics_summary_csv = "~{cellranger_raw_data_path}/~{sample.sample_id}.metrics_summary.csv"

		# NEW (BAM)
		String cellranger_genome_bam = "~{cellranger_raw_data_path}/~{sample.sample_id}.possorted_genome.bam"
		String cellranger_genome_bam_bai = "~{cellranger_raw_data_path}/~{sample.sample_id}.possorted_genome.bam.bai"

		if (cellranger_count_complete == "false") {
			call cellranger_count {
				input:
					sample_id = sample.sample_id,
					fastq_R1s = sample.fastq_R1s,
					fastq_R2s = sample.fastq_R2s,
					fastq_I1s = sample.fastq_I1s,
					fastq_I2s = sample.fastq_I2s,
					multimodal_sc_data = multimodal_sc_data,
					cellranger_reference_data = cellranger_reference_data,
					raw_data_path = cellranger_raw_data_path,
					workflow_info = workflow_info,
					billing_project = billing_project,
					container_registry = container_registry,
					zones = zones
			}
		}

		File raw_counts_output = select_first([cellranger_count.raw_counts, cellranger_raw_counts])
		File filtered_counts_output = select_first([cellranger_count.filtered_counts, cellranger_filtered_counts])
		File molecule_info_output = select_first([cellranger_count.molecule_info, cellranger_molecule_info])
		File metrics_summary_csv_output = select_first([cellranger_count.metrics_summary_csv, cellranger_metrics_summary_csv])

		# NEW (BAM)
		File genome_bam_output = select_first([cellranger_count.genome_bam, cellranger_genome_bam])
		File genome_bam_bai_output = select_first([cellranger_count.genome_bam_bai, cellranger_genome_bam_bai])

		String cellbender_report_html = "~{cellbender_raw_data_path}/~{sample.sample_id}.cellbender_report.html"
		String cellbender_removed_background_counts = "~{cellbender_raw_data_path}/~{sample.sample_id}.cellbender.h5"
		String cellbender_filtered_removed_background_counts = "~{cellbender_raw_data_path}/~{sample.sample_id}.cellbender_filtered.h5"
		String cellbender_cell_barcodes_csv = "~{cellbender_raw_data_path}/~{sample.sample_id}.cellbender_cell_barcodes.csv"
		String cellbender_graph_pdf = "~{cellbender_raw_data_path}/~{sample.sample_id}.cellbender.pdf"
		String cellbender_log = "~{cellbender_raw_data_path}/~{sample.sample_id}.cellbender.log"
		String cellbender_metrics_csv = "~{cellbender_raw_data_path}/~{sample.sample_id}.cellbender_metrics.csv"
		String cellbender_posterior_probability = "~{cellbender_raw_data_path}/~{sample.sample_id}.cellbend_posterior.h5"

		if (cellbender_remove_background_complete == "false") {
			call remove_technical_artifacts {
				input:
					sample_id = sample.sample_id,
					raw_counts = raw_counts_output,
					cellbender_fpr = cellbender_fpr,
					raw_data_path = cellbender_raw_data_path,
					workflow_info = workflow_info,
					billing_project = billing_project,
					container_registry = container_registry,
					zones = zones
			}
		}

		File report_html_output = select_first([remove_technical_artifacts.report_html, cellbender_report_html])
		File removed_background_counts_output = select_first([remove_technical_artifacts.removed_background_counts, cellbender_removed_background_counts])
		File filtered_removed_background_counts_output = select_first([remove_technical_artifacts.filtered_removed_background_counts, cellbender_filtered_removed_background_counts])
		File cell_barcodes_csv_output = select_first([remove_technical_artifacts.cell_barcodes_csv, cellbender_cell_barcodes_csv])
		File graph_pdf_output = select_first([remove_technical_artifacts.graph_pdf, cellbender_graph_pdf])
		File log_output = select_first([remove_technical_artifacts.log, cellbender_log])
		File metrics_csv_output = select_first([remove_technical_artifacts.metrics_csv, cellbender_metrics_csv])
		File posterior_probability_output = select_first([remove_technical_artifacts.posterior_probability, cellbender_posterior_probability])

		String preprocessed_adata_object = "~{adata_raw_data_path}/~{sample.sample_id}.cleaned_unfiltered.h5ad"

		if (initial_adata_object_complete == "false") {
			call counts_to_adata {
				input:
					sample_id = sample.sample_id,
					batch = select_first([sample.batch]),
					team_id = team_id,
					dataset_id = dataset_id,
					cellbender_counts = removed_background_counts_output,
					raw_data_path = adata_raw_data_path,
					workflow_info = workflow_info,
					billing_project = billing_project,
					container_registry = container_registry,
					zones = zones
			}
		}

		File preprocessed_adata_object_output = select_first([counts_to_adata.initial_adata_object, preprocessed_adata_object])
	}

	output {
		Array[Array[String]] project_sample_ids = project_sample_id

		Array[File] raw_counts = raw_counts_output
		Array[File] filtered_counts = filtered_counts_output
		Array[File] molecule_info = molecule_info_output
		Array[File] metrics_summary_csv = metrics_summary_csv_output

		# NEW (BAM)
		Array[File] genome_bam = genome_bam_output
		Array[File] genome_bam_bai = genome_bam_bai_output

		Array[File] report_html = report_html_output
		Array[File] removed_background_counts = removed_background_counts_output
		Array[File] filtered_removed_background_counts = filtered_removed_background_counts_output
		Array[File] cell_barcodes_csv = cell_barcodes_csv_output
		Array[File] graph_pdf = graph_pdf_output
		Array[File] log = log_output
		Array[File] metrics_csv = metrics_csv_output
		Array[File] posterior_probability = posterior_probability_output

		Array[File] initial_adata_object = preprocessed_adata_object_output
	}
}

task cellranger_count {
	input {
		String sample_id
		Array[File] fastq_R1s
		Array[File] fastq_R2s
		Array[File] fastq_I1s
		Array[File] fastq_I2s
		Boolean multimodal_sc_data
		File cellranger_reference_data
		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	String cellranger_arc_chemistry_flag = if multimodal_sc_data then "--chemistry=ARC-v1" else ""

	Int threads = 16
	Int mem_gb = 48
	Int disk_size = ceil((size(cellranger_reference_data, "GB") + size(flatten([fastq_R1s, fastq_R2s, fastq_I1s, fastq_I2s]), "GB")) * 4 + 50)

	command <<<
		set -euo pipefail

		mkdir cellranger_refdata
		tar -zxvf ~{cellranger_reference_data} -C cellranger_refdata --strip-components 1

		mkdir fastqs
		while read -r fastq || [[ -n "${fastq}" ]]; do
			validated_fastq_name=$(fix_fastq_names --fastq "${fastq}" --sample-id "~{sample_id}")
			ln -s "${fastq}" "fastqs/${validated_fastq_name}"
		done < <(cat \
			~{write_lines(fastq_R1s)} \
			~{write_lines(fastq_R2s)} \
			~{write_lines(fastq_I1s)} \
			~{write_lines(fastq_I2s)})

		cellranger count \
			--id=~{sample_id} \
			--transcriptome="$(pwd)/cellranger_refdata" \
			--fastqs="$(pwd)/fastqs" \
			--localcores ~{threads} \
			--localmem ~{mem_gb - 4} \
			~{cellranger_arc_chemistry_flag}

		mv ~{sample_id}/outs/raw_feature_bc_matrix.h5 ~{sample_id}.raw_feature_bc_matrix.h5
		mv ~{sample_id}/outs/filtered_feature_bc_matrix.h5 ~{sample_id}.filtered_feature_bc_matrix.h5
		mv ~{sample_id}/outs/molecule_info.h5 ~{sample_id}.molecule_info.h5
		mv ~{sample_id}/outs/metrics_summary.csv ~{sample_id}.metrics_summary.csv

		# NEW (BAM)
		mv ~{sample_id}/outs/possorted_genome_bam.bam ~{sample_id}.possorted_genome.bam
		mv ~{sample_id}/outs/possorted_genome_bam.bam.bai ~{sample_id}.possorted_genome.bam.bai

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{sample_id}.raw_feature_bc_matrix.h5" \
			-o "~{sample_id}.filtered_feature_bc_matrix.h5" \
			-o "~{sample_id}.molecule_info.h5" \
			-o "~{sample_id}.metrics_summary.csv" \
			-o "~{sample_id}.possorted_genome.bam" \
			-o "~{sample_id}.possorted_genome.bam.bai"
	>>>

	output {
		File raw_counts = "~{raw_data_path}/~{sample_id}.raw_feature_bc_matrix.h5"
		File filtered_counts = "~{raw_data_path}/~{sample_id}.filtered_feature_bc_matrix.h5"
		File molecule_info = "~{raw_data_path}/~{sample_id}.molecule_info.h5"
		File metrics_summary_csv = "~{raw_data_path}/~{sample_id}.metrics_summary.csv"
		File genome_bam = "~{raw_data_path}/~{sample_id}.possorted_genome.bam"
		File genome_bam_bai = "~{raw_data_path}/~{sample_id}.possorted_genome.bam.bai"
	}

	runtime {
		docker: "~{container_registry}/cellranger:7.1.0"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 40
		zones: zones
	}
}
