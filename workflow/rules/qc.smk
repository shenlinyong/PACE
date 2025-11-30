"""
Quality Control Rules for PACE
Author: Linyong Shen @ Northwest A&F University
"""

rule generate_qc_plot_and_summary:
    """Generate QC plots and summary"""
    input:
        enhPredictionsFull = os.path.join(RESULTS_DIR, "{biosample}", "Predictions", f"EnhancerPredictionsFull_{FILTERED_PREDICTION_FILE_FORMAT_TEMPLATE}.tsv"),
    params:
        output_dir = os.path.join(RESULTS_DIR, "{biosample}", "Metrics"),
        scripts_dir = SCRIPTS_DIR,
        sample_name = lambda wildcards: wildcards.biosample,
    conda:
        "../envs/pace-env.yml"
    output:
        qc_summary = os.path.join(RESULTS_DIR, "{biosample}", "Metrics", f"QCSummary_{FILTERED_PREDICTION_FILE_FORMAT_TEMPLATE}.tsv"),
        qc_plots = os.path.join(RESULTS_DIR, "{biosample}", "Metrics", f"QCPlots_{FILTERED_PREDICTION_FILE_FORMAT_TEMPLATE}.pdf")
    resources:
        mem_mb = determine_mem_mb
    log:
        os.path.join(RESULTS_DIR, "{biosample}", "logs", "qc.log")
    shell:
        """
        mkdir -p {params.output_dir}
        
        python {params.scripts_dir}/pace_metrics.py \
            --predictions {input.enhPredictionsFull} \
            --output_dir {params.output_dir} \
            --sample_name {params.sample_name} \
            2> {log}
        """
