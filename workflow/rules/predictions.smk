"""
PACE Prediction Rules
Author: Linyong Shen @ Northwest A&F University
"""

from functools import partial


def _get_hic_params(wildcards):
    """Get Hi-C parameters"""
    if "HiC_file" not in BIOSAMPLES_CONFIG.columns:
        return ""
    
    hic_file = BIOSAMPLES_CONFIG.loc[wildcards.biosample, "HiC_file"]
    
    if pd.isna(hic_file) or not hic_file:
        return ""
    
    hic_type = BIOSAMPLES_CONFIG.loc[wildcards.biosample, "HiC_type"] if "HiC_type" in BIOSAMPLES_CONFIG.columns else "hic"
    hic_resolution = BIOSAMPLES_CONFIG.loc[wildcards.biosample, "HiC_resolution"] if "HiC_resolution" in BIOSAMPLES_CONFIG.columns else 5000
    
    return f"--hic_file {hic_file} --hic_type {hic_type} --hic_resolution {hic_resolution}"


def _get_expression_params(wildcards):
    """Get expression parameters"""
    if "RNA_seq" not in BIOSAMPLES_CONFIG.columns:
        return ""
    
    expr_file = BIOSAMPLES_CONFIG.loc[wildcards.biosample, "RNA_seq"]
    
    if pd.isna(expr_file) or not expr_file:
        return ""
    
    use_weight = config.get('expression', {}).get('enabled', False)
    weight_method = config.get('expression', {}).get('weight_method', 'log')
    min_expression = config.get('expression', {}).get('min_expression', 1.0)
    
    params = f"--expression {expr_file} --min_expression {min_expression}"
    if use_weight:
        params += f" --use_expression_weight --weight_method {weight_method}"
    
    return params


rule create_predictions:
    """Generate PACE enhancer-gene predictions"""
    input:
        enhancers = os.path.join(RESULTS_DIR, "{biosample}", "Neighborhoods", "EnhancerList.txt"),
        genes = os.path.join(RESULTS_DIR, "{biosample}", "Neighborhoods", "GeneList.txt"),
    params:
        output_file = lambda wildcards: os.path.join(RESULTS_DIR, wildcards.biosample, "Predictions", "EnhancerPredictionsAllPutative.tsv.gz"),
        hic_params = _get_hic_params,
        expression_params = _get_expression_params,
        gamma = config['params_predict']['hic_gamma'],
        scale = config['params_predict']['hic_scale'],
        threshold = lambda wildcards: config['params_filter_predictions']['threshold'] or 0.02,
        scripts_dir = SCRIPTS_DIR,
    conda:
        "../envs/pace-env.yml"
    output:
        allPutative = os.path.join(RESULTS_DIR, "{biosample}", "Predictions", "EnhancerPredictionsAllPutative.tsv.gz"),
    resources:
        mem_mb = partial(determine_mem_mb, min_gb=20)
    log:
        os.path.join(RESULTS_DIR, "{biosample}", "logs", "predictions.log")
    shell:
        """
        python {params.scripts_dir}/pace_predict.py \
            --enhancers {input.enhancers} \
            --genes {input.genes} \
            --output {output.allPutative} \
            --hic_gamma {params.gamma} \
            --hic_scale {params.scale} \
            --threshold {params.threshold} \
            {params.hic_params} \
            {params.expression_params} \
            2> {log}
        """


rule filter_predictions:
    """Filter prediction results by threshold"""
    input:
        allPutative = os.path.join(RESULTS_DIR, "{biosample}", "Predictions", "EnhancerPredictionsAllPutative.tsv.gz"),
    params:
        score_column = config['params_filter_predictions']['score_column'],
        threshold = lambda wildcards: determine_threshold(wildcards.biosample),
        only_expressed = "--only_expressed" if config['params_filter_predictions'].get('only_expressed_genes', False) else "",
        scripts_dir = SCRIPTS_DIR,
    conda:
        "../envs/pace-env.yml"
    output:
        enhPredictionsFull = os.path.join(RESULTS_DIR, "{biosample}", "Predictions", f"EnhancerPredictionsFull_{FILTERED_PREDICTION_FILE_FORMAT_TEMPLATE}.tsv"),
        enhPredictionsSlim = os.path.join(RESULTS_DIR, "{biosample}", "Predictions", f"EnhancerPredictions_{FILTERED_PREDICTION_FILE_FORMAT_TEMPLATE}.tsv"),
    resources:
        mem_mb = determine_mem_mb
    log:
        os.path.join(RESULTS_DIR, "{biosample}", "logs", "filter_predictions.log")
    shell:
        """
        python {params.scripts_dir}/pace_filter.py \
            --predictions {input.allPutative} \
            --output {output.enhPredictionsFull} \
            --threshold {params.threshold} \
            --score_column {params.score_column} \
            {params.only_expressed} \
            2> {log}
        """
