"""
Utility functions and configuration loading for PACE
"""

import pandas as pd
import time
import random
from pandas.errors import EmptyDataError


class InvalidConfig(Exception):
    """Configuration error exception"""
    pass


wildcard_constraints:
    threshold=r"\d+\.\d+",
    separator=r".{0}|_",
    other_flags=r".{0}|[^0-9]+"

FILTERED_PREDICTION_FILE_FORMAT_TEMPLATE = "threshold{threshold}{separator}{other_flags}"
DEFAULT_THRESHOLD = 0.02
MAX_MEM_MB = 250 * 1000  # 250GB


def determine_mem_mb(wildcards, input, attempt, min_gb=8):
    """Dynamically calculate memory requirements"""
    input_size_mb = input.size_mb
    if ".gz" in str(input):
        input_size_mb *= 8
    attempt_multiplier = 2 ** (attempt - 1)
    mem_to_use_mb = attempt_multiplier * max(4 * input_size_mb, min_gb * 1000)
    return min(mem_to_use_mb, MAX_MEM_MB)


def make_paths_absolute(obj, base_path):
    """Convert relative paths to absolute paths"""
    if isinstance(obj, dict):
        for key, value in obj.items():
            obj[key] = make_paths_absolute(value, base_path)
    elif isinstance(obj, str):
        new_file = os.path.join(base_path, obj) if base_path else obj
        if os.path.exists(new_file):
            return new_file
    return obj


def determine_threshold(biosample):
    """Determine ABC score threshold"""
    config_threshold = config["params_filter_predictions"]["threshold"]
    if config_threshold:
        return config_threshold
    
    biosample_row = BIOSAMPLES_CONFIG[BIOSAMPLES_CONFIG["biosample"] == biosample].iloc[0]
    hic_type = biosample_row["HiC_type"]
    
    if hic_type is None:
        hic_type = "powerlaw"
    elif hic_type == "avg":
        hic_type = "avg_hic"
    elif hic_type == "hic":
        hic_type = "intact_hic"
    
    matching_row = ABC_THRESHOLDS[
        (ABC_THRESHOLDS["accessibility"] == biosample_row["default_accessibility_feature"])
        & (ABC_THRESHOLDS["has_h3k27ac"] == bool(biosample_row["H3K27ac"]))
        & (ABC_THRESHOLDS["hic_type"] == hic_type)
    ]
    
    if len(matching_row) == 0:
        print(f"Threshold not found for {biosample}, using default: {DEFAULT_THRESHOLD}")
        threshold = DEFAULT_THRESHOLD
    else:
        threshold = matching_row.iloc[0]["threshold"]
    
    return threshold


def determine_filtered_prediction_file_format(threshold, config):
    """Determine filtered prediction file format"""
    include_self_promoter = config['params_filter_predictions']['include_self_promoter']
    only_expressed_genes = config['params_filter_predictions']['only_expressed_genes']
    
    if include_self_promoter or only_expressed_genes:
        separator = '_'
        other_flags = []
        if include_self_promoter:
            other_flags.append('self_promoter')
        if only_expressed_genes:
            other_flags.append('only_expr_genes')
        other_flags = "__".join(other_flags)
    else:
        separator = ''
        other_flags = ''
    
    return FILTERED_PREDICTION_FILE_FORMAT_TEMPLATE.format(
        threshold=threshold, separator=separator, other_flags=other_flags
    )


def enable_retry(func, func_args={}, max_attempts=3, delay=0.5):
    """Execute function with retry logic"""
    for attempt in range(max_attempts):
        try:
            return func(**func_args)
        except Exception as e:
            if attempt == max_attempts - 1:
                raise
            sleep_time = delay + random.uniform(0, 0.5)
            time.sleep(sleep_time)
    return None


def load_biosamples_config(config):
    """Load sample configuration file"""
    biosamples_config = enable_retry(
        pd.read_csv,
        func_args={'filepath_or_buffer': config["biosamplesTable"], 'sep': "\t", 'comment': '#'}
    ).replace([np.nan], [None]).set_index("biosample", drop=False)
    
    biosamples_config["HiC_resolution"] = biosamples_config["HiC_resolution"].replace([None], [0]).astype(int)
    _validate_biosamples_config(biosamples_config)
    _configure_tss_and_gene_files(biosamples_config)
    
    return biosamples_config


def load_abc_thresholds(config):
    """Load ABC threshold file"""
    file = config["ref"]["abc_thresholds"]
    return pd.read_csv(file, sep='\t')


def get_accessibility_files(wildcards):
    """Get accessibility data files"""
    files = BIOSAMPLES_CONFIG.loc[wildcards.biosample, "DHS"] or BIOSAMPLES_CONFIG.loc[wildcards.biosample, "ATAC"]
    return files.split(",")


def _validate_accessibility_feature(row: pd.Series):
    """Validate accessibility data configuration"""
    if row["DHS"] and row["ATAC"]:
        raise InvalidConfig("Can only specify one of DHS or ATAC per sample")
    if not (row["DHS"] or row["ATAC"]):
        raise InvalidConfig("Must provide either DHS or ATAC accessibility file")


def _validate_hic_info(row: pd.Series):
    """Validate Hi-C data configuration"""
    if row["HiC_file"]:
        if not (row["HiC_type"] and row["HiC_resolution"]):
            raise InvalidConfig("Must provide HiC_type and HiC_resolution with HiC_file")
        if row["HiC_resolution"] != 5000:
            raise InvalidConfig("Only 5kb resolution Hi-C data is supported")


def _validate_biosamples_config(biosamples_config):
    """Validate sample configuration"""
    for _, row in biosamples_config.iterrows():
        _validate_hic_info(row)
        _validate_accessibility_feature(row)


def _configure_tss_and_gene_files(biosamples_config):
    """Configure TSS and gene files"""
    TSS_files = []
    gene_files = []
    
    for sample in biosamples_config['biosample']:
        tss_file = config['ref']['genome_tss']
        gene_file = config['ref']['genes']
        
        if biosamples_config.loc[sample, "alt_TSS"]:
            tss_file = biosamples_config.loc[sample, 'alt_TSS']
        if biosamples_config.loc[sample, "alt_genes"]:
            gene_file = biosamples_config.loc[sample, 'alt_genes']
        
        TSS_files.append(tss_file)
        gene_files.append(gene_file)
    
    biosamples_config["TSS"] = TSS_files
    biosamples_config["genes"] = gene_files
