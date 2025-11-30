"""
Neighborhood Analysis Rules for PACE
Author: Linyong Shen @ Northwest A&F University
"""

rule create_neighborhoods:
    """Create neighborhood analysis with multi-omics support"""
    input:
        candidateRegions = os.path.join(RESULTS_DIR, "{biosample}", "Peaks", "macs2_peaks.narrowPeak.sorted.candidateRegions.bed"),
    params:
        ATAC = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, "ATAC"] or '',
        DHS = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, "DHS"] if "DHS" in BIOSAMPLES_CONFIG.columns else '',
        default = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, 'default_accessibility_feature'],
        H3K27ac = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, "H3K27ac"] if "H3K27ac" in BIOSAMPLES_CONFIG.columns and pd.notna(BIOSAMPLES_CONFIG.loc[wildcards.biosample, "H3K27ac"]) else '',
        H3K4me1 = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, "H3K4me1"] if "H3K4me1" in BIOSAMPLES_CONFIG.columns and pd.notna(BIOSAMPLES_CONFIG.loc[wildcards.biosample, "H3K4me1"]) else '',
        H3K4me3 = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, "H3K4me3"] if "H3K4me3" in BIOSAMPLES_CONFIG.columns and pd.notna(BIOSAMPLES_CONFIG.loc[wildcards.biosample, "H3K4me3"]) else '',
        methylation = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, "methylation"] if "methylation" in BIOSAMPLES_CONFIG.columns and pd.notna(BIOSAMPLES_CONFIG.loc[wildcards.biosample, "methylation"]) else '',
        expression = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, "RNA_seq"] if "RNA_seq" in BIOSAMPLES_CONFIG.columns and pd.notna(BIOSAMPLES_CONFIG.loc[wildcards.biosample, "RNA_seq"]) else '',
        genes = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, 'genes'] if 'genes' in BIOSAMPLES_CONFIG.columns else config['ref']['genes'],
        chrom_sizes = config['ref']['chrom_sizes'],
        activity_method = config.get('activity_method', 'geometric_mean'),
        scripts_dir = SCRIPTS_DIR
    conda:
        "../envs/pace-env.yml"
    output:
        enhList = os.path.join(RESULTS_DIR, "{biosample}", "Neighborhoods", "EnhancerList.txt"),
        geneList = os.path.join(RESULTS_DIR, "{biosample}", "Neighborhoods", "GeneList.txt"),
        neighborhoodDirectory = directory(os.path.join(RESULTS_DIR, "{biosample}", "Neighborhoods")),
    resources:
        mem_mb = 32*1000
    log:
        os.path.join(RESULTS_DIR, "{biosample}", "logs", "neighborhoods.log")
    shell:
        """
        # Build accessibility argument
        if [ -n "{params.ATAC}" ]; then
            ACC_FILE="{params.ATAC}"
            ACC_TYPE="ATAC"
        else
            ACC_FILE="{params.DHS}"
            ACC_TYPE="DHS"
        fi
        
        # Build optional arguments
        EXTRA_ARGS=""
        if [ -n "{params.H3K27ac}" ]; then
            EXTRA_ARGS="$EXTRA_ARGS --H3K27ac {params.H3K27ac}"
        fi
        if [ -n "{params.H3K4me1}" ]; then
            EXTRA_ARGS="$EXTRA_ARGS --H3K4me1 {params.H3K4me1}"
        fi
        if [ -n "{params.H3K4me3}" ]; then
            EXTRA_ARGS="$EXTRA_ARGS --H3K4me3 {params.H3K4me3}"
        fi
        if [ -n "{params.methylation}" ]; then
            EXTRA_ARGS="$EXTRA_ARGS --methylation {params.methylation}"
        fi
        if [ -n "{params.expression}" ]; then
            EXTRA_ARGS="$EXTRA_ARGS --expression {params.expression}"
        fi
        
        python {params.scripts_dir}/pace_neighborhoods.py \
            --candidate_regions {input.candidateRegions} \
            --genes {params.genes} \
            --chrom_sizes {params.chrom_sizes} \
            --output_dir {output.neighborhoodDirectory} \
            --accessibility_file "$ACC_FILE" \
            --accessibility_type "$ACC_TYPE" \
            --activity_method {params.activity_method} \
            $EXTRA_ARGS \
            2> {log}
        """
