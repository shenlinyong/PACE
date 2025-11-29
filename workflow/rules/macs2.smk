"""
MACS2 Peak Calling Rules for PACE
"""

rule call_macs_peaks:
    """Call peaks using MACS2"""
    input:
        accessibility = get_accessibility_files,
    params:
        pval = config['params_macs']['pval'],
        genome_size = config['params_macs']['genome_size'],
    conda:
        "../envs/pace-env.yml"
    output:
        narrowPeak = os.path.join(RESULTS_DIR, "{biosample}", "Peaks", "macs2_peaks.narrowPeak")
    resources:
        mem_mb = determine_mem_mb
    log:
        os.path.join(RESULTS_DIR, "{biosample}", "logs", "macs2.log")
    shell:
        """
        # Determine MACS2 format based on file type
        if [[ "{input.accessibility}" == *tagAlign* ]]; then
            FORMAT="BED"
        else
            FORMAT="AUTO"
        fi

        macs2 callpeak \
            -f $FORMAT \
            -g {params.genome_size} \
            -p {params.pval} \
            -n macs2 \
            --shift -75 \
            --extsize 150 \
            --nomodel \
            --keep-dup all \
            --call-summits \
            --outdir {RESULTS_DIR}/{wildcards.biosample}/Peaks \
            -t {input.accessibility} \
            2> {log}
        """


rule generate_chrom_sizes_bed_file:
    """Generate chromosome sizes BED file"""
    input:
        chrom_sizes = config['ref']['chrom_sizes']
    output:
        chrom_sizes_bed = os.path.join(RESULTS_DIR, "tmp", os.path.basename(config['ref']['chrom_sizes']) + '.bed')
    resources:
        mem_mb = determine_mem_mb
    shell:
        """
        mkdir -p $(dirname {output.chrom_sizes_bed})
        awk 'BEGIN {{OFS="\\t"}} {{if (NF > 0) print $1,"0",$2 ; else print $0}}' {input.chrom_sizes} > {output.chrom_sizes_bed}
        """


rule sort_narrowpeaks:
    """Sort narrowPeak file"""
    input:
        narrowPeak = os.path.join(RESULTS_DIR, "{biosample}", "Peaks", "macs2_peaks.narrowPeak"),
        chrom_sizes_bed = os.path.join(RESULTS_DIR, "tmp", os.path.basename(config['ref']['chrom_sizes']) + '.bed')
    params:
        chrom_sizes = config['ref']['chrom_sizes']
    conda:
        "../envs/pace-env.yml"
    output:
        narrowPeakSorted = os.path.join(RESULTS_DIR, "{biosample}", "Peaks", "macs2_peaks.narrowPeak.sorted")
    resources:
        mem_mb = determine_mem_mb
    shell:
        """
        # Intersect with chromosome sizes to remove non-standard chromosomes
        bedtools intersect -u -a {input.narrowPeak} -b {input.chrom_sizes_bed} | \
        bedtools sort -faidx {params.chrom_sizes} -i stdin > {output.narrowPeakSorted}
        """
