"""
Candidate Enhancer Region Generation Rules for PACE
Author: Linyong Shen @ Northwest A&F University
"""

rule make_candidate_regions:
    """Generate candidate enhancer regions from MACS2 peaks"""
    input:
        narrowPeak = os.path.join(RESULTS_DIR, "{biosample}", "Peaks", "macs2_peaks.narrowPeak.sorted"),
    params:
        tss_regions = config['ref']['genome_tss'],
        chrom_sizes = config['ref']['chrom_sizes'],
        regions_blocklist = config['ref']['regions_blocklist'],
        peakExtendFromSummit = config['params_candidate']['peakExtendFromSummit'],
        nStrongestPeaks = config['params_candidate']['nStrongestPeaks'],
        scripts_dir = SCRIPTS_DIR,
    conda:
        "../envs/pace-env.yml"
    output:
        candidateRegions = os.path.join(RESULTS_DIR, "{biosample}", "Peaks", "macs2_peaks.narrowPeak.sorted.candidateRegions.bed")
    resources:
        mem_mb = determine_mem_mb
    log:
        os.path.join(RESULTS_DIR, "{biosample}", "logs", "candidate_regions.log")
    shell:
        """
        python {params.scripts_dir}/pace_candidate_regions.py \
            --narrowPeak {input.narrowPeak} \
            --chrom_sizes {params.chrom_sizes} \
            --output {output.candidateRegions} \
            --blacklist {params.regions_blocklist} \
            --tss_regions {params.tss_regions} \
            --peakExtendFromSummit {params.peakExtendFromSummit} \
            --nStrongestPeaks {params.nStrongestPeaks} \
            2> {log}
        """
