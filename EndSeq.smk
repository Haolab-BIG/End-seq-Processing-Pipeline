import os
from pathlib import Path

# Load configuration file
configfile: "config.yaml"

# Define output directories using absolute paths (created by Snakemake automatically)
OUTPUT_DIR = Path(config["outputdir"]).resolve()
QC_DIR = OUTPUT_DIR / "rawdata.qc"
TRIM_DIR = OUTPUT_DIR / "bam"
BAM_DIR = OUTPUT_DIR / "bam"
BW_DIR = OUTPUT_DIR / "bw"
FIGURE_DIR = OUTPUT_DIR / "figure"
PEAK_DIR = OUTPUT_DIR / "peak"
MULTIQC_DIR = OUTPUT_DIR / "multiqc"

# Define samples based on user-defined prefixes in config
def get_samples_raw_names():
    samples = []
    treatment_fq = config["treatment"]["fq"]
    treatment_basename = re.sub(r'\.f(ast)?q\.gz$', '', os.path.basename(treatment_fq))
    samples.append(treatment_basename)
    if "control" in config and config["control"]["fq"]:
        control_fq = config["control"]["fq"]
        control_basename = re.sub(r'\.f(ast)?q\.gz$', '', os.path.basename(control_fq))
        samples.append(control_basename)
    return samples
SAMPLES_RAW_NAMES = get_samples_raw_names() 
#
SAMPLES = [config["treatment"]["prefix"]]
if "control" in config and config["control"]["fq"]:
    SAMPLES.append(config["control"]["prefix"])


# Final target rule: defines all expected output files
rule all:
    input:
        # MultiQC report
        MULTIQC_DIR / "multiqc_report.html",
        # Peak calling result
        PEAK_DIR / f"{config['treatment']['prefix']}_peaks.broadPeak",
        # Quality control figures
        [QC_DIR / f"{s}_fastqc.html" for s in SAMPLES_RAW_NAMES],
        FIGURE_DIR / "BW_compare_PCA.pdf",
        FIGURE_DIR / "BW_compare_cor.pdf",
        FIGURE_DIR / "fingerprints.pdf",
        FIGURE_DIR / f"{config['treatment']['prefix']}.peak.pdf",
        # Deduplicated BAM files and their indexes
        [BAM_DIR / f"{s}.DeDup.bam" for s in SAMPLES],
        [BAM_DIR / f"{s}.DeDup.bam.bai" for s in SAMPLES],
        # BigWig files
        [BW_DIR / f"{s}.DeDup.bw" for s in SAMPLES]


# Rule 2: Quality control for raw gzipped fastq using FastQC
rule fastqc_raw:
    input:
        lambda wildcards: config["treatment"]["fq"] if wildcards.sample == SAMPLES_RAW_NAMES[0] else config["control"]["fq"]
    output:
        html = QC_DIR / "{sample}_fastqc.html",
        zip = QC_DIR / "{sample}_fastqc.zip"
    container: config["sif"]
    params:
        threads = config["threads"]
    shell:
        """
        fastqc -o {QC_DIR} -t {params.threads} -q {input}
        """

# Rule 3: Adapter trimming using cutadapt
rule trim_adapters:
    input:
        lambda wildcards: config["treatment"]["fq"] if wildcards.sample == config["treatment"]["prefix"] else config["control"]["fq"]
    output:
        temp(TRIM_DIR / "{sample}.trimmed.fastq.gz")
    container: config["sif"]
    params:
        threads = config["threads"],
        adapter_fa = config["adapterFa"]  # Absolute path from config
    shell:
        """
        cutadapt -q 25 -m 36 -e 0.1 -j {params.threads} \
            -a GGGGGGGG -a AAAAAAAA -g GGGGGGGG -g AAAAAAAA \
            -a file:{params.adapter_fa} --times 13 \
            -o {output} {input}
        """


# Rule 4: Quality control for trimmed fastq
rule fastqc_trimmed:
    input:
        TRIM_DIR / "{sample}.trimmed.fastq.gz"
    output:
        html = temp(TRIM_DIR / "{sample}.trimmed_fastqc.html"),
        zip = temp(TRIM_DIR / "{sample}.trimmed_fastqc.zip")
    container: config["sif"]
    shell:
        """
        fastqc -o {TRIM_DIR} --noextract -f fastq {input}
        """


# Rule 5: Align trimmed reads to reference genome using bowtie
rule bowtie_align:
    input:
        TRIM_DIR / "{sample}.trimmed.fastq.gz"
    output:
        bam = BAM_DIR / "{sample}.sort.bam",
        bowtie_stats = BAM_DIR / "{sample}.bowtie.stats"
    container: config["sif"]
    params:
        threads = config["threads"],
        ref_dir = config["referencedir"]  # Absolute path from config
    shell:
        """
        bowtie -S -p {params.threads} -x {params.ref_dir} \
            {input} --best --strata --all -m 1 -n 3 -l 50 \
            2> {output.bowtie_stats} | \
            samtools view -q 255 -bS - | \
            samtools sort -n -@ {params.threads} - -T {BAM_DIR}/{wildcards.sample} \
            -O BAM -o {output.bam}
        """


# Rule 6: Remove duplicate reads using samtools
rule remove_duplicates:
    input:
        BAM_DIR / "{sample}.sort.bam"
    output:
        dedup_bam = BAM_DIR / "{sample}.DeDup.bam",
        dedup_bai = BAM_DIR / "{sample}.DeDup.bam.bai",
        flagstat = BAM_DIR / "{sample}.flagstat.txt"
    container: config["sif"]
    params:
        threads = config["threads"]
    shell:
        """
        # Fixmate and sort reads before marking duplicates
        samtools fixmate -r -m {input} - | \
            samtools sort -@ {params.threads} - | \
            samtools markdup -r -s - {output.dedup_bam} \
            2> {BAM_DIR}/{wildcards.sample}.markdup.log
        
        # Generate flagstat report
        samtools flagstat {output.dedup_bam} > {output.flagstat}
        
        # Index the deduplicated BAM file
        samtools index {output.dedup_bam}
        
        # Remove intermediate sorted BAM to save space
        rm {input}
        """


# Rule 7: Convert BAM to BigWig for visualization
rule bam_to_bw:
    input:
        bam = BAM_DIR / "{sample}.DeDup.bam",
        bai = BAM_DIR / "{sample}.DeDup.bam.bai"
    output:
        BW_DIR / "{sample}.DeDup.bw"
    container: config["sif"]
    params:
        threads = config["threads"],
        bin_size = config["binSize"]
    shell:
        """
        bamCoverage -b {input.bam} -o {output} \
            --ignoreForNormalization chrM -p {params.threads} \
            --binSize {params.bin_size} --normalizeUsing CPM
        """


# Rule 8: Calculate BigWig correlations using deepTools
rule multi_bigwig_summary:
    input:
        [BW_DIR / f"{s}.DeDup.bw" for s in SAMPLES]
    output:
        temp(FIGURE_DIR / "BW_compare.npz")
    container: config["sif"]
    params:
        threads = config["threads"],
        labels = " ".join(SAMPLES)  # Use user-defined prefixes as labels
    shell:
        """
        multiBigwigSummary bins -b {input} --labels {params.labels} \
            -out {output} -p {params.threads}
        """


# Rule 9: Plot PCA and correlation heatmap from BigWig data
rule plot_pca_cor:
    input:
        FIGURE_DIR / "BW_compare.npz"
    output:
        pca = FIGURE_DIR / "BW_compare_PCA.pdf",
        cor = FIGURE_DIR / "BW_compare_cor.pdf"
    container: config["sif"]
    shell:
        """
        plotPCA -in {input} -o {output.pca}
        plotCorrelation -in {input} --corMethod pearson \
            --skipZeros -o {output.cor} --whatToPlot heatmap \
            --colorMap RdYlBu --plotNumbers
        """


# Rule 10: Plot fingerprint for alignment quality assessment
rule plot_fingerprint:
    input:
        [BAM_DIR / f"{s}.DeDup.bam" for s in SAMPLES]
    output:
        FIGURE_DIR / "fingerprints.pdf"
    container: config["sif"]
    params:
        labels = " ".join(SAMPLES)  # Use user-defined prefixes as labels
    shell:
        """
        plotFingerprint -b {input} --labels {params.labels} \
            --skipZeros --plotFile {output}
        """


# Rule 11: Call broad peaks using MACS3
rule call_peaks:
    input:
        treat_bam = BAM_DIR / f"{config['treatment']['prefix']}.DeDup.bam",
        control_bam = BAM_DIR / f"{config['control']['prefix']}.DeDup.bam" if ("control" in config and config["control"]["fq"]) else None
    output:
        broad_peak = PEAK_DIR / f"{config['treatment']['prefix']}_peaks.broadPeak",
        stats = PEAK_DIR / f"{config['treatment']['prefix']}.macs3.stats",
        stdout = PEAK_DIR / f"{config['treatment']['prefix']}.macs3.stdout"
    container: config["sif"]
    params:
        genome = config["g"],
        treatment_prefix = config["treatment"]["prefix"],
        peak_dir = PEAK_DIR,
        bowtie_dir = BAM_DIR
    shell:
        """
        if [ -z "{input.control_bam}" ]; then
            macs3 callpeak -t {input.treat_bam} \
                -f BAM -g {params.genome} --outdir {params.peak_dir} \
                -n {params.treatment_prefix} --broad --nomodel --nolambda \
                --llocal 100000 --keep-dup all -q 0.01 \
                2> >(tee -a {output.stats} >&2) \
                1> {output.stdout}
        else
            macs3 callpeak -t {input.treat_bam} -c {input.control_bam} \
                -f BAM -g {params.genome} --outdir {params.peak_dir} \
                -n {params.treatment_prefix} --broad --nomodel --nolambda \
                --keep-dup all -q 0.01 \
                2> >(tee -a {output.stats} >&2) \
                1> {output.stdout}
        fi
        """

# Rule 12: Compute matrix for peak enrichment analysis
rule compute_matrix:
    input:
        peaks = PEAK_DIR / f"{config['treatment']['prefix']}_peaks.broadPeak",
        bw_files = [BW_DIR / f"{s}.DeDup.bw" for s in SAMPLES]
    output:
        temp(FIGURE_DIR / f"{config['treatment']['prefix']}.peak.gz")
    container: config["sif"]
    params:
        threads = config["threads"]
    shell:
        """
        computeMatrix scale-regions -p {params.threads} \
            -S {input.bw_files} -R {input.peaks} \
            -a 5000 -b 5000 --regionBodyLength 1000 \
            --missingDataAsZero -o {output}
        """


# Rule 13: Plot heatmap for peak enrichment
rule plot_peak_heatmap:
    input:
        FIGURE_DIR / f"{config['treatment']['prefix']}.peak.gz"
    output:
        FIGURE_DIR / f"{config['treatment']['prefix']}.peak.pdf"
    container: config["sif"]
    shell:
        """
        plotHeatmap -m {input} -out {output} \
            --colorMap viridis --missingDataColor white \
            --heatmapHeight 12 --heatmapWidth 4 \
            --startLabel "Start" --endLabel "End"
        """


# Rule 14: Generate MultiQC report summarizing all QC metrics
rule multiqc_report:
    input:
        # FastQC reports (raw and trimmed)
        # [QC_DIR / f"{s}_fastqc.html" for s in SAMPLES],
        [TRIM_DIR / f"{s}.trimmed_fastqc.html" for s in SAMPLES],
        # Alignment statistics
        [BAM_DIR / f"{s}.bowtie.stats" for s in SAMPLES],
        [BAM_DIR / f"{s}.flagstat.txt" for s in SAMPLES],
        # Peak calling statistics
        PEAK_DIR / f"{config['treatment']['prefix']}.macs3.stats"
    output:
        MULTIQC_DIR / "multiqc_report.html"
    container: config["sif"]
    shell:
        """
        multiqc {TRIM_DIR}/* {BAM_DIR}/*.bowtie.stats \
            {BAM_DIR}/*flagstat.txt -o {MULTIQC_DIR} --force
        """


