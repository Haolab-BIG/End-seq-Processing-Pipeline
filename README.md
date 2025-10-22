# End-seq-Processing-Pipeline
This pipeline provides a fully containerized Singularity environment that bundles all required tools and dependencies. With a single command, the entire End-seq workflow—from raw FASTQ input through trimming, quality control, genome alignment, peak calling, fingerprinting, and PCA—can be executed reproducibly on any compatible system.

# Part I Workflow

Here stands an throughout workflow of End-seq data analysis.
<img width="2010" height="673" alt="End-seq" src="https://github.com/user-attachments/assets/d61be29f-fbf7-4f61-8e4c-1dd5d0bf8810" />

# Part II Requirements
1.  **Recommended Specs**:

      * 8-core CPU
      * 24 GB RAM

2.  **Singularity**: Must be installed on your system. Below are the detailed steps for installing on an Ubuntu 22.0.4 system. For other operating systems, please refer to the official installation guide: [https://docs.sylabs.io/guides/3.0/user-guide/installation.html](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)

      * **Step 1: Install System Dependencies**

        ```bash
        # Update package lists and install dependencies
        sudo apt-get update
        sudo apt-get install -y \
            build-essential \
            libseccomp-dev \
			libfuse3-dev \
            pkg-config \
            squashfs-tools \
            cryptsetup \
            curl wget git
        ```

      * **Step 2: Install Go Language**

        ```bash
        # Download and install Go
        wget https://go.dev/dl/go1.21.3.linux-amd64.tar.gz
        sudo tar -C /usr/local -xzvf go1.21.3.linux-amd64.tar.gz
        rm go1.21.3.linux-amd64.tar.gz

        # Configure Go environment variables and apply them
        echo 'export GOPATH=${HOME}/go' >> ~/.bashrc
        echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc
        source ~/.bashrc
        ```

      * **Step 3: Download, Build, and Install Singularity**

        ```bash
        # Note: The script navigates to /mnt/share/software. 
        # You can change this to your preferred directory for source code.
        cd /mnt/share/software

        # Download the Singularity CE source code
        wget https://github.com/sylabs/singularity/releases/download/v4.0.1/singularity-ce-4.0.1.tar.gz

        # Extract the archive and clean up
        tar -xvzf singularity-ce-4.0.1.tar.gz
        rm singularity-ce-4.0.1.tar.gz
        cd singularity-ce-4.0.1

        # Configure the build
        ./mconfig

        # Build Singularity (this can be time-consuming)
        cd builddir
        make

        # Install Singularity to the system
        sudo make install
        ```

      * **Step 4: Verify the Installation**

        ```bash
        # Check the installed version
        singularity --version

        # Display help information
        singularity -h
        ```

3.  **snakemake**: Snakemake must be installed on your system and requires a Python 3 distribution.

      ```bash
      pip install snakemake
      ```

4.  **Reference Data**: A directory containing bowtie index (Below are the detailed steps for the human hg38 genome. For other reference genomes, please download the corresponding files and replace them as needed).
      ```bash
      mkdir basement_data
      cd basement_data
      # Download Genome FASTA
      wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz
      # Unzip the files
      gunzip GRCh38.primary_assembly.genome.fa.gz
      gunzip gencode.v46.primary_assembly.annotation.gtf.gz
      # Remove scafford
      awk '/^>/ {p=0} /^>chr[0-9XYM]/ {p=1} p' GRCh38.primary_assembly.genome.fa > GRCh38.primary_assembly.genome.chr.fa
      # Build index
      mkdir hg38_chr_bowtie1_index
      singularity exec --cleanenv EndSeq.sif bowtie-build --threads 8 -f GRCh38.primary_assembly.genome.chr.fa ./hg38_chr_bowtie1_index/hg38_chr
      # Remove unnecessary files
      rm GRCh38.primary_assembly.genome.chr.fa
      rm GRCh38.primary_assembly.genome.fa
      ```

5.  **Required Files**:

      ```bash
      project_directory/
      ├── Scripts
            ├── config.yaml
            └── EndSeq.smk
      ├── Containers/
            └── EndSeq.sif
      ├── References/
            ├── illumina_adapter.fa
            └── hg38_chr_bowtie1_index
                  ├── hg38_chr.1.ebwt
                  ├── hg38_chr.2.ebwt
                  ├── hg38_chr.3.ebwt
                  ├── hg38_chr.4.ebwt
                  ├── hg38_chr.rev.1.ebwt
                  └── hg38_chr.rev.2.ebwt
      ```
      
      - **EndSeq.smk** — The main Snakemake workflow script.  
      - **config.yaml** — Configuration file containing paths, parameters, and sample information.  
        ⚠️ Must be located in the same directory as `EndSeq.smk`.
      - **EndSeq.sif** — Singularity container image with all required software and dependencies pre-installed.
      - **illumina_adapter.fa** — FASTA file containing Illumina adapter sequences; replace with your own if needed. 
      - **hg38_chr_bowtie1_index** — Reference genome index for Bowtie; replace with your preferred reference.

# Part III Running

   * **Example code**

      * **Step 1: Edit `config.yaml`**

        ```bash
        treatment:
          fq: "/project_directory/rawdata/Treatment.5w.fastq.gz"
          prefix: "Treatment"

        control:
          fq: "/project_directory/rawdata/Control.5w.fastq.gz"
          prefix: "Control"

        outputdir: "/project_directory/result"
        referencedir: "/project_directory/References/hg38_chr_bowtie1_index/hg38_chr"
        adapterFa: "/project_directory/References/illumina_adapter.fa"
        sif: "/project_directory/Containers/EndSeq.sif"
        threads: 8
        binSize: 10
        g: "hs"
        ```

      * **Step 2: run snakemake**

        ```bash
        snakemake -s pipeline.smk --cores 20 --use-singularity --singularity-args "--bind /project_directory:/project_directory"
        ```

   * **Command Parameters**

      - `treatment`:  Path to the treatment FASTQ file. For paired-end data, it is recommended to provide the path to the R1 file, gzipped, with User-defined prefix for all treatment output files (required)
      - `control`:    Path to the input control fastq. For paired-end data, it is recommended to provide the path to the R1 file, gzipped, with User-defined prefix for all treatment output files (optinal)
      - `outputdir`:    Path to the directory where the output will be stored (required)
      - `referencedir`: Path to the directory where bowtie reference build with prefix (required)
      - `adapterFa`:    Path to the adapter fasta (required)
      - `sif`:          Path to the singularity environment file (required)
      - `threads`:      Number of threads to use (optional, default: 8)
      - `binSize`:      Number of binsize to use (optional, default: 10)
      - `g`:            specise from macs3: hs (human); mm (mouse); ce (C. elegans); dm (Drosophila melanogaster); ...

# Part IV Output

   * **Output Structure**
      ```bash
      result/
      ├── bam
            ├── Control.bowtie.stats
            ├── Control.DeDup.bam
            ├── Control.DeDup.bam.bai
            ├── Control.flagstat.txt
            ├── Control.markdup.log
            ├── Treatment.bowtie.stats
            ├── Treatment.DeDup.bam
            ├── Treatment.DeDup.bam.bai
            ├── Treatment.flagstat.txt
            └── Treatment.markdup.log
      ├── bw/
            ├── Control.DeDup.bw
            └── Treatment.DeDup.bw
      ├── figure/
            ├── BW_compare_PCA.pdf
            ├── BW_compare_cor.pdf
            ├── fingerprints.pdf
            └── Treatment.peak.pdf
      ├── multiqc/
            ├── multiqc_data/
            └── multiqc_report.html
      └── peak/
            ├── Treatment.macs3.stats
            ├── Treatment.macs3.stdout
            ├── Treatment_peaks.broadPeak
            ├── Treatment_peaks.gappedPeak
            └── Treatment_peaks.xls
      ├── rawdata.qc/
            ├── Control.5w_fastqc.html
            ├── Control.5w_fastqc.zip
            ├── Treatment.5w_fastqc.html
            └── Treatment.5w_fastqc.zip
      ```
      
   * **Output Interpretation**

      - **`*.bowtie.stats`**

        - **Content**: Contains Bowtie alignment summary statistics, including the total number of reads processed, reads aligned, reads discarded, and uniquely mapped reads. It provides an overview of mapping quality and efficiency for each FASTQ file.
        - **Application**: Used to assess alignment quality and sequencing library performance. These statistics help in troubleshooting mapping issues, evaluating experiment success, and can be parsed by downstream tools like MultiQC for visualization and comparison across samples.
          
          <img width="490" height="105" alt="图片" src="https://github.com/user-attachments/assets/ad65b4e6-d210-4af1-8dbe-174f49c3a75e" />

      - **`*.DeDup.bam`**

        - **Content**: This is the main alignment file in Binary Alignment Map (BAM) format. It contains all the sequencing reads and their mapping coordinates on the reference genome. This version has had duplicate reads (PCR duplicates) removed. For more information please refer to: https://genome.ucsc.edu/goldenpath/help/bam.html.
        - **Application**: It's the primary evidence for read alignment and can be used for detailed inspection in genome browsers or for downstream analyses.

      - **`*.dedup.bam.bai`**

        - **Content**: This is the index file for the BAM file.
        - **Application**: It allows for fast random access to the BAM file, which is essential for visualization software (like IGV) to quickly load and display alignments for a specific genomic region without reading the entire file.

      - **`*.bw`**

        - **Content**: A BigWig file that represents the End-seq signal coverage across the genome. It shows the read density (how many reads cover each position) in a compressed format. For more information please refer to: https://genome.ucsc.edu/goldenpath/help/bigWig.html
        - **Application**: Primarily used for visualization. You can load this file into a genome browser (e.g., IGV, UCSC Genome Browser) to see a "signal track" that shows gene expression levels visually across chromosomes.

      - **`*.dedup.bam.bai`**

        - **Content**: This is the index file for the BAM file.
        - **Application**: It allows for fast random access to the BAM file, which is essential for visualization software (like IGV) to quickly load and display alignments for a specific genomic region without reading the entire file.

      - **`*.flagstat.txt`**

        - **Content**: Contains alignment statistics generated by samtools flagstat after removing duplicate reads. It reports the total number of reads, mapped reads, properly paired reads, singletons, and the number of duplicate reads removed, providing a summary of the final, deduplicated BAM file.
        - **Application**: Used to evaluate the quality of the deduplicated alignment, check library complexity, and ensure that downstream analyses (e.g., peak calling, coverage calculation) are based on high-quality, non-redundant reads.

          <img width="425" height="323" alt="图片" src="https://github.com/user-attachments/assets/dd07e6e2-9a7a-4f93-ae59-c4f1ab79726b" />

      - **`*.markdup.log`**

        - **Content**: Log file generated by `samtools markdup`, summarizing read duplication. It includes READ (total number of input reads), WRITTEN (reads retained after removing duplicates), EXCLUDED, EXAMINED, counts of PAIRED and SINGLE reads, as well as DUPLICATE SINGLE/PAIR and DUPLICATE TOTAL.
        - **Application**: Used to evaluate library complexity and duplication rate. A high WRITTEN/READ ratio indicates low duplication and good library complexity, while a low ratio suggests high PCR duplication or low-complexity sequencing.
        
		  <img width="380" height="105" alt="图片" src="https://github.com/user-attachments/assets/8aa4b3f7-ec00-4ab2-af1c-a2a7572d0386" />

      - **`multiqc_report`** : Open multiqc_report.html in a web browser to explore all sections interactively.

        - **General Statistics**: A combined table summarizing important metrics for each sample:
	  
          <img width="1632" height="229" alt="图片" src="https://github.com/user-attachments/assets/54a5e5fe-7e1f-45fc-b377-050426325d81" />

        - **FastQC**: Quality-control metrics on raw and trimmed reads, including 'Sequence Counts', 'Sequence Quality Histograms', 'Per Sequence Quality Scores', 'Per Base Sequence Content', 'Per Sequence GC Content', 'Per Base N Content', 'Sequence Length Distribution', 'Sequence Duplication Levels', 'Overrepresented sequences by sample', 'Top overrepresented sequences', 'Adapter Content':
	      
          - Sequence Quality Histograms: The mean quality value across each base position in the read.
	  
            <img width="1628" height="617" alt="图片" src="https://github.com/user-attachments/assets/d2911bd4-a13f-4d3f-8cce-c3e4bab3858f" />

          - Adapter Content: The cumulative percentage count of the proportion of your library which has seen each of the adapter sequences at each position.
	  
            <img width="1625" height="550" alt="图片" src="https://github.com/user-attachments/assets/1d3bb888-3bb0-424c-bf8d-b7323c54e960" />

        - **Samtools**: This module parses the output from samtools flagstat to report the percentage of total, mapped, and properly paired reads, providing a summary of alignment quality. Helps evaluate the effectiveness of deduplication and ensures that downstream analyses (e.g., peak calling, coverage profiling) are based on unique, non-redundant reads.
	  
          <img width="1636" height="722" alt="图片" src="https://github.com/user-attachments/assets/38efbb17-2596-46f3-a16a-637c885740ed" />

        - **Bowtie**: Alignment statistics such as total reads, uniquely mapped reads, and multi-mapping rates:
	  
          <img width="1630" height="578" alt="图片" src="https://github.com/user-attachments/assets/fb185ccb-58c3-4062-8a68-e9970df1046b" />

      - **`*BW_compare_PCA.pdf`**

        - **Content**: PDF file showing the principal component analysis (PCA) of BigWig signal profiles across multiple samples. It visualizes sample-to-sample similarity and variance based on genome-wide coverage or signal intensities.
        - **Application**: Used to assess the overall relationship between samples, detect outliers, and evaluate batch effects or experimental reproducibility in End-seq or ChIP-seq datasets.
          <img width="709" height="701" alt="PCA" src="https://github.com/user-attachments/assets/99c7ede0-6cd1-4fb5-8933-42a0b5612460" />
        
      - **`*BW_compare_cor.pdf`**

        - **Content**: PDF file showing a heatmap of pairwise correlations between samples based on BigWig signal profiles. It typically includes correlation values (Pearson) and visually represents sample similarity across the genome.
        - **Application**: Used to assess consistency and reproducibility between samples, identify outliers, and evaluate experimental quality in End-seq.

          <img width="860" height="766" alt="图片" src="https://github.com/user-attachments/assets/fa058e75-5d40-46df-a1fa-a2e6de583a21" />

      - **`*fingerprints.pdf`**

        - **Content**: PDF file generated by `plotFingerprint` showing the cumulative read coverage across the genome for each BAM file. It visualizes enrichment patterns and sequencing depth consistency among samples.
        - **Application**: In the fingerprint plot, a larger separation between treatment and control curves, indicates stronger enrichment and higher signal-to-noise ratio. Also, the fingerprint plot can help decide whether to call narrow peaks or broad peaks: Narrow peaks are appropriate when the signal is sharp and localized; Broad peaks are used when the signal spans wide genomic regions with diffuse enrichment, such as histone modifications.

          <img width="830" height="619" alt="图片" src="https://github.com/user-attachments/assets/8a4bc5cd-d2a8-4eff-974e-d08e86756f57" />

      - **`*..macs3.stats`**

        - **Content**: Contains summary statistics from MACS3 peak calling, including number of input reads, effective genome size, estimated fragment size, number of peaks called, and other runtime information.
        - **Application**: Used to check if MACS3 ran successfully and to detect any errors or warnings during the peak calling process.

      - **`*_peaks.broadPeak`**

        - **Content**: BED6+3 format file (similar to narrowPeak, but without the 10th column for peak summits). Only available when `--broad` is enabled. In broad peak mode, the peak summit isn’t called, so the 5th, 7th–9th columns are the mean values across the peak region. Can be loaded directly into UCSC Genome Browser with `--trackline`.
        - **Columns**:
        | Column | Description |
        |--------|-------------|
        | chrom  | Chromosome name |
        | start  | Start position of the broad peak (0-based) |
        | end    | End position of the broad peak (not inclusive) |
        | name   | Peak name or ID |
        | score  | Mean score across the broad peak (similar to narrowPeak 5th column) |
        | strand | Strand information (‘+’, ‘-’, or ‘.’ if not applicable) |
        | signalValue | Mean enrichment signal across the peak |
        | pValue | Mean -log10 p-value across the peak |
        | qValue | Mean -log10 q-value (FDR) across the peak |

        - **Application**: Used for visualization DNA breadk signals, with bigwig file.

          <img width="1134" height="350" alt="图片" src="https://github.com/user-attachments/assets/778c824c-8285-4278-91d1-cd165415d3ad" />

      - **`*_peaks.gappedPeak`**

        - **Content**: BED12+3 format file containing broad regions and narrow peaks within them. Only available when `--broad` is enabled. Can be loaded into UCSC Genome Browser. Columns 5, 7–9 may need adjustment if integrating with narrowPeak conventions.
        - **Columns**:
        | Column | Description |
        |--------|-------------|
        | chrom       | Chromosome name |
        | start       | Start of the broad region (0-based) |
        | end         | End of the broad region (not inclusive) |
        | name        | Peak name or ID |
        | score       | Score for display in UCSC browser (grey levels, similar to narrowPeak 5th column) |
        | strand      | Strand information (‘+’, ‘-’, or ‘.’) |
        | thickStart  | Start of the first narrow peak within the broad region |
        | thickEnd    | End of the first narrow peak within the broad region |
        | itemRgb     | RGB color for UCSC browser (0 uses default color) |
        | blockCount  | Number of blocks (including 1bp at start and end of broad regions) |
        | blockSizes  | Comma-separated lengths of each block |
        | blockStarts | Comma-separated start positions of each block relative to `start` |
        | foldChange  | Fold-change of enrichment within the peak |
        | -log10(pvalue) | -log10 p-value for the peak |
        | -log10(qvalue) | -log10 q-value (FDR) for the peak |

        - **Application**: Used to analyze subpeak structure, study internal peak features, or visualize complex enrichment patterns in broad regions.

      - **`*_peaks.xls`**

        - **Content**: Tab-delimited summary of all peaks called by MACS3 with detailed metrics.
        - **Columns**:
        | Column | Description |
        |--------|-------------|
        | chr    | Chromosome name |
        | start  | Peak start position (0-based) |
        | end    | Peak end position (not inclusive) |
        | length | Peak length (end - start) |
        | pileup | Maximum pileup (number of overlapping tags) at the peak |
        | -log10(pvalue) | -log10 of p-value for peak significance |
        | fold_enrichment | Fold enrichment of the peak over background |
        | -log10(qvalue) | -log10 of q-value (FDR) for peak significance |
        | name   | Peak name or ID |

        - **Application**: Used for quantitative peak analysis, filtering peaks by significance, fold enrichment, or integrating with downstream functional annotation, motif analysis, and visualization.

      - **`*.peak.pdf`**

        - **Content**: Heatmap visualizing read enrichment over peaks. Generated using `plotHeatmap` from deepTools with a `*_peaks.broadPeak` and `*bw` input. The heatmap shows signal intensity (color-coded, viridis colormap) across all peaks, with missing data represented in white. The height and width of the heatmap are set for clear visualization of peak patterns.
        - **Application**: Used to assess global enrichment patterns across peaks. Peaks with strong enrichment appear as high-intensity bands or curves; if the signal is higher than control samples, it indicates that peak calling was successful and represents true biological enrichment.
       
          <img width="495" height="746" alt="图片" src="https://github.com/user-attachments/assets/3859267d-c74a-4d3e-be3c-9c6c4f8a61e0" />

# Part V Video Tutorials

   Watch this video tutorial to see a complete walkthrough of running the pipeline:
   
   https://github.com/user-attachments/assets/30db9693-b61c-439c-bd5a-854ae7479763

# Reference
[1] Wong, Nancy et al. “END-seq: An Unbiased, High-Resolution, and Genome-Wide Approach to Map DNA Double-Strand Breaks and Resection in Human Cells.” Methods in molecular biology (Clifton, N.J.) vol. 2153 (2021): 9-31. doi:10.1007/978-1-0716-0644-5_2





