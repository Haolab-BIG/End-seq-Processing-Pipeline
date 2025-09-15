# End-seq-Processing-Pipeline
---
The End-seq analysis pipeline processes raw FASTQ data through a series of steps including adapter trimming, quality control, genome mapping, peak calling, fingerprint profiling, and principal component analysis.
# Part I Introduction
---
## i. Workflow
Here stands an throughout workflow of End-seq data analysis.
<img width="2010" height="673" alt="End-seq" src="https://github.com/user-attachments/assets/4688075e-b5c0-4f3f-b94a-3103a6fc384a" />

## ii. Features
This pipeline provides a fully containerized Singularity environment that bundles all required tools and dependencies. With a single command, the entire End-seq workflow—from raw FASTQ input through trimming, quality control, genome alignment, peak calling, fingerprinting, and PCA—can be executed reproducibly on any compatible system.

# Part II Requirements
---
1.  **Recommended System Configuration**:

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

3.  **Download Files**:

      * `run_Endseq.sh`
      * `End-seq.sif` (The Singularity container)
      * `illumina_adapter.fa`

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
      singularity exec --cleanenv End-seq.sif bowtie-build --threads 8 -f GRCh38.primary_assembly.genome.chr.fa ./hg38_chr_bowtie1_index/hg38_chr
      # Remove unnecessary files
      rm GRCh38.primary_assembly.genome.chr.fa
      rm GRCh38.primary_assembly.genome.fa
      ```
6.   **Required File Structure**
      ```bash
      basement_data/
      ├── End-seq.sif
      ├── hg38_chr_bowtie1_index/
            ├── hg38_chr.1.ebwt
            ├── hg38_chr.2.ebwt
            ├── hg38_chr.3.ebwt
            ├── hg38_chr.4.ebwt
            ├── hg38_chr.rev.1.ebwt
            └── hg38_chr.rev.2.ebwt
      ├── illumina_adapter.fa
      └── run_Endseq.sh
      ```

# Part III Running

   * **Example code**

      ```bash
      bash ./run_Endseq.sh --treatmentFq ./rawdata/Treatment.10w.fastq.gz \
                           --controlFq ./rawdata/Control.10w.fastq.gz \
                           --outputdir ./result \
                           --referencedir ./basement_data/hg38_chr_bowtie1_index/hg38_chr \
                           --adapterFa ./basement_data/illumina_adapter.fa \
                           --sif ./basement_data/End-seq.sif \
                           --threads 8 \
                           --binSize 10 \
                           --g hs
      ```
   * **Command Parameters**

      - `--treatmentFq`:  Path to the treatment fastq (required)
      - `--controlFq`:    Path to the input control fastq (optinal)
      - `--outputdir`:    Path to the directory where the output will be stored (required)
      - `--referencedir`: Path to the directory where bowtie reference build with prefix (required)
      - `--adapterFa`:    Path to the adapter fasta (required)
      - `--sif`:          Path to the singularity environment file (required)
      - `--threads`:      Number of threads to use (optional, default: 8)
      - `--binSize`:      Number of binsize to use (optional, default: 10)
      - `--g`:            specise from macs3: hs (human); mm (mouse); ce (C. elegans); dm (Drosophila melanogaster); ...

# Part IV Output

   * **Output Structure**
      ```bash
      result/
      ├── bam
            ├── Control.10w.bowtie.stats
            ├── Control.10w.DeDup.bam
            ├── Control.10w.DeDup.bam.bai
            ├── Control.10w.flagstat.txt
            ├── Control.10w.markdup.log
            ├── Treatment.10w.bowtie.stats
            ├── Treatment.10w.DeDup.bam
            ├── Treatment.10w.DeDup.bam.bai
            ├── Treatment.10w.flagstat.txt
            └── Treatment.10w.markdup.log
      ├── bw/
            ├── Control.10w.DeDup.bw
            └── Treatment.10w.DeDup.bw
      ├── figure/
            ├── BW_compare.pdf
            ├── BW_compare_cor.pdf
            ├── fingerprints.pdf
            └── Treatment.10w.peak.pdf
      ├── multiqc/
            ├── multiqc_data/
            └── multiqc_report.html
      └── peak/
            ├── Treatment.10w.macs3.stats
            ├── Treatment.10w_peaks.broadPeak
            ├── Treatment.10w_peaks.gappedPeak
            └── Treatment.10w_peaks.xls
      ```
   * **Output Interpretation**

      - **`*.bowtie.stats`**

        - **Content**: Contains Bowtie alignment summary statistics, including the total number of reads processed, reads aligned, reads discarded, and uniquely mapped reads. It provides an overview of mapping quality and efficiency for each FASTQ file.
        - **Application**: Used to assess alignment quality and sequencing library performance. These statistics help in troubleshooting mapping issues, evaluating experiment success, and can be parsed by downstream tools like MultiQC for visualization and comparison across samples.

      - **`*.DeDup.bam`**

        - **Content**: This is the main alignment file in Binary Alignment Map (BAM) format. It contains all the sequencing reads and their mapping coordinates on the reference genome. This version has had duplicate reads (PCR duplicates) removed. For more information please refer to: https://genome.ucsc.edu/goldenpath/help/bam.html.
        - **Application**: It's the primary evidence for read alignment and can be used for detailed inspection in genome browsers or for downstream analyses.

      - **`*.dedup.bam.bai`**

        - **Content**: This is the index file for the BAM file.
        - **Application**: It allows for fast random access to the BAM file, which is essential for visualization software (like IGV) to quickly load and display alignments for a specific genomic region without reading the entire file.

      - **`*.bw`**

        - **Content**: A BigWig file that represents the End-seq signal coverage across the genome. It shows the read density (how many reads cover each position) in a compressed format.
        - **Application**: Primarily used for visualization. You can load this file into a genome browser (e.g., IGV, UCSC Genome Browser) to see a "signal track" that shows gene expression levels visually across chromosomes.

        *\<p align="center"\> [IGV photo] \</p\>*

      - **`*.dedup.bam.bai`**

        - **Content**: This is the index file for the BAM file.
        - **Application**: It allows for fast random access to the BAM file, which is essential for visualization software (like IGV) to quickly load and display alignments for a specific genomic region without reading the entire file.

      - **`*.flagstat.txt`**

        - **Content**: Contains alignment statistics generated by samtools flagstat after removing duplicate reads. It reports the total number of reads, mapped reads, properly paired reads, singletons, and the number of duplicate reads removed, providing a summary of the final, deduplicated BAM file.
        - **Application**: Used to evaluate the quality of the deduplicated alignment, check library complexity, and ensure that downstream analyses (e.g., peak calling, coverage calculation) are based on high-quality, non-redundant reads.

      - **`*.markdup.log`**

        - **Content**: Log file generated by samtools markdup, summarizing read duplication information. It includes the total number of reads examined, reads written to the deduplicated BAM file, excluded reads, paired and single reads, as well as the number of duplicate reads detected and removed. This reflects the degree of PCR duplication in the library.
        - **Application**: Used to assess library complexity and duplication rate. Helps evaluate the effectiveness of deduplication and ensures that downstream analyses (e.g., peak calling, coverage profiling) are based on unique, non-redundant reads.

      - **`*BW_compare.pdf`**

        - **Content**: PDF file showing the principal component analysis (PCA) of BigWig signal profiles across multiple samples. It visualizes sample-to-sample similarity and variance based on genome-wide coverage or signal intensities.
        - **Application**: Used to assess the overall relationship between samples, detect outliers, and evaluate batch effects or experimental reproducibility in End-seq or ChIP-seq datasets.

	  *\<p align="center"\> [PCA photo \</p\<img width="709" height="701" alt="PCA" src="https://github.com/user-attachments/assets/22a29886-cbb8-4b0c-abe7-f33909b89350" />
>*

      - **`*.dedup.bam.bai`**

        - **Content**: This is the index file for the BAM file.
        - **Application**: It allows for fast random access to the BAM file, which is essential for visualization software (like IGV) to quickly load and display alignments for a specific genomic region without reading the entire file.

      - **`*.dedup.bam.bai`**

        - **Content**: This is the index file for the BAM file.
        - **Application**: It allows for fast random access to the BAM file, which is essential for visualization software (like IGV) to quickly load and display alignments for a specific genomic region without reading the entire file.

      - **`*.dedup.bam.bai`**

        - **Content**: This is the index file for the BAM file.
        - **Application**: It allows for fast random access to the BAM file, which is essential for visualization software (like IGV) to quickly load and display alignments for a specific genomic region without reading the entire file.

      - **`*.dedup.bam.bai`**

        - **Content**: This is the index file for the BAM file.
        - **Application**: It allows for fast random access to the BAM file, which is essential for visualization software (like IGV) to quickly load and display alignments for a specific genomic region without reading the entire file.

      - **`*.dedup.bam.bai`**

        - **Content**: This is the index file for the BAM file.
        - **Application**: It allows for fast random access to the BAM file, which is essential for visualization software (like IGV) to quickly load and display alignments for a specific genomic region without reading the entire file.

      - **`*.dedup.bam.bai`**

        - **Content**: This is the index file for the BAM file.
        - **Application**: It allows for fast random access to the BAM file, which is essential for visualization software (like IGV) to quickly load and display alignments for a specific genomic region without reading the entire file.

      - **`*.dedup.bam.bai`**

        - **Content**: This is the index file for the BAM file.
        - **Application**: It allows for fast random access to the BAM file, which is essential for visualization software (like IGV) to quickly load and display alignments for a specific genomic region without reading the entire file.

      - **`*.dedup.bam.bai`**

        - **Content**: This is the index file for the BAM file.
        - **Application**: It allows for fast random access to the BAM file, which is essential for visualization software (like IGV) to quickly load and display alignments for a specific genomic region without reading the entire file.

      - **`*.dedup.bam.bai`**

        - **Content**: This is the index file for the BAM file.
        - **Application**: It allows for fast random access to the BAM file, which is essential for visualization software (like IGV) to quickly load and display alignments for a specific genomic region without reading the entire file.

      - **`*.dedup.bam.bai`**

        - **Content**: This is the index file for the BAM file.
        - **Application**: It allows for fast random access to the BAM file, which is essential for visualization software (like IGV) to quickly load and display alignments for a specific genomic region without reading the entire file.



#### Aggregate Result Files

  - **`deg_results.txt`**

      - **Content**: A tab-separated text file containing the results of the differential expression analysis from DESeq2. Each row corresponds to a gene, and columns typically include:
          - `baseMean`: Average normalized count across all samples.
          - `log2FoldChange`: The logarithm (base 2) of the fold change between the 'Treated' and 'Control' conditions. A positive value means the gene is upregulated in the 'Treated' group; a negative value means it is downregulated.
          - `lfcSE`: The standard error of the `log2FoldChange` estimate.
          - `stat`: The Wald statistic.
          - `pvalue`: The raw p-value for the statistical test.
          - `padj`: The p-value adjusted for multiple testing (e.g., using Benjamini-Hochberg correction).
      - **Application**: This is the final result file. You can filter this file based on `padj` (e.g., `padj < 0.05`) and `log2FoldChange` thresholds to obtain a list of statistically significant differentially expressed genes.

  - **`normalized_counts.txt`**

      - **Content**: A tab-separated text file containing a matrix of normalized expression counts. Rows represent genes, and columns represent samples. These counts are adjusted for differences in sequencing depth between libraries, making them comparable across samples.
          - `row names`: Gene name.
          - `Control_Rep1`: The normalized count for Control_Rep1.
          - `Control_Rep2`: The normalized count for Control_Rep2.
          - `Treated_Rep1`: The normalized count for Treated_Rep1.
          - `Treated_Rep2`: The normalized count for Treated_Rep2.
      - **Application**: This matrix is essential for downstream analyses and visualizations beyond simple DEG lists. It can be used as input for generating heatmaps, performing principal component analysis (PCA) to check for sample clustering, or conducting gene set enrichment analysis.

  - **`multiqc_report`** : Open multiqc_report.html in a web browser to explore all sections interactively.

      - **General Statistics**: A combined table summarizing important metrics for each sample:
	  
	  *\<p align="center"\> [qc report photo \</p\>*

      - **FastQC**: Quality-control metrics on raw and trimmed reads, including 'Sequence Counts', 'Sequence Quality Histograms', 'Per Sequence Quality Scores', 'Per Base Sequence Content', 'Per Sequence GC Content', 'Per Base N Content', 'Sequence Length Distribution', 'Sequence Duplication Levels', 'Overrepresented sequences by sample', 'Top overrepresented sequences', 'Adapter Content':
	  
	  Sequence Quality Histograms: The mean quality value across each base position in the read.
	  
	  *\<p align="center"\> [qc report photo \</p\>*
	  
	  Adapter Content: The cumulative percentage count of the proportion of your library which has seen each of the adapter sequences at each position.
	  
	  *\<p align="center"\> [qc report photo \</p\>*
	  
      - **Cutadapt**: Reports the number of reads and bases trimmed for adapters and quality:
	  
	  *\<p align="center"\> [qc report photo \</p\>*
	  
      - **STAR**: Alignment statistics such as total reads, uniquely mapped reads, and multi-mapping rates:
	  
	  *\<p align="center"\> [qc report photo \</p\>*  
	  
      - **featureCounts**: Gene-level quantification results, including total counts and assignment rates:
	  
	  *\<p align="center"\> [qc report photo \</p\>*  
	  
      - **Application**: This is the first file you should check to assess the overall quality of your sequencing data and the alignment process. It helps identify problematic samples (e.g., low alignment rate, high duplication) early on.

After the pipeline completes, the output directory will contain several files and directories. Below is a detailed explanation of what each file is and how it can be used.

-----

### Align Mode Output

```
./project_results/
├── Control_Rep1/
│   ├── Control_Rep1.dedup.bam         # Final processed BAM file
│   ├── Control_Rep1.dedup.bam.bai     # BAM index file
│   └── Control_Rep1.bw                # BigWig signal track
├── Control_Rep2/
│   ├── Control_Rep2.dedup.bam         # Final processed BAM file
│   ├── Control_Rep2.dedup.bam.bai     # BAM index file
│   └── Control_Rep2.bw                # BigWig signal track
├── Treated_Rep1/
│   ├── Treated_Rep1.dedup.bam         # Final processed BAM file
│   ├── Treated_Rep1.dedup.bam.bai     # BAM index file
│   └── Treated_Rep1.bw                # BigWig signal track
├── Treated_Rep2/
│   ├── Treated_Rep2.dedup.bam         # Final processed BAM file
│   ├── Treated_Rep2.dedup.bam.bai     # BAM index file
│   └── Treated_Rep2.bw                # BigWig signal track
├── multiqc_report/
│   └── multiqc_report.html            # Aggregated QC report for all samples
├── deg_results.txt                    # Differential expression gene list from DESeq2
└── normalized_counts.txt              # Normalized counts matrix from DESeq2
```






