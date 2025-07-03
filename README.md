# End-seq-Processing-Pipeline
This is the pipeline of downstream analysis for End-seq data, including quality contrl (QC) and peak calling.
# Part I Introduction
## i. Workflow
Here stands an throughout workflow of End-seq data analysis.
![图片](https://github.com/user-attachments/assets/02d8de2f-e907-4273-8fe5-68f4d5783ab5)

As illustrated in the figure,
(i) Solid-line boxes represent the files we input/output;
(ii) Dashed-line boxes represent the QC results after processing at each step.

We will proceed by structuring our workflow according to (ii).

## ii. File Structure
```
End-seq
├─ 1.rawdata
│    ├─ wm_13AID_AUX_shCTRL_S29_L001_R1_001.fastq.gz
│    ├─ wm_13AID_AUX_shCTRL_S29_L001_R2_001.fastq.gz
│    ├─ wm_13AID_AUX_shRAD21_S30_L001_R1_001.fastq.gz
│    └─ wm_13AID_AUX_shRAD21_S30_L001_R2_001.fastq.gz
├─ 2.trim
│    ├─ wm_13AID_AUX_shCTRL_S29_L001_R1_trimmed.fastq.gz
│    └─ wm_13AID_AUX_shRAD21_S30_L001_R1_trimmed.fastq.gz
├─ 3.bowtie 
│    ├─ wm_13AID_AUX_shCTRL_S29_L001.bowtie.stats
│    ├─ wm_13AID_AUX_shCTRL_S29_L001.sort.bam
│    ├─ wm_13AID_AUX_shCTRL_S29_L001.sort.bam.bai
│    ├─ wm_13AID_AUX_shCTRL_S29_L001_marked_dup_metrics.txt
│    ├─ wm_13AID_AUX_shCTRL_S29_L001_removeDup.bam
│    ├─ wm_13AID_AUX_shCTRL_S29_L001_removeDup.bam.bai
│    ├─ wm_13AID_AUX_shRAD21_S30_L001.bowtie.stats
│    ├─ wm_13AID_AUX_shRAD21_S30_L001.sort.bam
│    ├─ wm_13AID_AUX_shRAD21_S30_L001.sort.bam.bai
│    ├─ wm_13AID_AUX_shRAD21_S30_L001_marked_dup_metrics.txt
│    ├─ wm_13AID_AUX_shRAD21_S30_L001_removeDup.bam
│    └─ wm_13AID_AUX_shRAD21_S30_L001_removeDup.bam.bai
├─ 4.bw
│    ├─ wm_13AID_AUX_shCTRL_S29_L001_removeDup_rmChrM.bw
│    └─ wm_13AID_AUX_shRAD21_S30_L001_removeDup_rmChrM.bw
└─ 5.macs3
```
## iii. Software Required
```
cutadapt
fastqc
bowtie1
picard
deeptools
macs3
```
# Part II Codes
## Preprocessing
In this section, you will convert the raw FASTQ files into bam files, which can be used for subsequent analysis.
```
#The location of rawdata
datadir=/1.rawdata
#the QC of raw data
FastQCdir=/1.rawdata
#Trimmed data and QC
trimedfadir=/2.trim
#Mapped data to the human genome
bowtieoutdir=/3.bowtie
hg38ref=XXX
#obtain bigwig files
bwoutdir=/4.bw
cd ${datadir}
fastqc -o ${FastQCdir} -t 8 -q *.fastq.gz
##
for fileName in `ls -1 *.fastq.gz |sed 's/_R1/\t/g'|sed 's/_R2/\t/g'|cut -f 1|sort|uniq`;
do
	echo "cut adapters and sequence quality filters"
	cutadapt -q 25 -m 36 -e 0.1 -j 30 -a G{8} -a A{8} -g G{8} -g A{8} -a file:illumina_adapter.fa --times 13 -o ${trimedfadir}/${fileName}_R1_trimmed.fastq.gz ${datadir}/${fileName}_R1_001.fastq.gz
	echo "the second quality control of these files"
	fastqc -o ${trimedfadir} --noextract -f fastq ${trimedfadir}/${fileName}_R1_trimmed.fastq.gz
	echo "mapped the data to the genome"
	bowtie -S -p 20 -x ${hg38ref} ${trimedfadir}/${fileName}_R1_trimmed.fastq.gz --best --strata --all -m 1 -n 3 -l 50  2> ${bowtieoutdir}/${fileName}.bowtie.stats  | samtools view -q 255 -bS - | samtools sort - -o ${bowtieoutdir}/${fileName}.sort.bam
	samtools index ${bowtieoutdir}/${fileName}.sort.bam
	echo "rm duplicate reads"
	java -jar /mnt/share/software/picard/picard.jar  MarkDuplicates I=${bowtieoutdir}/${fileName}.sort.bam O=${bowtieoutdir}/${fileName}_removeDup.bam M=${bowtieoutdir}/${fileName}_marked_dup_metrics.txt REMOVE_DUPLICATES=true
        echo "obtain BW files"
	samtools index ${bowtieoutdir}/${fileName}_removeDup.bam
        bamCoverage -b ${bowtieoutdir}/${fileName}_removeDup.bam -o ${bowtieoutdir}/${fileName}_removeDup_rmChrM.bw  --ignoreForNormalization chrM -p 10 --binSize 1 --ignoreDuplicates --normalizeUsing CPM
done
```
## Check the End-seq quality
```
bwoutdir=/4.bw
FigureDir=/4.bw
echo "Check the correlation among replicates in HCT116 cells"
multiBigwigSummary bins -b ${bwoutdir}/wm_13AID_AUX_shCTRL_S29_L001_removeDup_rmChrM.bw ${bwoutdir}/wm_13AID_AUX_shRAD21_S30_L001_removeDup_rmChrM.bw --labels AUD_shCTRL AUD_shRAD21 -out ${FigureDir}/BW_compare_HCT116.npz -p 30
plotPCA -in ${FigureDir}/BW_compare_HCT116.npz -o ${FigureDir}/BW_compare_HCT116_PCA.pdf
plotCorrelation  -in ${FigureDir}/BW_compare_HCT116.npz --corMethod pearson --skipZeros -o ${FigureDir}/BW_compare_HCT116_cor.pdf --whatToPlot heatmap --colorMap RdYlBu --plotNumbers

echo "check the End-seq quality"
bowtieoutdir=/3.bowtie
plotFingerprint -b ${bowtieoutdir}/wm_13AID_AUX_shRAD21_S30_L001_removeDup.bam ${bowtieoutdir}/wm_13AID_AUX_shCTRL_S29_L001_removeDup.bam --labels AUD_shCTRL AUD_shRAD21 --skipZeros --plotFile ${FigureDir}/fingerprints.png

```
## Peak calling with QC
The plotFingerprint curve showed a sharp increase toward the upper right corner, indicating strong signal enrichment in a small subset of genomic bins, which is characteristic of narrow peaks. 
END-seq data typically displays sharp, localized enrichment patterns consistent with narrow peaks; thus, it is standard practice to use MACS3 in narrow peak mode, and plotFingerprint can help confirm signal enrichment but is not required to decide peak type.
```
echo "peak calling"
bwoutdir=/4.bw
bowtieoutdir=/3.bowtie
pcoutdir=/5.macs3
for namep in 13AID_AUX_shCTRL_S29 13AID_AUX_shRAD21_S30;do
  macs3 callpeak -t ${bowtieoutdir}/wm_${namep}_L001_removeDup.bam -f BAM -g hs --outdir ${pcoutdir} -n ${namep} --nomodel --nolambda --llocal 100000 --keep-dup all -q 0.01 2> >(tee -a ${pcoutdir}/${namep}.macs1.stats >&2) 1> ${pcoutdir}/${namep}.macs1.stdout
  ## peak_enrichment.pdf
  computeMatrix scale-regions -p 20 -S ${bwoutdir}/wm_${namep}_L001_removeDup_rmChrM.bw \
                              -R ${pcoutdir}/${namep}_peaks.narrowPeak -a 5000 -b 5000 --regionBodyLength 0 \
                              --missingDataAsZero -o ${pcoutdir}/${namep}.peak1.5kb.gz
  plotHeatmap -m ${pcoutdir}/${namep}.peak1.5kb.gz -out ${pcoutdir}/${namep}.peak1.5kb.pdf --colorMap viridis --missingDataColor white --heatmapHeight 12 --heatmapWidth 4
  ## CTCF_enrichment.pdf
  CTCFfile=XXX ### path to CTCF bw
  computeMatrix scale-regions -p 20 -S ${bwoutdir}/${CTCFfile} \
                              -R ${pcoutdir}/${namep}_peaks.narrowPeak -a 5000 -b 5000 --regionBodyLength 0 \
                              --missingDataAsZero -o ${pcoutdir}/${namep}.CTCF.5kb.gz
  plotHeatmap -m ${pcoutdir}/${namep}.CTCF.5kb.gz -out ${pcoutdir}/${namep}.CTCF.5kb.pdf --colorMap viridis --missingDataColor white --heatmapHeight 12 --heatmapWidth 4
  ## Top2_enrichment.pdf
  Top2file=XXX ### path to Top2 bw
  computeMatrix scale-regions -p 20 -S ${bwoutdir}/${Top2file} \
                              -R ${pcoutdir}/${namep}_peaks.narrowPeak -a 5000 -b 5000 --regionBodyLength 0 \
                              --missingDataAsZero -o ${pcoutdir}/${namep}.Top2.5kb.gz
  plotHeatmap -m ${pcoutdir}/${namep}.Top2.5kb.gz -out ${pcoutdir}/${namep}.Top2.5kb.pdf --colorMap viridis --missingDataColor white --heatmapHeight 12 --heatmapWidth 4
  ## TSS_enrichment.pdf
  TSSfile=XXX ### path to Nacent RNA-seq bw
  computeMatrix scale-regions -p 20 -S ${bwoutdir}/${TSSfile} \
                              -R ${pcoutdir}/${namep}_peaks.narrowPeak -a 5000 -b 5000 --regionBodyLength 0 \
                              --missingDataAsZero -o ${pcoutdir}/${namep}.TSS.5kb.gz
  plotHeatmap -m ${pcoutdir}/${namep}.TSS.5kb.gz -out ${pcoutdir}/${namep}.TSS.5kb.pdf --colorMap viridis --missingDataColor white --heatmapHeight 12 --heatmapWidth 4
  ## H3K27ac_enrichment.pdf
  H3K27acfile=XXX ### path to Nacent RNA-seq bw
  computeMatrix scale-regions -p 20 -S ${bwoutdir}/${H3K27acfile} \
                              -R ${pcoutdir}/${namep}_peaks.narrowPeak -a 5000 -b 5000 --regionBodyLength 0 \
                              --missingDataAsZero -o ${pcoutdir}/${namep}.H3K27ac.5kb.gz
  plotHeatmap -m ${pcoutdir}/${namep}.H3K27ac.5kb.gz -out ${pcoutdir}/${namep}.H3K27ac.5kb.pdf --colorMap viridis --missingDataColor white --heatmapHeight 12 --heatmapWidth 4
done
```
















