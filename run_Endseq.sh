#!/bin/bash


print_help() {
  echo "Usage: $0 --treatmentFq <treatment_fastq> --controlFq <control_fastq> --outputdir <output_directory> --referencedir <reference_index_directory> --adapterFa <illumina_adapter> --sif <singularity_environment> --threads <number_of_threads> --binSize <binSize_for_bw> --g <hs>"
  echo "--treatmentFq        : Path to the treatment fastq"
  echo "--controlFq          : (optinal) Path to the input control fastq"
  echo "--outputdir          : path to the directory where the output will be stored"
  echo "--referencedir       : path to the directory where bowtie reference build with prefix"
  echo "--adapterFa          : path to the adapter fasta"
  echo "--sif                : path to the singularity environment file"
  echo "--threads            : (optinal) Number of threads to use, default 8"
  echo "--binSize            : (optinal) Number of binsize to use, default 10"
  echo "--g                  : Specise, default hs"
}

if [[ "$#" -eq 0 || "$1" == "--help" ]]; then
  print_help
  exit 0
fi

treatmentFq=""
controlFq=""
outputdir=""
referencedir=""
adapterFa=""
sif=""
threads=8
binSize=10
g=hs

while [[ $# -gt 0 ]]; do
  case $1 in
    --treatmentFq) 
      if [[ -z "$2" ||! -f "$2" ]]; then
        echo "Error: fastq file $2 not found!"
        exit 1
      fi
      treatmentFq="$2"; shift 2 ;;
    --controlFq)
      if [[ -z "$2" ]]; then
        controlFq=""
        shift 1
      else
        if [[ ! -f "$2" ]]; then
          echo "Error: fastq file $2 not found!"
          exit 1
        fi
        controlFq="$2"
        shift 2
      fi;;
    --outputdir) 
      if [[ -z "$2" ]]; then
        echo "Error: output directory $2 not found!"
        exit 1
      fi
      outputdir="$2"; shift 2 ;;
    --referencedir) 
      if [[ -z "$2" ]]; then
        echo "Error: reference directory $2 not found!"
        exit 1
      fi
      referencedir="$2"; shift 2 ;;
    --adapterFa) 
      if [[ -z "$2" || ! -f "$2" ]]; then
        echo "Error: adapter fa file $2 not found!"
        exit 1
      fi
      adapterFa="$2"; shift 2 ;;
    --sif) 
      if [[ -z "$2" || ! -f "$2" ]]; then
        echo "Error: singularity sif file $2 not found!"
        exit 1
      fi
      sif="$2"; shift 2 ;;
    --threads) 
      threads="$2"; shift 2 ;;
    --binSize) 
      binSize="$2"; shift 2 ;;
    --g) 
      g="$2"; shift 2 ;;
    *) 
      echo "Unknown option $1"; print_help; exit 1 ;;
  esac
done

if [[ -z "$treatmentFq" || -z "$outputdir" || -z "$referencedir" || -z "$adapterFa" || -z "$sif" ]]; then
  echo "Error: Missing required parameters"
  print_help
  exit 1
fi

########### mkdir
SING_EXEC1="singularity exec --cleanenv $sif"
$SING_EXEC1 mkdir -p ${outputdir}
$SING_EXEC1 mkdir -p ${outputdir}/qc
FastQCdir=${outputdir}/qc
$SING_EXEC1 mkdir -p ${outputdir}/trim
trimedfadir=${outputdir}/trim
$SING_EXEC1 mkdir -p ${outputdir}/bam
bowtieoutdir=${outputdir}/bam
$SING_EXEC1 mkdir -p ${outputdir}/bw
bwoutdir=${outputdir}/bw
$SING_EXEC1 mkdir -p ${outputdir}/figure
FigureDir=${outputdir}/figure
$SING_EXEC1 mkdir -p ${outputdir}/peak
pcoutdir=${outputdir}/peak
$SING_EXEC1 mkdir -p ${outputdir}/multiqc
multiqcdir=${outputdir}/multiqc

########### singularity command

[[ -n "$treatmentFq" ]] && treatmentFq=$(readlink -f "$treatmentFq")
[[ -n "$controlFq" ]] && controlFq=$(readlink -f "$controlFq")
[[ -n "$outputdir" ]] && outputdir=$(readlink -f "$outputdir")
[[ -n "$referencedir" ]] && referencedir=$(readlink -f "$referencedir")
[[ -n "$adapterFa" ]] && adapterFa=$(readlink -f "$adapterFa")
[[ -n "$sif" ]] && sif=$(readlink -f "$sif")
bind_dirs=()
[[ -n "$treatmentFq" ]] && bind_dirs+=("$(dirname "$treatmentFq")")
[[ -n "$controlFq" ]] && bind_dirs+=("$(dirname "$controlFq")")
[[ -n "$outputdir" ]] && bind_dirs+=("$outputdir")
[[ -n "$referencedir" ]] && bind_dirs+=("$(dirname "$referencedir")")
[[ -n "$adapterFa" ]] && bind_dirs+=("$(dirname "$adapterFa")")
bind_dirs_unique=($(printf "%s\n" "${bind_dirs[@]}" | sort -u))
SING_EXEC="singularity exec --cleanenv"
for dir in "${bind_dirs_unique[@]}"; do
    real_dir=$(readlink -f "$dir")
    [[ -n "$real_dir" ]] && SING_EXEC+=" -B $real_dir:$real_dir"
done
SING_EXEC+=" $sif"

########### gzip
if [[ "$treatmentFq" != *.gz ]]; then
    echo "gzip... ${treatmentFq}"
    $SING_EXEC gzip "$treatmentFq"
    treatmentFq=$treatmentFq.gz
fi
if [[ -n "$controlFq" && "$controlFq" != *.gz ]]; then
    echo "gzip... ${controlFq}"
     $SING_EXEC gzip "$controlFq"
    controlFq=$controlFq.gz
fi

########### processing
files=()
[[ -n "$treatmentFq" ]] && files+=("$treatmentFq")
[[ -n "$controlFq" ]] && files+=("$controlFq")

for fqpath in "${files[@]}";
do
	fileName=$(basename "$fqpath" | sed -E 's/\.f(ast)?q\.gz$//')
	echo "run fastqc"
	$SING_EXEC fastqc -o ${FastQCdir} -t ${threads} -q $fqpath
	echo "cut adapters and sequence quality filters"
	$SING_EXEC cutadapt -q 25 -m 36 -e 0.1 -j ${threads} -a G{8} -a A{8} -g G{8} -g A{8} -a file:${adapterFa} --times 13 -o ${trimedfadir}/${fileName}.trimmed.fastq.gz $fqpath
	echo "the second quality control of these files"
	$SING_EXEC fastqc -o ${trimedfadir} --noextract -f fastq ${trimedfadir}/${fileName}.trimmed.fastq.gz
	echo "mapped the data to the genome"
	$SING_EXEC bowtie -S -p ${threads} -x ${referencedir} ${trimedfadir}/${fileName}.trimmed.fastq.gz --best --strata --all -m 1 -n 3 -l 50  2> ${bowtieoutdir}/${fileName}.bowtie.stats  | samtools view -q 255 -bS - | samtools sort -n -@ 20 - -T ${fileName} -O BAM -o ${bowtieoutdir}/${fileName}.sort.bam
	echo "rm duplicate reads"
	$SING_EXEC samtools fixmate -r -m ${bowtieoutdir}/${fileName}.sort.bam -| samtools sort -@ ${threads} - | samtools markdup -r -s - ${bowtieoutdir}/${fileName}.DeDup.bam 2> ${bowtieoutdir}/${fileName}.markdup.log
	$SING_EXEC samtools flagstat ${bowtieoutdir}/${fileName}.DeDup.bam > ${bowtieoutdir}/${fileName}.flagstat.txt
	$SING_EXEC samtools index ${bowtieoutdir}/${fileName}.DeDup.bam
	rm ${bowtieoutdir}/${fileName}.sort.bam
	echo "obtain BW files"
	$SING_EXEC samtools index ${bowtieoutdir}/${fileName}.DeDup.bam
	$SING_EXEC bamCoverage -b ${bowtieoutdir}/${fileName}.DeDup.bam -o ${bwoutdir}/${fileName}.DeDup.bw  --ignoreForNormalization chrM -p ${threads} --binSize ${binSize} --normalizeUsing CPM
done

########### correlation
echo "Check the correlation"
bwfiles=("${bwoutdir}"/*.bw)
labels=($(basename -s .bw "${bwfiles[@]}"))
$SING_EXEC multiBigwigSummary bins -b "${bwfiles[@]}" --labels "${labels[@]}" -out ${FigureDir}/BW_compare.npz -p ${threads}
$SING_EXEC plotPCA -in ${FigureDir}/BW_compare.npz -o ${FigureDir}/BW_compare.pdf
$SING_EXEC plotCorrelation -in ${FigureDir}/BW_compare.npz --corMethod pearson --skipZeros -o ${FigureDir}/BW_compare_cor.pdf --whatToPlot heatmap --colorMap RdYlBu --plotNumbers

echo "check the End-seq quality"
bamfiles=("${bowtieoutdir}"/*.DeDup.bam)
$SING_EXEC plotFingerprint -b "${bamfiles[@]}" --labels "${labels[@]}" --skipZeros --plotFile ${FigureDir}/fingerprints.pdf

echo "peak calling"
namep=$(basename "$treatmentFq" | sed -E 's/\.f(ast)?q\.gz$//')
if [[ -z "$controlFq" ]]; then
	$SING_EXEC macs3 callpeak -t ${bowtieoutdir}/${namep}.DeDup.bam -f BAM -g ${g} --outdir ${pcoutdir} -n ${namep} --broad --nomodel --nolambda --llocal 100000 --keep-dup all -q 0.01 2> >(tee -a ${pcoutdir}/${namep}.macs3.stats >&2) 1> ${pcoutdir}/${namep}.macs3.stdout
else
	namec=$(basename "$controlFq" | sed -E 's/\.f(ast)?q\.gz$//')
	$SING_EXEC macs3 callpeak -t ${bowtieoutdir}/${namep}.DeDup.bam -c ${bowtieoutdir}/${namec}.DeDup.bam -f BAM -g ${g} --outdir ${pcoutdir} -n ${namep} --broad --nomodel --nolambda --keep-dup all -q 0.01 2> >(tee -a ${pcoutdir}/${namep}.macs3.stats >&2) 1> ${pcoutdir}/${namep}.macs3.stdout
fi

echo "peak enrichment"
bwfiles=(${bwoutdir}/*.bw)
$SING_EXEC computeMatrix scale-regions -p ${threads} -S ${bwfiles[@]} \
                              -R ${pcoutdir}/${namep}_peaks.broadPeak -a 5000 -b 5000 --regionBodyLength 1000 \
                              --missingDataAsZero -o ${FigureDir}/${namep}.peak.gz
$SING_EXEC plotHeatmap -m ${FigureDir}/${namep}.peak.gz -out ${FigureDir}/${namep}.peak.pdf --colorMap viridis --missingDataColor white --heatmapHeight 12 --heatmapWidth 4 --startLabel "Start" --endLabel "End"

########### multiqc
$SING_EXEC multiqc ${FastQCdir}/* ${trimedfadir}/* ${bowtieoutdir}/*.bowtie.stats ${bowtieoutdir}/*flagstat.txt -o ${multiqcdir} --force

########### remove
rm -r ${FastQCdir}
rm -r ${trimedfadir}
rm ${FigureDir}/${namep}.peak.gz
rm ${pcoutdir}/${namep}.macs3.stdout

echo "Finished."

