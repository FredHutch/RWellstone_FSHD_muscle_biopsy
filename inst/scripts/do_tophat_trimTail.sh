#!/bin/bash
#./do_tophat.sh
#SBATCH -n6 -t 1-0 -p campus --mail-type=ALL --mail-user=cwon2@fhcrc.org -A tapscott_s

sampleDir=$1
sampleName=$2
pkgDir=$3
tophatOutDir=$pkgDir/tophat
trimmedDir=$pkgDir/trimmed
fastqcDir=$pkgDir/fastqc
filteredDir=$pkgDir/filtered

echo $sampleName

bowtieVersion=/home/solexa/apps/bowtie/bowtie2-2.2.6
tophatVersion=/home/solexa/apps/tophat/tophat-2.1.0.Linux_x86_64
gtfFile=/shared/biodata/ngs/Reference/iGenomes/Homo_sapiens/UCSC/hg38/Annotation/Archives/archive-2015-08-14-08-18-15/Genes/genes.gtf
genomeBuild=/shared/biodata/ngs/Reference/iGenomes/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome
maxIntronLength=500000
innerDist=80
segmentLength=25
libraryType=fr-unstranded

#umask 0002
export PATH=/usr/kerberos/bin:/usr/local/bin:/bin:/usr/bin:/usr/X11R6/bin:/opt/moab/bin:$bowtieVersion:$tophatVersion:/home/solexa/apps/samtools/samtools-0.1.19:/home/solexa/apps/FastQC:/home/solexa/apps/trim_galore:/app/cutadapt/1.1/bin

#remove /app/bin from the PATH, as it will default to running bowtie 0.12.8
export PATH=${PATH/\/app\/bin:}

#define some directories
mkdir -p $filteredDir
filterPerSample=$filteredDir/$sampleName
mkdir $filterPerSample
mkdir -p $fastqcDir
mkdir -p $trimmedDir
mkdir -p $tophatOutDir
bamDir=$tophatOutDir/bam
mkdir -p $bamDir
tmpDir=$tophatOutDir/$sampleName
mkdir -p $tmpDir

# (1) by pass filter
#cd $sampleDir
#for i in *fastq.gz
#    do
#        i2=${i//.gz/}
#        zgrep -A 3 '^@.*[^:]*:N:[^:]*:' $i | zgrep -v '^\-\-$' > $filterPerSample/#$i2
#done
#gzip -r $filterPerSample

# (2) zcat all the fastq.gz 
cd $filterPerSample
fqFile=$trimmedDir/$sampleName.filterd.fastq
zcat $(ls *.fastq.gz) > $fqFile

# (3) trim adapter and tail (TruSeq3-SE.fa is saved at $trimmedDir)
cd $trimmedDir
trimFq=$trimmedDir/$sampleName.trim.fastq
java -jar /home/solexa/apps/trimmomatic/Trimmomatic-0.32/trimmomatic-0.32.jar SE -threads 6 $fqFile $trimFq ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 TRAILING:3 MINLEN:36 
rm $fqFile

# (5) fastqc again
fastqc -t 6 $trimFq -o $fastqcDir --casava
touch $fastqcDir/$sampleName\_fastqcDone.txt

# (6) tophat alignment
tophat --mate-inner-dist $innerDist --num-threads 4 -G $gtfFile --library-type $libraryType -I $maxIntronLength --segment-length $segmentLength --no-coverage-search -o $tmpDir $genomeBuild $trimFq

mv $tmpDir/accepted_hits.bam $bamDir/$sampleName.bam

# sort bam
cd $bamDir
samtools sort -@ 4 $sampleName.bam $sampleName.bam.sorted
mv $sampleName.bam.sorted.bam $sampleName.bam
samtools index $sampleName.bam
touch $sampleName.tophatDone.txt
exit 0


