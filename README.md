# EEB5300
Files for collaborative project on northern sandlance population genetics

# Data processing for low coverage whole genome sequencing 

Scripts for going from raw .fastq files to processed and quality-checked .bam files for downstream analysis


## Before you start

#### Create a project directory

As a first step, you should create a project directory (referred to as `BASEDIR` in certain scripts), with the following subdirectories:

  * `adapter_clipped` adapter clipped fastq files
  
  * `bam` bam files
  
  * `nohups` nohup log files
  
  * `sample_lists` sample tables and sample lists
   
  * `markdowns` markdown files tracking the workflow; it is recommended that you first create Rmd files using Rstudio, and knit these into GitHub markdown format

The following subdirectories are optional, depending on your workflow:

  * `angsd` angsd output files
  
  * `fastq` raw fastq files (if you don't already have it backed up on the backup folder)
  
  * `species_stats` fastq species detector output
  
  * `fastqc` FastQC output
  
  * `qual_filtered` quality filtered and/or polyG trimming
    
  * `scripts` scripts specific to your projects (e.g. merging certain bam files)

It is recommended that you make this project directory a separate GitHub repository. This way, you can track all the changes in your project directory. To make GitHub ignore certain directories (e.g. `adapter_clipped` and `bam`, since these will contain many large files), create a `.gitignore` file in your project directory, in which you specify the paths that should be ignored by git. 

#### Prepare sample lists

A sample list is a list of the prefixes of raw fastq files. These should be unique for each fastq file, and the rest of the names have to be the same for all raw fastq files. No header should be included in this list. 

You can have more than one sample list for each project; by doing so you will be able to run some of the scripts in parallel. When you have different suffixes for your raw fastq files (e.g when sequences come from different sequencing platforms), you will have to create multiple sample lists. 

Here is an example of a sample list: 

#### Prepare a sample table

A sample table is a **tab deliminated** table that includes relevant information for all fastq files. It should include the following six columns, strictly in this order:

  * `prefix` prefix of raw fastq files; it should be the union of all your sample lists
  
  * `lane_number` lane number; each sequencing lane should be assigned a different number
  
  * `seq_id` sequence ID，this can be the same thing as sample ID or lane ID and it does not matter except for when different libraries were prepared out of the same sample and were run in the same lane. In this case, seq_id should be used to distinguish these.
  
  * `sample_id` sample ID
  
  * `population` population name
  
  * `data_type` data type; there can only be two possible entries: `pe` (for paired-end data) or `se` (for single end data)

It is important to make sure that the combination of lane_number, seq_id, and sample_id has to be unique for each fastq file. 

It is recommended that you have one single sample table per project. 

An example of a sample table: 

An example of how a sample table can be created from a list of fastq file names and the original sample infomation sheet: 

#### Install programs used

If you are not working on the Therkildsen server, you might need to intall the following programs to your machine. The program paths will then need to be inputted to each script (see script headers).
 
 * `Trimmomatic` http://www.usadellab.org/cms/?page=trimmomatic
 * `fastp` https://github.com/OpenGene/fastp 
 * `bowtie2` http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
 * `Picard` https://broadinstitute.github.io/picard/
 * `bamUtil` https://github.com/statgen/bamUtil
 * `GenomeAnalysisTK-3.7` https://software.broadinstitute.org/gatk/documentation/version-history.php?id=8692&page=3
 * `BBMap` https://github.com/BioInfoTools/BBMap

## Demultiplex

If your fastq file has not been demultiplex and if the barcodes are part of the fastq headers, use the following line to demultiplex. Replace the items in quotes with appropriate paths and names. An example of this workflow can be found [here]. 

Each line in the barcode list should be in the format of `i7sequence+i5sequence`, and an example of such a barcode list can be found [here].

See [demuxbyname.sh](https://github.com/BioInfoTools/BBMap/blob/master/sh/demuxbyname.sh) for details. 

``` bash
nohup bash /programs/bbmap-38.45/demuxbyname.sh \
in="Path to gzipped fastq files" \
in2="Path to gzipped fastq files if your have pair-end reads" \
out="Suffix of output names; this should start with a percentage sign (%)" \
out2="Suffix of output names if your have pair-end reads; this should start with a percentage sign (%)" \
outu=unknown_barcode_1.fastq.gz \
outu2=unknown_barcode_2.fastq.gz \
prefixmode=f \
names="Path to a list of barcode sequences" \
> "Path to the nohup file" &
```

After this is run, make sure to rename your demulitplexed fastq files so that the file names do not contain a plus sign (`+`), which can be a problem in downstream analyses. 

## Get started


1. Clip adapters using [adapter_clipping.sh]

2. Quality filtering using [quality_filtering.sh]

3. Build bowtie reference index using [build_bowtie_ref_index.sh]. This only needs to be done once with the same reference genome.

4. Map to reference, sort, and quality filter using [low_coverage_mapping.sh]

5. Merge samples that were sequenced multiple times on different lanes (or within the same lane, with different barcodes). You should write your own script to do this

6. Deduplicate (all samples) and clip overlapping read pairs (pair-end only) 


## Optional steps

1. Check contamination or species composition in fastq files using [run_species_detector.sh]

2. Count fastq files using [count_fastq.sh]

3. Count bam files before merging using [count_bam_unmerged.sh]

4. Count bam files after merging using [count_bam_merged.sh]

5. Count per position depth using [count_depth_per_position_per_sample.sh] and summarize these using [summarize_depth_per_position.R]


## Scripts for data processing - One lane given per "job" 

Demultiplex shell

```{bash}
nohup bash /programs/bbmap-37.41/demuxbyname.sh \
in="Path to gzipped fastq files" \
in2="Path to gzipped fastq files if your have pair-end reads" \
out="Suffix of output names; this should start with a percentage sign (%)" \
out2="Suffix of output names if your have pair-end reads; this should start with a percentage sign (%)" \
outu=unknown_barcode_1.fastq.gz \
outu2=unknown_barcode_2.fastq.gz \
prefixmode=f \
names="Path to a list of barcode sequences" \
> "Path to the nohup file" &

```

Demultiplex job
```{bash}
#!/bin/bash
#SBATCH --job-name=demultiplexL5
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=5G
#SBATCH --mail-user=lucas.jones@uconn.edu
#SBATCH -o demultiplex_%j.out
#SBATCH -e demultiplex_%j.err


bash /isg/shared/apps/bbmap/37.41/demuxbyname.sh \
in=/labs/Baumann/BaumannLab/SLGenetics/128.120.88.251/X202SC20041425-Z01-F001/raw_data/Pool2/Pool2_CKDL200153435-1B_HC25JCCX2_L5_1.fq.gz \
in2=/labs/Baumann/BaumannLab/SLGenetics/128.120.88.251/X202SC20041425-Z01-F001/raw_data/Pool2/Pool2_CKDL200153435-1B_HC25JCCX2_L5_2.fq.gz \
out=%_lane_5_R1.fastq.gz \
out2=%_lane_5_R2.fastq.gz \
outu=unknown_barcode_1.fastq.gz \
outu2=unknown_barcode_2.fastq.gz \
prefixmode=f \
names=/labs/Baumann/BaumannLab/SLGenetics/128.120.88.251/X202SC20041425-Z01-F001/raw_data/Pool2/newpool2.txt
```

Adapter Clipping shell

```#!/bin/bash

## This script is used to clip adapters. It can process both paired end and single end data. 
SAMPLELIST=$1 
SAMPLETABLE=$2 
RAWFASTQDIR=$3 
BASEDIR=$4 
RAWFASTQSUFFIX1=$5 
RAWFASTQSUFFIX2=$6 
ADAPTERS=$7 


## Loop over each sample
for SAMPLEFILE in `cat $SAMPLELIST`; do
	
	
	SAMPLE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 4`
	SEQ_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 3`
	LANE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 2`
	SAMPLE_SEQ_ID=$SAMPLE_ID'_'$SEQ_ID'_'$LANE_ID  # When a sample has been sequenced in multiple lanes, we need to be able to identify the files from each run uniquely
	
	
	DATATYPE=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 6`
	
	
	RAWFASTQ_ID=$RAWFASTQDIR$SAMPLEFILE
	SAMPLEADAPT=$BASEDIR'adapter_clipped/'$SAMPLE_SEQ_ID
	
	
	if [ $DATATYPE = pe ]; then
		java -jar /isg/shared/apps/Trimmomatic/0.39/trimmomatic-0.39.jar PE -threads 18 -phred33 $RAWFASTQ_ID$RAWFASTQSUFFIX1 $RAWFASTQ_ID$RAWFASTQSUFFIX2 $SAMPLEADAPT'_adapter_clipped_f_paired.fastq.gz' $SAMPLEADAPT'_adapter_clipped_f_unpaired.fastq.gz' $SAMPLEADAPT'_adapter_clipped_r_paired.fastq.gz' $SAMPLEADAPT'_adapter_clipped_r_unpaired.fastq.gz' 'ILLUMINACLIP:'$ADAPTERS':2:30:10:1:true'
	
	elif [ $DATATYPE = se ]; then
		java -jar /isg/shared/apps/Trimmomatic/0.39/trimmomatic-0.39.jar SE -threads 18 -phred33 $RAWFASTQ_ID$RAWFASTQSUFFIX1 $SAMPLEADAPT'_adapter_clipped_se.fastq.gz' 'ILLUMINACLIP:'$ADAPTERS':2:30:10'
	fi
	
done
```

Adapter clipping job

```#!/bin/bash
#SBATCH --job-name=fulladapter_clipping.sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 14
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=10G
#SBATCH --mail-user=lucas.jones@uconn.edu
#SBATCH -o fulladapter_clipping_%j.out
#SBATCH -e fulladapter_clipping_%j.err

echo `hostname`

bash /labs/Baumann/BaumannLab/SLGenetics/Scripts/adapter_clipping.sh /labs/Baumann/BaumannLab/SLGenetics/Misc/AdapterClippingLane3List.tsv /labs/Baumann/BaumannLab/SLGenetics/Misc/AdapterClippingLane3Table.tsv /labs/Baumann/BaumannLab/SLGenetics/128.120.88.251/X202SC20041425-Z01-F001/raw_data/Pool1/Demultiplex2/ /labs/Baumann/BaumannLab/SLGenetics/adapter_clipped/ _R1.fastq.gz _R2.fastq.gz /labs/Baumann/BaumannLab/SLGenetics/Misc/NexteraPE_NT.fa 


bash /labs/Baumann/BaumannLab/SLGenetics/Scripts/adapter_clipping.sh /labs/Baumann/BaumannLab/SLGenetics/Misc/AdapterClippingLane4List.tsv /labs/Baumann/BaumannLab/SLGenetics/Misc/AdapterClippingLane4Table.tsv /labs/Baumann/BaumannLab/SLGenetics/128.120.88.251/X202SC20041425-Z01-F001/raw_data/Pool1/Demultiplex2/ /labs/Baumann/BaumannLab/SLGenetics/adapter_clipped/ _R1.fastq.gz _R2.fastq.gz /labs/Baumann/BaumannLab/SLGenetics/Misc/NexteraPE_NT.fa 

bash /labs/Baumann/BaumannLab/SLGenetics/Scripts/adapter_clipping.sh /labs/Baumann/BaumannLab/SLGenetics/Misc/AdapterClippingLane5List.tsv /labs/Baumann/BaumannLab/SLGenetics/Misc/AdapterClippingLane5Table.tsv /labs/Baumann/BaumannLab/SLGenetics/128.120.88.251/X202SC20041425-Z01-F001/raw_data/Pool2/Demultiplexpool2/ /labs/Baumann/BaumannLab/SLGenetics/adapter_clipped/ _R1.fastq.gz _R2.fastq.gz /labs/Baumann/BaumannLab/SLGenetics/Misc/NexteraPE_NT.fa 
```

Edit Masurca Configuration file 
```{bash}
module load R/3.6.3
R
install.packages(tidyverse)
library(tidyverse)
sr_config <-"

DATA

PE= pe 470 52.5  /labs/Baumann/BaumannLab/SLGenetics/128.120.88.251/X202SC20041425-Z01-F001/raw_data/deepseq/128.120.88.251/X202SC20050719-Z01-F001/raw_data//s608_CSFP200004436-1a_H5TLWDSXY_L1_1.fq.gz  /labs/Baumann/BaumannLab/SLGenetics/128.120.88.251/X202SC20041425-Z01-F001/raw_data/deepseq/128.120.88.251/X202SC20050719-Z01-F001/raw_data/s608/s608_CSFP200004436-1a_H5TLWDSXY_L1_2.fq.gz

END

PARAMETERS

EXTEND_JUMP_READS=0
GRAPH_KMER_SIZE = auto
USE_LINKING_MATES = 1
USE_GRID=0
GRID_QUEUE=all.q
GRID_BATCH_SIZE=300000000 
LHE_COVERAGE=30
LIMIT_JUMP_COVERAGE = 300
CA_PARAMETERS =  cgwErrorRate=0.1
KMER_COUNT_THRESHOLD = 1
CLOSE_GAPS=1
NUM_THREADS = 32
JF_SIZE = 16000000000
SOAP_ASSEMBLY=0

END"
write_lines(sr_config, "/labs/Baumann/BaumannLab/SLGenetics/genome/sr_config.txt")

```

Transcript assembly 
```{bash}
module load MaSuRCA

/isg/shared/apps/MaSuRCA/3.3.3/bin/masurca \
/labs/Baumann/BaumannLab/SLGenetics/genome/sr_config.txt
```

Build Bowtie reference 

```{bash}
#!/bin/bash
module load samtools/1.9
## This script is used to build the bow tie reference index. 
# Run this only when working with a new reference that has not been formatted for bowtie2

REFERENCE=$1 # path to reference fasta file and file name,
REFNAME=$2 # reference name to add to output files,

REFBASENAME="${REFERENCE%.*}"

## First create .bai and .dict files if they haven't been created
if [ ! -f $REFERENCE'.fai' ] ; then
	samtools faidx $REFERENCE
fi

if [ ! -f $REFBASENAME'.dict' ] ; then
	java -jar /isg/shared/apps/picard/2.2.1/picard-2.2.1.jar CreateSequenceDictionary R=$REFERENCE O=$REFBASENAME'.dict'
fi

## Build the reference index
bowtie2-build $REFERENCE $REFBASENAME
```
Mapping Shell

```#!/bin/bash
module load bowtie2/2.3.5.1
module load samtools/1.9
## This script is used to quality filter and trim poly g tails. It can process both paired end and single end data.
SAMPLELIST=$1 # Path to a list of prefixes of the raw fastq files. It should be a subset of the the 1st column of the sample table. 
SAMPLETABLE=$2 # Path to a sample table where the 1st column is the prefix of the raw fastq files. The 4th column is the sample ID, the 2nd column is the lane number, and the 3rd column is sequence ID. The combination of these three columns have to be unique. The 6th column should be data type, which is either pe or se. 
FASTQDIR=$3 # Path to the directory where fastq file are stored. 
BASEDIR=$4 # Path to the base directory where adapter clipped fastq file are stored in a subdirectory titled "adapter_clipped" and into which output files will be written to separate subdirectories.
FASTQSUFFIX1=$5 # Suffix to fastq files. Use forward reads with paired-end data. 
FASTQSUFFIX2=$6 # Suffix to fastq files. Use reverse reads with paired-end data. 
MAPPINGPRESET=$7 # The pre-set option to use for mapping in bowtie2 (very-sensitive for end-to-end (global) mapping [typically used when we have a full genome reference], very-sensitive-local for partial read mapping that allows soft-clipping [typically used when mapping genomic reads to a transcriptome]
REFERENCE=$8 # Path to reference fasta file and file name,
REFNAME=$9 # Reference name to add to output files, 

## Loop over each sample
for SAMPLEFILE in `cat $SAMPLELIST`; do

	## Extract relevant values from a table of sample, sequencing, and lane ID (here in columns 4, 3, 2, respectively) for each sequenced library
	SAMPLE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 4`
	SEQ_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 3`
	LANE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 2`
	SAMPLE_SEQ_ID=$SAMPLE_ID'_'$SEQ_ID'_'$LANE_ID

	## Extract data type from the sample table
	DATATYPE=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 6`

	## The input and output path and file prefix
	SAMPLETOMAP=$FASTQDIR$SAMPLE_SEQ_ID
	SAMPLEBAM=$BASEDIR'bam/'$SAMPLE_SEQ_ID

	## Define platform unit (PU), which is the lane number
	PU=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 2`

	## Define reference base name
	REFBASENAME="${REFERENCE%.*}"

	## Map reads to the reference
	# Map the paired-end reads
	if [ $DATATYPE = pe ]; then
	# We ignore the reads that get orphaned during adapter clipping because that is typically a very small proportion of reads. If a large proportion of reads get orphaned (loose their mate so they become single-end), these can be mapped in a separate step and the resulting bam files merged with the paired-end mapped reads.
		bowtie2 -q --phred33 --$MAPPINGPRESET -p 16 -I 0 -X 1500 --fr --rg-id $SAMPLE_SEQ_ID --rg SM:$SAMPLE_ID --rg LB:$SAMPLE_ID --rg PU:$PU --rg PL:ILLUMINA -x $REFBASENAME -1 $SAMPLETOMAP$FASTQSUFFIX1 -2 $SAMPLETOMAP$FASTQSUFFIX2 -S $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'.sam'

	# Map the single-end reads
	elif [ $DATATYPE = se ]; then
		bowtie2 -q --phred33 --$MAPPINGPRESET -p 16 --rg-id $SAMPLE_SEQ_ID --rg SM:$SAMPLE_ID --rg LB:$SAMPLE_ID --rg PU:$PU --rg PL:ILLUMINA -x $REFBASENAME -U $SAMPLETOMAP$FASTQSUFFIX1 -S $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'.sam'

	fi

	## Convert to bam file for storage
	samtools view -bS -F 4 $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'.sam' > $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'.bam'
	rm $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'.sam'

	## Filter the mapped reads
	# Filter bam files to remove poorly mapped reads (non-unique mappings and mappings with a quality score < 20) -- do we want the quality score filter??
	samtools view -h -q 20 $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'.bam' | samtools view -buS - | samtools sort -o $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'_minq20_sorted.bam'

done
```

Mapping Job
```
#!/bin/bash
#SBATCH --job-name=low_coverage_mapping3.job
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=10G
#SBATCH --mail-user=lucas.jones@uconn.edu
#SBATCH -o mapping3_%j.out
#SBATCH -e mapping3_%j.err

echo `hostname`

bash /labs/Baumann/BaumannLab/SLGenetics/Scripts/low_coverage_mapping.sh \
/labs/Baumann/BaumannLab/SLGenetics/Misc/AdapterClippingLane3List.tsv \
/labs/Baumann/BaumannLab/SLGenetics/Misc/AdapterClippingLane3Table.tsv \
/labs/Baumann/BaumannLab/SLGenetics/adapter_clipped/ \
/labs/Baumann/BaumannLab/SLGenetics/ \
_adapter_clipped_f_paired.fastq.gz \
_adapter_clipped_r_paired.fastq.gz \
very-sensitive \
/labs/Baumann/BaumannLab/SLGenetics/genome/CA/final.genome.scf.fasta \
ADubius \
```

Merge Using Samtools and Sample list - paths to files (same samples over multiple lanes to be added together)

Deduplicate Shell

```
#!/bin/bash
module load bamutil/1.0.7
module load samtools/1.9
module load picard/2.2.1

## This script is used to deduplicate bam files and clipped overlapping read pairs for paired end data. It can process both paired end and single end data.
BAMLIST=$1 # Path to a list of merged, deduplicated, and overlap clipped bam files. Full paths should be included. 
SAMPLETABLE=$2 # Path to a sample table where the 1st column is the prefix of the MERGED bam files. The 4th column is the sample ID, the 2nd column is the lane number, and the 3rd column is sequence ID. The 5th column is population name and 6th column is the data type.
BASEDIR=$3 # Path to the base directory where adapter clipped fastq file are stored in a subdirectory titled "adapter_clipped" and into which output files will be written to separate subdirectories.
REFNAME=$4 # Reference name to add to output files, 

## Loop over each sample
for SAMPLEBAM in `cat $BAMLIST`; do

	## Extract the file name prefix for this sample
	SAMPLEPREFIX=`echo $SAMPLEBAM | sed 's/_bt2_.*//' | sed -e 's#.*/bam/\(\)#\1#'`

	## Remove duplicates and print dupstat file
	# We used to be able to just specify picard.jar on the CBSU server, but now we need to specify the path and version
        java -jar /isg/shared/apps/picard/picard-tools-2.2.1/picard.jar MarkDuplicates I=$SAMPLEBAM O=$BASEDIR'bam/'$SAMPLEPREFIX'_bt2_'$REFNAME'_minq20_sorted_dedup.bam' M=$BASEDIR'bam/'$SAMPLEPREFIX'_bt2_'$REFNAME'_minq20_sorted_dupstat.txt' VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true

	## Extract data type from the merged sample table
	DATATYPE=`grep -P "${SAMPLEPREFIX}\t" $SAMPLETABLE | cut -f 6`

	if [ $DATATYPE != se ]; then
		## Clip overlapping paired end reads (only necessary for paired end data)
		/isg/shared/apps/bamutil/bamUtil_1.0.7/bamUtil/bin/bam clipOverlap --in $BASEDIR'bam/'$SAMPLEPREFIX'_bt2_'$REFNAME'_minq20_sorted_dedup.bam' --out $BASEDIR'bam/'$SAMPLEPREFIX'_bt2_'$REFNAME'_minq20_sorted_dedup_overlapclipped.bam' --stats
	fi

done
```
Deduplicate Job

```
#!/bin/bash
#SBATCH --job-name=deduplicate.job
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=32G
#SBATCH --mail-user=lucas.jones@uconn.edu
#SBATCH -o deduplicate_%j.out
#SBATCH -e deduplicate_%j.err

bash /labs/Baumann/BaumannLab/SLGenetics/Scripts/deduplicate_clipoverlap.sh \
/labs/Baumann/BaumannLab/SLGenetics/Misc/bam_list_merged.tsv \
/labs/Baumann/BaumannLab/SLGenetics/Misc/sample_table_merged.tsv \
/labs/Baumann/BaumannLab/SLGenetics/ \
ADubius \
```

ANGSD SNP calling Shell
```
#!/bin/bash
## This script is used to call SNPs using angsd

BAMLIST=$1 # Path to textfile listing bamfiles to include in global SNP calling with absolute paths
BASEDIR=$2 # Path to the base directory where adapter clipped fastq file are stored in a subdirectory titled "adapter_clipped" and into which output files will be written to separate subdirectories. 
REFERENCE=$3 # Path to reference genome
MINDP=$4 # Minimum depth filter
MAXDP=$5 # Maximum depth filter
MININD=$6 # Minimum individual filter
MINQ=$7 # Minimum quality filter
MINMAF=$8 #Minimum minor allele frequency filter


module load samtools/1.9
module load angsd/0.933-102

## Extract the name of the bam list (excluding path and suffix)
BAMLISTNAME=`echo $BAMLIST | sed 's/\..*//' | sed -e 's#.*/\(\)#\1#'`

## Build base name of output files
OUTBASE=$BAMLISTNAME'_mindp'$MINDP'_maxdp'$MAXDP'_minind'$MININD'_minq'$MINQ

## Call SNPs
/isg/shared/apps/angsd/0.933-102/angsd -b $BAMLIST -anc $REFERENCE -out $BASEDIR'angsd/'$OUTBASE -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doPost 1 -doBcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -doIBS 1 -makematrix 1 -doCov 1 -P 10 -SNP_pval 1e-6 -setMinDepth $MINDP -setMaxDepth $MAXDP -minInd $MININD -minQ $MINQ -minMaf $MINMAF >& $BASEDIR'angsd/'$OUTBASE'.log'

## Create a SNP list to use in downstream analyses
gunzip -c $BASEDIR'angsd/'$OUTBASE'.mafs.gz' | cut -f 1,2,3,4 | tail -n +2 > $BASEDIR'angsd/global_snp_list_'$OUTBASE'.txt'
/isg/shared/apps/angsd/0.933-102/angsd sites index $BASEDIR'angsd/global_snp_list_'$OUTBASE'.txt'

## Also make it in regions format for downstream analyses
cut -f 1,2 $BASEDIR'angsd/global_snp_list_'$OUTBASE'.txt' | sed 's/\t/:/g' > $BASEDIR'angsd/global_snp_list_'$OUTBASE'.regions'

## Lastly, extract a list of chromosomes/LGs/scaffolds for downstream analysis
cut -f1 $BASEDIR'angsd/global_snp_list_'$OUTBASE'.txt' | sort | uniq > $BASEDIR'angsd/global_snp_list_'$OUTBASE'.chrs'
```
ANGSD SNP calling Job
```
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 24
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=20G
#SBATCH --mail-user=lucas.jones@uconn.edu
#SBATCH -o angsd_global_snp_calling%j.out
#SBATCH -e angsd_global_snp_calling%j.err


bash /labs/Baumann/BaumannLab/SLGenetics/Scripts/angsd_global_snp_calling.sh \
/labs/Baumann/BaumannLab/SLGenetics/bam/bam_list_dedup_overlapclipped_final.txt \
/labs/Baumann/BaumannLab/SLGenetics/ \
/labs/Baumann/BaumannLab/SLGenetics/genome/CA/final.genome.scf.fasta \
93 \
370 \
77 \
20 \
0.01 \
```

PCAngsd Shell

```
#!/bin/bash
## This script is used to run PCA using pcangsd. It can be used to run individual-based PCA, estimate selection, inbreeding coefficient, kinship, admixture, and others. The input is a beagle formatted genotype likelihood file.
## See https://github.com/Rosemeis/pcangsd for details

module load pcangsd/0.99
module load python/2.7.9
module load python/3.8.1

BASEDIR=$1 # Path to the base directory where adapter clipped fastq file are stored in a subdirectory titled "adapter_clipped" and into which output files will be written to separate subdirectories. An example for the Greenland cod data is: /workdir/cod/greenland-cod/
BEAGLE=$2 # Path to the beagle formatted genotype likelihood file
MINMAF=$3 # Minimum allele frequency filter
ANALYSIS=$4 # Type of analysis with pcangsd. It can be one of the following: pca, selection, inbreedSites, kinship, admix
MINE=$5 # Minimum number of eigenvectors to use in the modelling of individual allele frequencies (relevant for admix)
MAXE=$6 # Maximum number of eigenvectors to use in the modelling of individual allele frequencies (relevant for admix)

PREFIX=`echo $BEAGLE | sed 's/\..*//' | sed -e 's#.*/\(\)#\1#'`

if [ $ANALYSIS = pca ]; then
	python3 /isg/shared/apps/pcangsd/0.99/pcangsd.py -beagle $BEAGLE -minMaf $MINMAF -threads 16 -o $BASEDIR'angsd/pcangsd_'$PREFIX

elif [ $ANALYSIS = selection ]; then
	python3 /isg/shared/apps/pcangsd/0.99/pcangsd.py -beagle $BEAGLE -selection -minMaf $MINMAF -threads 16 -o $BASEDIR'angsd/pcangsd_'$PREFIX -sites_save

elif [ $ANALYSIS = inbreedSites ]; then
	python3 /isg/shared/apps/pcangsd/0.99/pcangsd.py -beagle $BEAGLE -inbreedSites -minMaf $MINMAF -threads 16 -o $BASEDIR'angsd/pcangsd_'$PREFIX -sites_save

elif [ $ANALYSIS = kinship ]; then
	python3 /isg/shared/apps/pcangsd/0.99/pcangsd.py -beagle $BEAGLE -kinship -minMaf $MINMAF -threads 16 -o $BASEDIR'angsd/pcangsd_'$PREFIX

elif [ $ANALYSIS = admix ]; then
	for ((E = $MINE; E <= $MAXE; E++)); do
		echo $E
		python3 /isg/shared/apps/pcangsd/0.99/pcangsd.py -beagle $BEAGLE -admix -e $E -minMaf $MINMAF -threads 16 -o $BASEDIR'angsd/pcangsd_'$PREFIX'_e'$E
	done

fi
```
PCAngsd Job

```
#!/bin/bash
#SBATCH --job-name=run_pcangsd.job
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=100G
#SBATCH --mail-user=lucas.jones@uconn.edu
#SBATCH -o run_pcangsd_%j.out
#SBATCH -e run_pcangsd_%j.err

bash /labs/Baumann/BaumannLab/SLGenetics/analysis/scripts/run_pcangsd.sh \
/labs/Baumann/BaumannLab/SLGenetics/ \
/labs/Baumann/BaumannLab/SLGenetics/angsd/trimmed_bam_list_dedup_overlapclipped_final_mindp93_maxdp370_minind77_minq20.beagle.gz \
0.05 \
pca \
1 \
6 \
```

Plot in R
```
#### PCA containing outliers ####
#Note because of formatting with cbind V2 correlates to PC1, V3 correlates to PC2, etc
cov <- as.matrix(read.table("~/Genomics/pcangsd_trimmed_bam_list_dedup_overlapclipped_final_mindp93_maxdp370_minind77_minq20.cov", header = F))
pop <- c(pca_sample_table$V2)
mme.pca <- eigen(cov) #perform the pca using the eigen function. 
eigenvectors = mme.pca$vectors #extract eigenvectors 
pca.vectors = as.data.frame(cbind(pop, eigenvectors)) #combine with our population assignments
df = type_convert(pca.vectors)
pca = ggplot(data = df, aes(x=V2, y=V3)) + geom_point()+
geom_point(aes(color=pop)) #(shape=pop)) #alpha=0.25)
pca_color <- pca + scale_color_manual("Population", values=c("violet","purple", "lightsalmon1", "plum4", "midnightblue", "grey1", "royalblue", "red", "skyblue","firebrick4", "wheat4", "limegreen", "snow4", "black"))+
  xlab("PC1") + ylab("PC2") 
ggsave(filename = "~/Genomics/pcangsd_trimmed_pca.pdf", plot = pca_color)

#### PCA removing outliers ####
#Note because of formatting with cbind V2 correlates to PC1, V3 correlates to PC2, etc
nocov <- as.matrix(read.table("~/Genomics/pcangsd_trimmed_bam_list_dedup_overlapclipped_final_mindp93_maxdp370_minind77_minq20.cov", header = F))[-c(256:279),-c(256:279)]
mme.pca <- eigen(cov)#perform the pca using the eigen function. 
eigenvectors = mme.pca$vectors #extract eigenvectors 
eigenvectors1 <- eigenvectors
pop1 <- c(pca_sample_table_final$V2)
pca.vectors1 = as.data.frame(cbind(pop1, eigenvectors1)) #combine with our population assignments
df1 = type_convert(pca.vectors1)
pca1 = ggplot(data = df1, aes(x=V2, y=V3))+ geom_point()+
   geom_point(aes(color=pop1)) #alpha=0.25)
pca_color1 <- pca1 + scale_color_manual("Population", values=c("violet", "lightsalmon1","plum4", "grey1", "royalblue", "red", "skyblue", "firebrick4", "wheat4", "limegreen"))+
  xlab("PC1") + ylab("PC2") 
ggsave(filename = "~/Genomics/pcangsd_trimmed_pca1.pdf", plot = pca_color1)

pca2 = ggplot(data = df1, aes(x=V3, y=V5))+ geom_point()+
  geom_point(aes(color=pop1)) #alpha=0.25)
pca_color2 <- pca2 + scale_color_manual("Population", values=c("violet", "lightsalmon1","plum4", "grey1", "royalblue", "red", "skyblue", "firebrick4", "wheat4", "limegreen"))+
  xlab("PC2") + ylab("PC4") 
ggsave(filename = "~/Genomics/pcangsd_trimmed_pca2.pdf", plot = pca_color2)

pca3 = ggplot(data = df1, aes(x=V3, y=V4))+ geom_point()+
  geom_point(aes(color=pop1)) #alpha=0.25)
pca_color3 <- pca3 + scale_color_manual("Population", values=c("violet", "lightsalmon1","plum4", "grey1", "royalblue", "red", "skyblue", "firebrick4", "wheat4", "limegreen"))+
  xlab("PC2") + ylab("PC3") 
ggsave(filename = "~/Genomics/pcangsd_trimmed_pca3.pdf", plot = pca_color3)
```
![image](https://user-images.githubusercontent.com/65695212/117350665-d48a2b00-ae7a-11eb-8016-b04af148ddb4.png)

![image](https://user-images.githubusercontent.com/65695212/117350588-bb817a00-ae7a-11eb-9043-e00ce51756b8.png)
![image](https://user-images.githubusercontent.com/65695212/117350686-d94edf00-ae7a-11eb-8536-e71e6705639b.png)


![image](https://user-images.githubusercontent.com/65695212/117350588-bb817a00-ae7a-11eb-9043-e00ce51756b8.png)
