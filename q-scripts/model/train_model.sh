#!/bin/bash

#q    require $CELL_LINE $METHOD $MODELS_DIR
#q    require $BED_FILE $BAM_GLOB $REFERENCE $CHROMS
#q    require $CORE_PRJT_DIR $ARRAY_FORMAT $NAME_COLUMN $ARRAY_TYPE $GENOME_FASTA

#$    -N  train_$CELL_LINE
#$    -wd $MODELS_DIR
#$    -l  vf=4G
#$    -q  wilsonte_lab.q

#SBATCH --job-name=train_$CELL_LINE
#SBATCH --chdir=$MODELS_DIR
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=8G 
#SBATCH --time=24:00:00
#SBATCH --account=$SLURM_ACCOUNT

echo "training $CELL_LINE model using method $METHOD"
echo $BED_FILE
echo $BAM_GLOB
echo $REFERENCE
echo

# set the data file based on method
DATA_FILE=$BED_FILE
if [ "$METHOD" = "bam" ]; then
    DATA_FILE=$BAM_GLOB
fi

msvtools train \
-M $METHOD \
-a $CORE_PRJT_DIR \
-F $ARRAY_FORMAT \
-N $NAME_COLUMN \
-A $ARRAY_TYPE \
-G $GENOME_FASTA \
-z 5000 \
-L 1000 \
-q 5 \
-r $REFERENCE \
-c 5 \
-C $CHROMS \
-o $MODELS_DIR \
-d $CELL_LINE $DATA_FILE 
checkPipe

echo
echo "done"


# usage:  msvtools <command> [options] <input file(s)> [chr:start-end]

# train       create a baseline copy number/zygosity model for a genome

# Input Options
#     -M,--method <str>        **REQUIRED** method used to construct model (bam|M|F)
#     -m,--max-mem <int>       maximum RAM bytes to use when sorting [1000000000]
#     -T,--tmp-dir <str>       temporary directory (must exist) [/tmp]

# Array Options
#     -a,--array-dir <str>     **REQUIRED** directory containing array files
#     -F,--array-format <str>  array format, i.e. vendor (Illumina) [Illumina]
#     -N,--name-column <str>   name of column to use as sample names [Illumina=>DNA_ID]
#     -A,--array-type <str>    **REQUIRED** text descriptor of the array probe set
#     -G,--genome-fasta <str>  **REQUIRED** path to the reference genome fasta matching the array

# BAM Options
#     -z,--bin-size <int>      use this bin size for paired-end bam method [5000]
#     -L,--max-TLen <int>      max allowed insert size for paired-end bam method [1000]
#     -q,--min-qual <int>      minimum MAPQ of one read for a pair to be considered [5]
#     -r,--reference <str>     copy number reference region as chrom:start-end/copy_number
#     -c,--max-copy-number <int>maximum copy number allowed in the HMM [4]

# Output Options
#     -C,--chromosomes <str>   **REQUIRED** comma-delimited list of all chromosomes to be analyzed
#     -o,--output-dir <str>    **REQUIRED** output directory where files will be placed (must exist)
#     -d,--model-name <str>    **REQUIRED** name for the copy number model (e.g. a cell line)
