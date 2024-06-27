#!/bin/bash

#q    require $PROJECT $PROJECT_DIR $SAMPLES $INDEXED_DIR

#$    -N  index_$PROJECT
#$    -wd $PROJECT_DIR
#$    -l  vf=2G
#$    -t  1-$N_SAMPLES

#$    -q  wilsonte_lab.q

#SBATCH --job-name=index_$PROJECT
#SBATCH --chdir=$PROJECT_DIR
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=16G 
#SBATCH --time=24:00:00
#SBATCH --account=$SLURM_ACCOUNT

getTaskObject SAMPLE $SAMPLES
echo "indexing $SAMPLE Illumina array data"
echo

msvtools index \
-i $PROJECT_DIR \
-F Illumina \
-p $PROJECT \
-N DNA_ID \
-s $SAMPLE \
-o $INDEXED_DIR \
-m 1000000000
checkPipe

echo
echo "done"
