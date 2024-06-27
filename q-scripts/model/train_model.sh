#!/bin/bash

#q    require $CELL_LINE $METHOD $MODELS_DIR
#q    require $BAM_GLOB $REFERENCE $CHROMS

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
echo $BAM_GLOB
echo $REFERENCE
echo

msvtools train \
-M $METHOD \
-z 5000 \
-L 1000 \
-q 5 \
-r $REFERENCE \
-c 5 \
-C $CHROMS \
-o $MODELS_DIR \
-d $CELL_LINE $BAM_GLOB 
checkPipe

echo
echo "done"
