#!/bin/bash

#q    require $EXPERIMENT $EXPERIMENT_DIR $CELL_CLONES

#$    -N  assemble_$EXPERIMENT
#$    -wd $EXPERIMENT_DIR
#$    -l  vf=16G
#$    -q  wilsonte_lab.q,all.q

#SBATCH --job-name=assemble_$EXPERIMENT
#SBATCH --chdir=$EXPERIMENT_DIR
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=16G 
#SBATCH --time=24:00:00
#SBATCH --account=$SLURM_ACCOUNT

echo "assembling jpg viewer for all samples in experiment $EXPERIMENT"
echo

# collect the lists of samples and aliases
SAMPLES=""
ALIASES=""
for CELL_CLONE in $CELL_CLONES; do
    SMP_FILE=$EXPERIMENT_DIR/compare/msvtools.compare.samples.$CELL_CLONE.txt
    SMPS=`head -n1 $SMP_FILE`
    ALSS=`tail -n1 $SMP_FILE`
    SAMPLES="$SAMPLES $SMPS"
    ALIASES="$ALIASES $ALSS"
done
export SAMPLES
export ALIASES

HTML_FILE=$EXPERIMENT_DIR/summary_images.html
perl $SLAVES_DIR/assemble/assemble.pl > $HTML_FILE
checkPipe

echo
echo "done"
