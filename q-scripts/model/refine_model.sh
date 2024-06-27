#!/bin/bash

#q    require $EXPERIMENT $CELL_LINE $CELL_CLONE $SAMPLES
#q    require $MODELS_DIR
#q    require $CORE_PRJT_DIR $CORE_PRJT $ARRAY_FORMAT $NAME_COLUMN
#q    require $CHROMS $EXP_MODELS_DIR

#$    -N  refine_$EXPERIMENT\_$CELL_CLONE
#$    -wd $EXP_MODELS_DIR
#$    -l  vf=24G
#$    -q  wilsonte_lab.q,all.q

#SBATCH --job-name=refine_$EXPERIMENT\_$CELL_CLONE
#SBATCH --chdir=$EXP_MODELS_DIR
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=24G 
#SBATCH --time=24:00:00
#SBATCH --account=$SLURM_ACCOUNT

echo "refining copy number model using array data"
echo "EXPERIMENT = $EXPERIMENT"
echo "CELL_LINE  = $CELL_LINE"
echo "CELL_CLONE = $CELL_CLONE"
echo "SAMPLES = "
echo $SAMPLES
echo

SAMPLES=`echo "$SAMPLES" | sed 's/ /,/g'`

echo "----------------------------------------------------"
echo "refinement round 1"
echo "----------------------------------------------------"
msvtools refine \
-D $MODELS_DIR \
-t $CELL_LINE \
-a $CORE_PRJT_DIR \
-p $CORE_PRJT \
-F $ARRAY_FORMAT \
-N $NAME_COLUMN \
-s $SAMPLES \
-C $CHROMS \
-o $EXP_MODELS_DIR \
-d $CELL_CLONE \
-c 5
checkPipe
echo "----------------------------------------------------"
echo

echo "----------------------------------------------------"
echo "refinement round 2"
echo "----------------------------------------------------"
msvtools refine \
-r $EXP_MODELS_DIR/msvtools.refine.probes.$CELL_CLONE.bed.bgz \
-a $CORE_PRJT_DIR \
-p $CORE_PRJT \
-F $ARRAY_FORMAT \
-N $NAME_COLUMN \
-s $SAMPLES \
-C $CHROMS \
-o $EXP_MODELS_DIR \
-d $CELL_CLONE"_Round2" \
-c 5
checkPipe
echo "----------------------------------------------------"
echo

echo "done"

#refine      refine the copy number/zygosity model using a set of array samples
#
#Model Options
#    -D,--model-dir <str>     directory containing training model files
#    -d,--model-name <str>    **REQUIRED** name for the copy number model (e.g. a cell line)
#    -t,--train-name <str>    name of the model used to train this analysis [model-name]
#    -r,--refine-file <str>   full path to 'refine.probes' file used to train this analysis
#
#Array Options
#    -a,--array-dir <str>     **REQUIRED** directory containing array files
#    -F,--array-format <str>  array format, i.e. vendor (Illumina) [Illumina]
#    -p,--project <str>       **REQUIRED** project, name, exclusive of dates, e.g. Prjt_273
#    -N,--name-column <str>   name of column to use as sample names [Illumina=>DNA_ID]
#    -s,--samples <chr>       **REQUIRED** comma-delimited list of samples to use from project
#
#Output Options
#    -C,--chromosomes <str>   **REQUIRED** comma-delimited list of all chromosomes to be analyzed
#    -o,--output-dir <str>    **REQUIRED** output directory where files will be placed (must exist)
#    -c,--max-copy-number <int>maximum copy number allowed in the HMM [4]
#    -T,--tmp-dir <str>       temporary directory (must exist) [/tmp]
#    -m,--max-mem <int>       maximum RAM bytes to use when sorting [1000000000]
