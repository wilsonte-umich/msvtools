#!/bin/bash

#q    require $EXPERIMENT $CELL_LINE $CELL_CLONE $SAMPLES $N_SAMPLES
#q    require $CORE_PRJT_DIR $CORE_PRJT $ARRAY_FORMAT $NAME_COLUMN
#q    require $EXP_MODELS_DIR $EXP_SAMPLES_DIR

#$    -N  segment_$EXPERIMENT\_$CELL_CLONE
#$    -wd $EXP_SAMPLES_DIR
#$    -l  vf=16G
#$    -t  1-$N_SAMPLES
#$    -q  wilsonte_lab.q
##,all.q

#SBATCH --job-name=segment_$EXPERIMENT\_$CELL_CLONE
#SBATCH --chdir=$EXP_SAMPLES_DIR
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=16G 
#SBATCH --time=24:00:00
#SBATCH --array=1-$N_SAMPLES
#SBATCH --account=$SLURM_ACCOUNT

echo "segmenting array data based on refined copy number model"
getTaskObject SAMPLE $SAMPLES
echo "EXPERIMENT = $EXPERIMENT"
echo "CELL_LINE  = $CELL_LINE"
echo "CELL_CLONE = $CELL_CLONE"
echo "SAMPLE     = $SAMPLE"
echo

msvtools segment \
-r $EXP_MODELS_DIR/msvtools.refine.probes.$CELL_CLONE"_Round2.bed.bgz" \
-a $CORE_PRJT_DIR \
-p $CORE_PRJT \
-F $ARRAY_FORMAT \
-N $NAME_COLUMN \
-s $SAMPLE \
--bad-probe-freq  2e-3 \
--transition-prob 1e-5 \
--preservation    '0.98,0.90,0.5,0.1' \
-o $EXP_SAMPLES_DIR
checkPipe

echo
echo "done"

#segment     find genome spans different than the baseline model, per sample
#
#Model Options
#    -r,--refine-file <str>   **REQUIRED** full path to 'refine.probes' file used to train this analysis
#
#Array Options
#    -a,--array-dir <str>     **REQUIRED** directory containing array files
#    -F,--array-format <str>  array format, i.e. vendor (Illumina) [Illumina]
#    -p,--project <str>       **REQUIRED** project, name, exclusive of dates, e.g. Prjt_273
#    -N,--name-column <str>   name of column to use as sample names [Illumina=>DNA_ID]
#    -s,--sample <str>        **REQUIRED** sample name for the input data
#
#Segmentation Options
#    -z,--extreme-zyg <dbl>   largest zygosity informative to copy number [0.85]
#    -l,--lrr-lim <str>       range of LRR values used for model fitting [-1.5,1.0]
#    -b,--bad-probe-freq <dbl>assume this fraction of probes give unpredictable values [1e-3]
#    -t,--transition-prob <dbl>HMM transition probability [1e-4]
#    -P,--preservation <str>  penalities for CN changes from model [0.99,0.95,0.9,0.8]
#
#Output Options
#    -o,--output-dir <str>    **REQUIRED** output directory where files will be placed (must exist)
#    -T,--tmp-dir <str>       temporary directory (must exist) [/tmp]
#    -m,--max-mem <int>       maximum RAM bytes to use when sorting [1000000000]

