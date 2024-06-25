#!/bin/bash

#q    require $EXPERIMENT $CELL_LINE $CELL_CLONE $SAMPLES
#q    require $EXP_SAMPLES_DIR $EXP_COMPARE_DIR

#$    -N  compare_$EXPERIMENT\_$CELL_CLONE
#$    -wd $EXP_COMPARE_DIR
#$    -l  vf=4G

#$    -q  wilsonte_lab.q
##,all.q

echo "comparing called SVs across samples"
echo "EXPERIMENT = $EXPERIMENT"
echo "CELL_LINE  = $CELL_LINE"
echo "CELL_CLONE = $CELL_CLONE"
echo "SAMPLES = "
echo $SAMPLES
echo

SAMPLES=`echo "$SAMPLES" | sed 's/ /,/g'`

msvtools compare \
-i $EXP_SAMPLES_DIR \
-s $SAMPLES \
-o $EXP_COMPARE_DIR \
-G $CELL_CLONE
checkPipe

echo
echo "done"

#compare     compare found SVs/CNVs among many samples to determine uniqueness
#
#Input Options
#    -i,--input-dir <str>     **REQUIRED** input directory containing input files (e.g. from previous command)
#    -s,--samples <chr>       **REQUIRED** comma-delimited list of samples to compare
#
#Output Options
#    -o,--output-dir <str>    **REQUIRED** output directory where files will be placed (must exist)
#    -z,--size-threshold <int>make two files of large and small SVs, above and below z bp [2000000]
#    -G,--group-name <str>    **REQUIRED** group name for set of input samples/libraries
#    -T,--tmp-dir <str>       temporary directory (must exist) [/tmp]
