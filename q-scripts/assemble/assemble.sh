#!/bin/bash

#q    require $EXPERIMENT $EXPERIMENT_DIR $CELL_CLONES

#$    -N  assemble_$EXPERIMENT
#$    -wd $EXPERIMENT_DIR
#$    -l  vf=2G

#$    -q  wilsonte_lab.q
##,all.q

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
