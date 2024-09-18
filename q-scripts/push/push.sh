#!/bin/bash

#q    require $EXPERIMENTS_DIR $EXPERIMENT $PUSH_SERVER $PUSH_DIR $PUSH_USER $PUSH_KEY 

#$    -N  push_$EXPERIMENT
#$    -wd $EXPERIMENTS_DIR
#$    -l  vf=4G
#$    -q  wilsonte_lab.q,all.q

#SBATCH --job-name=push_$EXPERIMENT
#SBATCH --chdir=$EXPERIMENTS_DIR
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=4G 
#SBATCH --time=12:00:00
#SBATCH --account=$SLURM_ACCOUNT

SRC_EXPERIMENT_DIR=${EXPERIMENTS_DIR}/${EXPERIMENT}
DEST_EXPERIMENT_DIR=${PUSH_DIR}/experiments/${EXPERIMENT}
cd ${SRC_EXPERIMENT_DIR} # since this script is often called with --execute

ACCESS="-i ${PUSH_KEY} -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null"
SERVER="${PUSH_USER}@${PUSH_SERVER}"
SSH="ssh ${ACCESS} ${SERVER}"
SCP="scp ${ACCESS}"

make_sub_dir () {
    SUB_DIR=$1
    $SSH mkdir -p ${DEST_EXPERIMENT_DIR}/${SUB_DIR}
}
push_files () {
    SUB_DIR=$1
    SOURCE_GLOB=$2
    RECURSIVE=$3
    echo $PWD
    $SCP ${RECURSIVE} ${SUB_DIR}/${SOURCE_GLOB} ${SERVER}:${DEST_EXPERIMENT_DIR}/${SUB_DIR}
    checkPipe
}

echo "pushing files to web server instance for visualization"

echo "  compare"
make_sub_dir compare
echo "    txt"
push_files compare msvtools.compare.samples*.txt
echo "    bed"
push_files compare msvtools.compare.SVs*.bed

echo "  samples"
make_sub_dir samples
echo "    txt"
push_files samples msvtools.segment.probes.*.bgz*
echo "    plots"
push_files samples plots -r

echo "  summary"
$SCP summary_images.html ${SERVER}:${DEST_EXPERIMENT_DIR}

echo "setting permission on folder"
$SSH sudo chgrp -R www-data ${DEST_EXPERIMENT_DIR}

echo "done"
