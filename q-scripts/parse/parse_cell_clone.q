
# get reference samples for this cell clone
$SAMPLES run $PARSE_SAMPLES Sample_ID CONTROL
dieIf run [ "$SAMPLES" = "" ] && echo "No Control Samples!"

# create the copy number model from reference arrays for this clone
##qsub $SLAVES_DIR/model/train_model.sh # NA, use a prior model for cell line
qsub $SLAVES_DIR/model/refine_model.sh

# get all samples for this cell clone
$SAMPLES run $PARSE_SAMPLES Sample_ID
$N_SAMPLES run echo $SAMPLES | wc -w

# segment each sample
qsub $SLAVES_DIR/segment/segment_sample.sh

# compare samples of a cell clone to find unique CNVs
qsub $SLAVES_DIR/compare/compare_samples.sh
