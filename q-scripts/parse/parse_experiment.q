
# initialize the experiment directory
$EXPERIMENT_DIR     $EXPERIMENTS_DIR/$EXPERIMENT
$EXP_MODELS_DIR     $EXPERIMENT_DIR/models
$EXP_SAMPLES_DIR    $EXPERIMENT_DIR/samples
$EXP_COMPARE_DIR    $EXPERIMENT_DIR/compare
invoke file/create.q $DIR $EXPERIMENT_DIR $EXP_MODELS_DIR $EXP_SAMPLES_DIR $EXP_COMPARE_DIR

# get the cell line for this experiment (generally expect only one)
$CELL_LINES  run $PARSE_SAMPLES Cell_Line

# cascade to queue all jobs by cell line
invoke $SLAVES_DIR/parse/parse_cell_line.q $CELL_LINE $CELL_LINES

# make the composite viewer for displaying sample jpg summary images
$SAMPLES run $PARSE_SAMPLES Sample_ID
qsub $SLAVES_DIR/assemble/assemble.sh
