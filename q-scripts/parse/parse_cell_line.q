
# get all cell clones for this cell line
$CELL_CLONES run $PARSE_SAMPLES Cell_Clone

# cascade to queue all jobs by cell clone
invoke $SLAVES_DIR/parse/parse_cell_clone.q $CELL_CLONE $CELL_CLONES

preserve $CELL_CLONES $CELL_CLONES
