
Folder 'msvtools/q-scripts' carries q master scripts 
for running msvtools analysis using the 
[q-pipeline-manager](https://github.com/wilsonte-umich/q-pipeline-manager).

The q utility assists with integrated pipeline execution and submission
to a job scheduler on a shared cluster server. 
It is not necessary to use q to run msvtools, but it is a convenient 
way of organizing and executing work.

These q worker scripts are best used by making calls to:
- model/train_model.sh (once per cell line)
- parse/parse_experiment.q (once per experimental sample set)

Experiment parsing to job submission is achieved by reading a
`sample_table.csv` file that defines a set of experimental samples.
See `parse/README.md` for details on constructing a proper sample table.

Then, include the following lines in the indicated q master scripts:

```sh
# environment.q
$PARSE_SAMPLES  perl /path/to/msvtools/q-scripts/parse/parse_sample_table.pl
```

```sh
# job-file.q
$EXPERIMENTS  run $PARSE_SAMPLES Experiment
invoke $SLAVES_DIR/parse/parse_experiment.q $EXPERIMENT $EXPERIMENTS
```

See the q documentation and q-scripts for additional details.
