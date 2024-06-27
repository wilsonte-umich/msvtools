
The following is an example of a valid `sample_table.csv` file
used to describe a set of samples for integrated analysis.

```csv
Experiment,Cell_Line,Cell_Clone,Control,Sample_ID,Sample_Name
exp1,cellLineA,A_WT,Yes,S01,WT_ctrl_1
exp1,cellLineA,A_WT,Yes,S02,WT_ctrl_2
exp1,cellLineA,A_WT,No,S03,WT_treated_1
exp1,cellLineA,A_WT,No,S04,WT_treated_2
exp1,cellLineA,A_MUT,Yes,S05,MUT_ctrl_1
exp1,cellLineA,A_MUT,Yes,S06,MUT_ctrl_2
exp1,cellLineA,A_MUT,No,S07,MUT_treated_1
exp1,cellLineA,A_MUT,No,S08,MUT_treated_2
exp2,etc.
```

A sample file generally describes one array run, i.e., a set of 
array samples processed together at the same time.

Each row after the required header 
describes one sample analyzed by one array, e.g., one bead chip.
Different rows are typically multiple independent sub-clones derived from
a single Cell_Clone as defined below.

Columns are:

### Experiment

A short text label identifying a set of samples considered to be 
part of the same overall experiment, to be analyzed together. 
One array run might be used to generate data for one or more experiments.

### Cell_Line

The cell line or other source identifier for the samples in an an experiment.
Typically, one experiment uses just one cell line, but could use more.

### Cell_Clone

The specific parental clone of the Cell_Line that gave rise to a set
of related sub-clonal samples. Typical examples of cell clones are wild-type and a mutant 
derivative that was re-cloned after some genetic manipulation. The purpose
of identifying a parent cell clone for each sample is that samples derived from the same clone 
may contain shared non-mosaic CNVs that are not present in a different parent clone.
`msvtools` seeks to find mosaic, not shared clonal CNVs.

### Control

Either Yes or No, where Yes means that this sample should be used to refine
the reference copy number model of the Cell_Line to update it to match
a specific Cell_Clone. Some control samples must always be present for each Cell_Clone.
Typically, untreated samples less likely to carry mosaic CNVs are declared as controls.

### Sample_ID

The name identifier for a sample as it is found in the Illumina report files.
Specifically, the sample table Sample_ID value must be found in:
- the FinalReport Name column
- the DNAReport DNA_ID column

where column DNA_ID must be declared in your job master file:

```sh
# job-file.q
$NAME_COLUMN    DNA_ID
```

If you have an older/different formatting of the DNAReport files, change DNA_ID
to the appropriate column header name.

### Sample_Name

A name meaningful to the user in the context of the experiment. This values is not
used for parsing array information and does not need to match any array reports, etc.
