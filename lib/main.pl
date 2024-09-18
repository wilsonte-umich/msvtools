use strict;
use warnings;
use Cwd(qw(abs_path));
			
#========================================================================
# 'main.pl' is the command-line interpreter and help generator
#========================================================================

#========================================================================
# define variables
#------------------------------------------------------------------------
# global program information
#------------------------------------------------------------------------
use vars qw($utility $libDir);
require "$libDir/command.pl";
require "$libDir/file.pl";
require "$libDir/formats.pl";
require "$libDir/common.pl";
our $error = "error";
our ($command, @args) = @ARGV;
$| = 1;
#------------------------------------------------------------------------
# commands and options
#------------------------------------------------------------------------
our %options; # collects the actual options specified by user on command line
#------------------------------------------------------------------------
my @mdlName = ("d", "<str>", "name for the copy number model (e.g. a cell line)");
my @trnName = ("t", "<str>", "name of the model used to train this analysis [model-name]");
my @refFile = ("r", "<str>", "full path to 'refine.probes' file used to train this analysis");
my @input   = ("i", "<str>", "input directory containing input files (e.g. from previous command)");
my @file    = ("f", "<str>", "path to the input file(s) [STDIN]");
my @output  = ("o", "<str>", "output directory where files will be placed (must exist)");
my @outFile = ("o", "<str>", "output file path");

my @format  = ("F", "<str>", "array format, i.e. vendor (Illumina) [Illumina]");
my @project = ("p", "<str>", "project, name, exclusive of dates, e.g. Prjt_273");
my @nameCol = ("N", "<str>", "name of column to use as sample names [Illumina=>DNA_ID]");
my @sample  = ("s", "<str>", "sample name for the input data");
my @tmpDir  = ("T", "<str>", "temporary directory (must exist) [/tmp]");
my @maxMem  = ("m", "<int>", "maximum RAM bytes to use when sorting [1000000000]");

my @exclude = ("x", "<str>", "path to a BED file containing genomic regions to exclude");
my @skipChrom=("k", "<str>", "comma-delimited list of chromosomes to skip in SV calling [none]");

my @library = ("l", "<str>", "name of sequencing library (a sample may have many libraries)");
my @group   = ("G", "<str>", "group name for set of input samples/libraries");

my @plotDir = ("P", "<str>", "target directory for summary plots (must exist) [plots not created]");
my @repeat  = ("r", undef,   "repeat all input data to STDOUT (requires --output/-o)");
my @unique  = ("u", undef,   "only report variants called in exactly one sample");
my @strict  = ("U", undef,   "only report variants with matching pairs in one sample");
my @minSize = ("z", "<int>", "only report variants >= z bp on the same chromosome [0]");
my @maxSize = ("Z", "<int>", "only report variants <= Z bp on the same chromosome [1000000000]");
my @chrom   = ("c", "<str>", "chromosome to analyze in this execution");
my @chroms  = ("C", "<str>", "comma-delimited list of all chromosomes to be analyzed");
#------------------------------------------------------------------------
our %commands =  ( # 0=allowed, 1=required
#------------------------------------------------------------------------
# list the samples in an array project/run
#------------------------------------------------------------------------
    list => [
        \&msvtools_list,
        "list the samples in an array project/run, in project order", {
		'input-dir'=>       [1, 1,  @input, "Input Options"],
        'array-format' =>   [0, 2,  @format],     
        'project' =>        [1, 3,  @project],
        'name-column' =>    [0, 4,  @nameCol],
    }],
#    index => [
#        \&msvtools_index,
#        "reformat sample data files and index for random access/browser display", {
#		'input-dir'=>       [1, 1,  @input, "Input Options"],
#        'array-format' =>   [0, 2,  @format],     
#        'project' =>        [1, 3,  @project],
#        'name-column' =>    [0, 4,  @nameCol],
#        'sample' =>         [1, 5,  @sample],
#        'output-dir'=>      [1, 11, @output, "Output Options"],
#        'max-mem'=>         [0, 13, @maxMem],
#        #'plot-dir'=>        [0, 19, @plotDir],
#    }],
#------------------------------------------------------------------------
# create a baseline copy number/zygosity model for all genome regions
#------------------------------------------------------------------------
    train => [
        \&msvtools_train,
        "create a baseline copy number/zygosity model for a genome", {
        'method' =>         [1, 1,  "M", "<str>", "method used to construct model (bam|bed6|M|F)", "Input Options"],
        'max-mem'=>         [0, 2,  @maxMem],
        'tmp-dir'=>         [0, 3,  @tmpDir],
        'array-dir'=>       [1, 10, "a", "<str>", "directory containing array files", "Array Options"], 
        'array-format' =>   [0, 11, @format],
        'name-column' =>    [0, 13, @nameCol],
        'array-type'=>      [1, 14, "A", "<str>", "text descriptor of the array probe set"],
        'genome-fasta'=>    [1, 15, "G", "<str>", "path to the reference genome fasta matching the array"],
        'bin-size'=>        [0, 21, "z", "<int>", "use this bin size for paired-end bam method [5000]", "BAM Options"],
        'max-TLen'=>        [0, 22, "L", "<int>", "max allowed insert size for paired-end bam method [1000]"],
        'min-qual'=>        [0, 23, "q", "<int>", "minimum MAPQ of one read for a pair to be considered [5]"],
		'reference'=>       [0, 24, "r", "<str>", "copy number reference region as chrom:start-end/copy_number"],
		'max-copy-number'=> [0, 25, "c", "<int>", "maximum copy number allowed in the HMM [4]"],
        'chromosomes'=>     [1, 31, @chroms, "Output Options"],     
        'output-dir'=>      [1, 32, @output],
        'model-name'=>      [1, 33, @mdlName],
    }],
    refine => [ # NB: can repeat refine sequentially to establish convergence
        \&msvtools_refine,
        "refine the copy number/zygosity model using a set of array samples", {
		'model-dir'=>       [1, 1,  "D", "<str>", "directory containing training model files", "Model Options"], 
        'model-name'=>      [1, 2,  @mdlName],
        'train-name'=>      [0, 3,  @trnName],
        'refine-file'=>     [0, 4,  @refFile],
        'array-dir'=>       [1, 10, "a", "<str>", "directory containing array files", "Array Options"], 
        'array-format' =>   [0, 11, @format],     
        'project' =>        [1, 12, @project],
        'name-column' =>    [0, 13, @nameCol],
        'array-type'=>      [1, 14, "A", "<str>", "text descriptor of the array probe set"],
        'samples' =>        [1, 15, "s", "<chr>", "comma-delimited list of samples to use from project"],
        'chromosomes'=>     [1, 20, @chroms, "Output Options"],
        'output-dir'=>      [1, 21, @output],
        'max-copy-number'=> [0, 23, "c", "<int>", "maximum copy number allowed in the HMM [4]"],
        'transition-prob'=> [0, 24, "t", "<dbl>", "HMM transition probability [1e-4]"],
        'tmp-dir'=>         [0, 30,  @tmpDir],
        'max-mem'=>         [0, 33, @maxMem],
    }],
#------------------------------------------------------------------------
# use the baseline copy number/zygosity model to call SV change states per sample
#------------------------------------------------------------------------	
    segment => [
        \&msvtools_segment,
        "find genome spans different than the baseline model, per sample", {
        'refine-file'=>     [1, 1,  @refFile, "Model Options"],
        'array-dir'=>       [1, 10, "a", "<str>", "directory containing array files", "Array Options"], 
        'array-format' =>   [0, 11, @format],     
        'project' =>        [1, 12, @project],
        'name-column' =>    [0, 13, @nameCol],
        'sample' =>         [1, 14, @sample],
        'extreme-zyg'=>     [0, 21, "z", "<dbl>", "largest zygosity informative to copy number [0.85]", "Segmentation Options"],  
        'lrr-lim'=>         [0, 22, "l", "<str>", "range of LRR values used for model fitting [-1.5,1.0]"],
        'bad-probe-freq'=>  [0, 23, "b", "<dbl>", "assume this fraction of probes give unpredictable values [1e-3]"],
        'transition-prob'=> [0, 24, "t", "<dbl>", "HMM transition probability [1e-4]"],
        'preservation'=>    [0, 25, "P", "<str>", "penalities for CN changes from model [0.99,0.95,0.9,0.8]"],
        'output-dir'=>      [1, 31, @output, "Output Options"],
        'tmp-dir'=>         [0, 40, @tmpDir],
        'max-mem'=>         [0, 43, @maxMem],
    }],
#------------------------------------------------------------------------
# use a set SV change states for many samples to determine uniqueness
#------------------------------------------------------------------------	
    compare => [
        \&msvtools_compare,
        "compare found SVs/CNVs among many samples to determine uniqueness", {
        'input-dir'=>       [1, 1,  @input, "Input Options"],
        'samples' =>        [1, 2, "s", "<chr>", "comma-delimited list of samples to compare"],
        'output-dir'=>      [1, 11, @output, "Output Options"],
        'size-threshold'=>  [0, 12, "z", "<int>", "make two files of large and small SVs, above and below z bp [2000000]"],
        'group-name'=>      [1, 13, @group],
        'tmp-dir'=>         [0, 14, @tmpDir],
    }],
    
#------------------------------------------------------------------------
# compare 
#------------------------------------------------------------------------	 
    
#------------------------------------------------------------------------
# option handling
#------------------------------------------------------------------------	
);
our %longOptions; # for converting short options to long
foreach my $command(keys %commands){
    %{$longOptions{$command}} =
        map { $commands{$command}[2]{$_}[2] => $_ } keys %{$commands{$command}[2]}; 
}
#------------------------------------------------------------------------
# break utility commands into managable help chunks 
#------------------------------------------------------------------------
sub reportCommandChunks {
    reportCommandChunk('Project Utilities', qw(list));    
    reportCommandChunk('Create Genome Model', qw(train refine));
    reportCommandChunk('Find SVs/CNVs', qw(segment));
    reportCommandChunk('Compare SVs/CNVs', qw(compare));
}
#========================================================================

1;
