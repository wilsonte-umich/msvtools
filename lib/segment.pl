use strict;
use warnings;

# define variables
use vars qw($version $utility $error $libDir
            $TRUE $FALSE
            %prbCol);
require "$libDir/Illumina.pl";
require "$libDir/formats.pl";
our ($refineFile,
     $tmpDir, $maxMem,
     $arrayDir, $arrayFmt, $project, $nameCol, $sample,
     $extZyg, $lrrLim, $badPrbFreq, $transProb, $preservation,
     $outputDir, $plotDir);
our ($inputDir, $valType, %chroms);
my ($mdlH, $datH, %mdl);
my ($nMdl, $nArrIn, $nArrMatch) = (0, 0, 0);

# manage options
sub setOptions_segment {
    setOptionValue(\$refineFile,  'refine-file');
    setOptionValue(\$arrayDir,    'array-dir');
    setOptionValue(\$arrayFmt,    'array-format', 'Illumina');
    setOptionValue(\$project,     'project');
    setOptionValue(\$nameCol,     'name-column',  'DNA_ID');
    setOptionValue(\$sample,      'sample');
    setOptionValue(\$extZyg,      'extreme-zyg',     0.85);
    setOptionValue(\$lrrLim,      'lrr-lim',         '-1.5,1.0');
    setOptionValue(\$badPrbFreq,  'bad-probe-freq',  1e-3);
    setOptionValue(\$transProb,   'transition-prob', 1e-4);
    setOptionValue(\$preservation,'preservation',    '0.99,0.95,0.9,0.8');
    setOptionValue(\$outputDir,   'output-dir');
    setOptionValue(\$tmpDir,      'tmp-dir',      '/tmp');
    setOptionValue(\$maxMem,      'max-mem',         1000000000);  
    #-------------------------------------
    -e $refineFile or die "$error: file not found: $refineFile\n";
    -d $arrayDir or die "$error: directory not found: $arrayDir\n";
    -d $outputDir or die "$error: directory not found: $outputDir\n";
    -d $tmpDir or die "$error: $tmpDir does not exist or is not a directory\n";
    $maxMem =~ m|\D| and die "$error: max-mem must be an integer number of bytes\n";
}
# main execution block
sub msvtools_segment {

    # initialize
    setOptions_segment();
    print STDERR "$utility segment: " . getTime(), "\n";
    
    # # load the trained model CN (and read count) for each probe
    # print STDERR "loading model refinement file\n";
    # loadProbeModel();
    # print STDERR join("\t", commify($nMdl), "probes recovered with data from model refinement file"), "\n";

    # # load probe data for target sample
    # print STDERR "calculating probe median values\n";
    # $inputDir = $arrayDir;
    # my $datFile = getOutFile('data', $sample, 'gz');
    # openOutputStream($datFile, \$datH, $TRUE);
    # print $datH join("\t", qw(
    #     CHROM START POS PRB_NAME IS_NA STRAND
    #     LRR BAF 
    #     AS_IN CN_IN INF_IN
    #     LRR_MED BAF_MED
    # )), "\n";
    # indexIllumina(\&getProbeData); # not sorted yet, occurs in R at HMM
    # closeHandles($datH);
    # print STDERR join("\t", commify($nArrIn),    "probes recovered from array"), "\n";
    # print STDERR join("\t", commify($nArrMatch), "probes matched the model refinement file"), "\n";

###############################
my $datFile = getOutFile('data', $sample, 'gz');

    # execute the HMM via R
    %mdl = ();
    $ENV{SEGMENTING_ACTION} = "SEGMENT_SAMPLE";
    $ENV{GC_BIAS_FILE} = $refineFile;
    $ENV{GC_BIAS_FILE} =~ s/refine.probes/refine.gc_bias/;
    $ENV{GC_BIAS_FILE} =~ s/bed.bgz/rds/;
    segmentArray($sample, $datFile, "$outputDir/plots", "__NULL__",
                 $extZyg, $lrrLim, $badPrbFreq, $transProb, $preservation);

#################
die "REMOVE THIS LINE\n";

}

sub loadProbeModel {
    print STDERR "$utility segment: loading model: $refineFile\n";
    openInputStream($refineFile, \$mdlH, undef, undef, $TRUE);
    my $header = <$mdlH>;
    while (my $prb = <$mdlH>){
        chomp $prb;
        my @f = split("\t", $prb);
        my $prbName = $f[$prbCol{PRB_NAME}];
        my $rejectProbe;
        my @prbData = map { 
            my $val = $f[$prbCol{$_}];
            (!defined($val) or $val eq "NA") and $rejectProbe++;
            $val;
        } qw(AS_OUT CN_OUT INF_OUT 
             LRR BAF);
        $rejectProbe and next;
        $nMdl++;
        $mdl{$prbName} = \@prbData;
        $chroms{$f[$prbCol{CHROM}]}++;
    }
    closeHandles($mdlH);
}

sub getProbeData {
    my ($f, $hdr) = @_;
    $nArrIn++;
    my $prbName = $$f[$$hdr{PRB_NAME}];
    $mdl{$prbName} or return;
    $nArrMatch++;
    print $datH join("\t",
        $$f[$$hdr{CHROM}], $$f[$$hdr{POS}]-1, $$f[$$hdr{POS}], $prbName, -1, '+',
        $$f[$$hdr{LRR}], $$f[$$hdr{BAF}], # this sample's values, could be NA even if other samples have data
        @{$mdl{$prbName}}
    ), "\n";
}

1;
