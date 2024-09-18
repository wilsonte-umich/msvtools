use strict;
use warnings;

# define variables
use vars qw($version $utility $error $libDir
            $TRUE $FALSE
            %prbCol);
require "$libDir/exclude_subs.pl";
require "$libDir/Illumina.pl";
require "$libDir/formats.pl";
our ($modelDir, $modelName, $trainName, $refineFile,
     $tmpDir, $maxMem,
     $arrayDir, $arrayFmt, $project, $nameCol, $arrayType, $samples,
     $chroms, $outputDir, $plotDir, $maxCN,
     $transProb);
our ($inputDir, $sample, $valType);
our (@samples, @chroms, %chroms, $minQual);
my ($mdlH, $datH);
our (%data, %coord, %fracGC);
my ($nPrbTotal, $nPrbTrained, $nPrbUsed) = (0, 0, 0);

# manage options
sub setOptions_refine {
    setOptionValue(\$modelDir,   'model-dir');
    setOptionValue(\$modelName,  'model-name');
    setOptionValue(\$trainName,  'train-name');
    setOptionValue(\$refineFile, 'refine-file'); # an input training file, not an output
    setOptionValue(\$arrayDir,   'array-dir');
    setOptionValue(\$arrayFmt,   'array-format', 'Illumina');
    setOptionValue(\$project,    'project');
    setOptionValue(\$nameCol,    'name-column',  'DNA_ID');
    setOptionValue(\$arrayType,  'array-type');
    setOptionValue(\$samples,    'samples');
    setOptionValue(\$chroms,     'chromosomes');
    setOptionValue(\$outputDir,  'output-dir');
    setOptionValue(\$maxCN,      'max-copy-number', 4);
    setOptionValue(\$tmpDir,     'tmp-dir',      '/tmp');
    setOptionValue(\$maxMem,     'max-mem',         1000000000); 
    setOptionValue(\$transProb,  'transition-prob', 1e-4); 
    #-------------------------------------
    -d $modelDir or die "$error: directory not found: $modelDir\n";
    $refineFile and (-e $refineFile or die "$error: file not found: $refineFile\n");
    $refineFile or $trainName or die "either train-name or refine-file is required\n";
    $trainName or $trainName = $modelName;
    -d $arrayDir or die "$error: directory not found: $arrayDir\n";
    @samples = sort {$a cmp $b} split(",", $samples);
    @chroms = sort {$a cmp $b} split(",", $chroms);
    %chroms = map { $_ => 1 } @chroms;    
    -d $outputDir or die "$error: directory not found: $outputDir\n";    
    $maxCN =~ m|\D| and die "$error: max-mem must be an integer number\n";    
    -d $tmpDir or die "$error: $tmpDir does not exist or is not a directory\n";
    $maxMem =~ m|\D| and die "$error: max-mem must be an integer number of bytes\n";
    $minQual = $arrayFmt eq 'Illumina' ? 0.15 : 0;
}
# main execution block
sub msvtools_refine {

    # initialize
    setOptions_refine();
    print STDERR "$utility refine: " . getTime(), "\n";

    # load the percent GC for each probe
    print STDERR "loading probe fraction GC\n";
    my $gcFile ="$modelDir/msvtools.train.probes.$arrayType.gc.bed.gz";
    open my $gcH, "-|", "zcat $gcFile" or die "could not open: $gcFile\n";
    my $header = <$gcH>;
    my ($nGcProbes, $nKeptGcProbes) = (0, 0);
    while (my $probe = <$gcH>){
        $nGcProbes++;
        chomp $probe;
        my ($prbName, $gcVal) = split("\t", $probe);
        $gcVal or next;
        $gcVal eq "NA" and next;
        $nKeptGcProbes++;
        $fracGC{$prbName} = $gcVal;
    }
    close $gcH;
    print STDERR join("\t", commify($nGcProbes), "probes recovered from fraction GC file"), "\n";
    print STDERR join("\t", commify($nKeptGcProbes), "probes kept from fraction GC file"), "\n";
    print STDERR join("\t", commify(scalar(keys %fracGC)), "unique probes kept from fraction GC file"), "\n";

    # load the trained model CN (and read count) for each probe
    if($refineFile){ # this is a 2nd round refinement, use 1st round file for training
        print STDERR "loading previous model refinement\n";
        loadProbeModel();
        print STDERR join("\t", commify($nPrbUsed), "probes recovered from previous model refinement"), "\n";
        $ENV{SEGMENTING_ACTION} = "MODEL_REFINE_2";
    } else {
        print STDERR "loading training model\n";
        foreach my $list(['model', 'CN'],
                         ['bins',  'RC']){
            $inputDir = $modelDir;
            my $mdlFile = getInFile('train', $$list[0], $trainName, 'bed.bgz');
            -e $mdlFile or next; # bins file not present for bed6 or sexed training method
            print STDERR "using: $mdlFile\n";
            initializeExclude($mdlFile);
            $inputDir = $arrayDir;
            $sample = $samples[0];
            $valType = $$list[1];
            ($nPrbTotal, $nPrbTrained, $nPrbUsed) = (0, 0, 0);
            indexIllumina(\&getProbeTrainedVal);
            print STDERR join("\t", commify($nPrbTotal),   "probes recovered from indexing array"), "\n";
            print STDERR join("\t", commify($nPrbTrained), "index probes matched the training model"), "\n";
            print STDERR join("\t", commify($nPrbUsed),    "index probes probes had a reported $valType"), "\n";
        }
        $ENV{SEGMENTING_ACTION} = "MODEL_REFINE_1";
    }

    # load probe data for all reference samples
    print STDERR "loading individual sample arrays\n";
    $inputDir = $arrayDir;
    foreach my $sample_(@samples){
        $sample = $sample_;
        indexIllumina(\&getProbeData)
    }
    
    # save table of CN, LRR, BAF correlations
    print STDERR "calculating probe median values\n";
    my $datFile = getOutFile('data', $modelName, 'gz');
    openOutputStream($datFile, \$datH, $TRUE);
    print $datH join("\t", qw(
        CHROM START POS PRB_NAME IS_NA STRAND
        LRR BAF 
        CN_IN
        FRAC_GC
    )), "\n";
    my $nProbesWithData = 0;
    foreach my $prbName(keys %data){ # propagate all probes known to model
        my $d = $data{$prbName};
        $$d{LRR} or next;
        $$d{BAF} or next;
        $$d{CN}  or next;
        $fracGC{$prbName} or next;
        $nProbesWithData++;
        print $datH join("\t",
            @{$coord{$prbName}}, $prbName, -1, '+',
            (map { median(@{$$d{$_}}) } qw(LRR BAF CN)),
            $fracGC{$prbName}
        ), "\n";
    }
    closeHandles($datH);
    print STDERR join("\t", commify($nProbesWithData), "probes had CN, LRR, BAF and GC data to support model fitting"), "\n";

# ##########################
# $ENV{SEGMENTING_ACTION} = "MODEL_REFINE_1";
# my $datFile = getOutFile('data', $modelName, 'gz');

    # execute the HMM via R
    %data = ();
    $inputDir = $modelDir;
    my $rDataFile = $refineFile ? "__NULL__" : getInFile('train', 'model', $trainName, 'RData');
    $ENV{GC_BIAS_FILE} = getOutFile('gc_bias', $modelName, 'rds');
    segmentArray($modelName, $datFile, "$outputDir/plots", $rDataFile,
                 undef, undef, undef, $transProb);
}

# load a 'train' file to train a 1st round refinement
sub getProbeTrainedVal { 
    my ($f, $hdr) = @_;
    $nPrbTotal++;
    my @coord = ($$f[$$hdr{CHROM}], $$f[$$hdr{POS}]-1, $$f[$$hdr{POS}]);
    my $val = checkExclude(@coord);
    if(defined $val){
        $nPrbTrained++;  # only keep probes known to training model (e.g. those that could be sequenced)
        my $prbName = $$f[$$hdr{PRB_NAME}];
        $coord{$prbName} = \@coord;
        if($val ne 'NA'){
            $nPrbUsed++;
            $data{$prbName}{$valType} = [$val]; 
        }
    }
}

# load a 1st round 'refine.probes' file to train a 2nd round refinement
sub loadProbeModel {
    print STDERR "$utility refine: loading model: $refineFile\n";
    openInputStream($refineFile, \$mdlH, undef, undef, $TRUE);
    my $header = <$mdlH>;
    while (my $prb = <$mdlH>){
        chomp $prb;
        my @f = split("\t", $prb);
        my $prbName = $f[$prbCol{PRB_NAME}];
        $coord{$prbName} = [$f[$prbCol{CHROM}], $f[$prbCol{POS}]-1, $f[$prbCol{POS}]];
        $nPrbUsed++;
        $f[$prbCol{CN_OUT}] ne 'NA' and $data{$prbName}{CN} = [$f[$prbCol{CN_OUT}]];
    }
    closeHandles($mdlH); 
}

# load array data from one of a series of training samples
sub getProbeData {
    my ($f, $hdr) = @_;
    my $prbName = $$f[$$hdr{PRB_NAME}];
    $data{$prbName} or return;
    $$f[$$hdr{LRR}] ne 'NA' and push @{ $data{$prbName}{LRR} }, $$f[$$hdr{LRR}];
    $$f[$$hdr{BAF}] ne 'NA' and push @{ $data{$prbName}{BAF} }, $$f[$$hdr{BAF}];
}

1;
