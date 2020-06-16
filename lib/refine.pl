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
     $arrayDir, $arrayFmt, $project, $nameCol, $samples,
     $chroms, $outputDir, $plotDir, $maxCN);
our ($inputDir, $sample, $valType);
our (@samples, @chroms, %chroms, $minQual);
my ($mdlH, $datH);
my (%data, %coord);
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
    setOptionValue(\$samples,    'samples');
    setOptionValue(\$chroms,     'chromosomes');    
    setOptionValue(\$outputDir,  'output-dir');
    setOptionValue(\$maxCN,      'max-copy-number', 4);
    setOptionValue(\$tmpDir,     'tmp-dir',      '/tmp');
    setOptionValue(\$maxMem,     'max-mem',         1000000000);  
    #-------------------------------------
    ($modelDir or $refineFile) or die "$error: either model-dir, train-file or refine-file must be provided\n";    
    $modelDir   and (-d $modelDir   or die "$error: directory not found: $modelDir\n");
    $refineFile and (-e $refineFile or die "$error: file not found: $refineFile\n");
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
    
    # load the trained model CN (and read count) for each probe
    if($refineFile){ # this is a 2nd round refinement, use 1st round file for training
        loadProbeModel();
        print STDERR "$nPrbUsed probes extracted from previous model refinement\n";
    } else {
        foreach my $list(['model', 'CN'],
                         ['bins',  'RC']){
            $inputDir = $modelDir;
            my $mdlFile = getInFile('train', $$list[0], $trainName, 'bed.bgz');
            -e $mdlFile or next; # bins file not present for sexed training method
            print STDERR "using: $mdlFile\n";            
            initializeExclude($mdlFile);        
            $inputDir = $arrayDir;
            $sample = $samples[0];
            $valType = $$list[1];
            ($nPrbTotal, $nPrbTrained, $nPrbUsed) = (0, 0, 0);
            indexIllumina(\&getProbeTrainedVal);
            print STDERR "$nPrbTotal probes in indexing array\n";
            print STDERR "$nPrbTrained probes matched the training model\n";
            print STDERR "$nPrbUsed probes had a usable $valType\n";
        }        
    }
    
    # load probe data for all reference samples
    $inputDir = $arrayDir;    
    foreach my $sample_(@samples){
        $sample = $sample_;
        indexIllumina(\&getProbeData)
    }
    
    # save table of CN, LRR, BAF correlations
    my $datFile = getOutFile('data', $modelName, 'gz');
    openOutputStream($datFile, \$datH, $TRUE);
    print $datH join("\t", qw(
        CHROM START POS PRB_NAME IS_NA STRAND
        RC LRR BAF 
        CN_IN
    )), "\n";
    foreach my $probe(keys %coord){ # propagate all probes known to model
        my $d = $data{$probe};
        print $datH join("\t",
            @{$coord{$probe}}, $probe, -1, '+',
            map {
                defined $$d{$_} ? median(@{$$d{$_}}) : "NA";
            } qw(RC LRR BAF
                 CN),
        ), "\n";   
    
    }
    closeHandles($datH);

    # execute the HMM via R
    %data = ();
    $inputDir = $modelDir;
    my $rDataFile = $refineFile ? "__NULL__" : getInFile('train', 'model', $trainName, 'RData');
    segmentArray($modelName, $datFile, "$outputDir/plots", $rDataFile);
}

# load a 'train' file to train a 1st round refinement
sub getProbeTrainedVal { 
    my ($f, $hdr) = @_;
    $nPrbTotal++;
    my @coord = ($$f[$$hdr{CHROM}], $$f[$$hdr{POS}]-1, $$f[$$hdr{POS}]);
    my $val = checkExclude(@coord);
    if(defined $val){
        $nPrbTrained++;  # only keep probes known to training model (e.g. those that could be sequenced)
        $coord{$$f[$$hdr{PRB_NAME}]} = \@coord;
        if($val ne 'NA'){
            $nPrbUsed++;
            $data{$$f[$$hdr{PRB_NAME}]}{$valType} = [$val]; 
        }
    }
}

# load a 1st round 'refine.probes' file to train a 2nd round refinement
sub loadProbeModel {
    print STDERR "$utility refine: loading model: $refineFile\n";
    openInputStream($refineFile, \$mdlH, undef, undef, $TRUE);
    my $header = <$mdlH>;
    while (my $prb = <$mdlH>){
        $nPrbUsed++;
        chomp $prb;
        my @f = split("\t", $prb);
        $coord{$f[$prbCol{PRB_NAME}]} = [$f[$prbCol{CHROM}], $f[$prbCol{POS}]-1, $f[$prbCol{POS}]];
        $f[$prbCol{CN_OUT}] ne 'NA' and
            $data{$f[$prbCol{PRB_NAME}]}{CN} = [$f[$prbCol{CN_OUT}]];
        $f[$prbCol{RC}] ne 'NA' and
            $data{$f[$prbCol{PRB_NAME}]}{RC} = [$f[$prbCol{RC}]];
    }
    closeHandles($mdlH); 
}

# load array data from one of a series of training samples
sub getProbeData {
    my ($f, $hdr) = @_;
    $$f[$$hdr{LRR}] ne 'NA' and push @{ $data{$$f[$$hdr{PRB_NAME}]}{LRR} }, $$f[$$hdr{LRR}];
    $$f[$$hdr{BAF}] ne 'NA' and push @{ $data{$$f[$$hdr{PRB_NAME}]}{BAF} }, $$f[$$hdr{BAF}];
    #push @{ $data{$$f[$$hdr{PRB_NAME}]}{ZYG} }, abs($$f[$$hdr{BAF}] - 0.5) + 0.5;
}

1;
